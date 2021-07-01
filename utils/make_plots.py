import os

os.sys.path.append('./utils/')
from utils import functions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Add docstrings at some point this week

class MakePlots():
    def __init__(self, so):
        font = {'size': 8}
        matplotlib.rc('font', **font)
        self.info = so
        plt.ion()

    def plot_example_spectra(self, savefig=False):
        fig, axs = plt.subplots(2, 2, sharex=True)
        axs[0][0].plot(self.info.tel.v, self.info.tel.s, color='blue')
        axs[0][0].set(ylabel='Spectrum',
                      title="Telluric model")

        axs[0][1].plot(self.info.exo.v, self.info.exo.depth, color='maroon')
        axs[0][1].set(ylabel='Depth',
                      title="Exoplanet model")

        axs[1][0].plot(self.info.stel.vraw, self.info.stel.sraw, color='yellow')
        axs[1][0].set(ylabel='Spectrum',
                      title="Stellar model")

        axs[1][1].plot(self.info.hirax.wavegrid, self.info.hirax.hfp[0], color='black')
        axs[1][1].set(ylabel='Bandpass',
                      title="HIRAX model")
        fig.text(0.525, 0.02, 'wavelength, nm', ha='center')
        plt.tight_layout()
        if savefig:
            fig.savefig('./figures/example_spectra.png')

    def plot_no_doppler(self, savefig=False):
        xaxis = self.info.hirax.center_lam
        flux, flux_err = functions.calc_flux(self.info)
        no_exo_flux, no_exo_flux_err = functions.calc_flux(self.info, exoplanet=False)
        ratio = flux / no_exo_flux
        error_in_ratio = ratio * np.sqrt((flux_err / flux) ** 2 + (no_exo_flux_err / no_exo_flux) ** 2)
        fig, ax = plt.subplots()
        ax.errorbar(xaxis, ratio, yerr=error_in_ratio, fmt='.', color='black')
        ax.plot(xaxis, ratio, '-', alpha=0.3, color='blue')
        ax.set_xlabel("Center wavelengths, nm")
        ax.set_ylabel("Flux ratio")
        ax.set_title("Ratio of flux with and without exoplanet vs Center wavelengths")
        if savefig:
            fig.savefig('./figures/no_doppler_effect_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot.png', dpi=600)

    def plot_exo_doppler(self, savefig=False):
        xaxis = self.info.hirax.center_lam
        flux, flux_err = functions.calc_flux_exo_doppler(self.info)
        no_exo_flux, no_exo_flux_err = functions.calc_flux_exo_doppler(self.info, exoplanet=False)
        ratios = [fl / no_exo_fl for fl, no_exo_fl in zip(flux, no_exo_flux)]
        error_in_ratio = [r * np.sqrt((fr / f) ** 2 + (efr / ef) ** 2) for r, f, ef, fr, efr in
                          zip(ratios, flux, no_exo_flux, flux_err, no_exo_flux_err)]

        fig, ax = plt.subplots()
        labels = ['-100 km/s', '-50 km/s', '0 km/s', '50 km/s', '100 km/s']
        for ratio, yerr, label in zip(ratios, error_in_ratio, labels):
            ax.errorbar(xaxis, ratio, yerr=yerr, fmt='.', color="black", alpha=0.3)
            ax.plot(xaxis, ratio, '-', label=f"{label}")
        plt.xlabel("Center wavelengths, nm")
        plt.ylabel("Flux ratio")
        plt.title(
            "Ratio of flux with and without exoplanet vs Center wavelengths\n(including effects of doppler shift)")
        plt.legend()
        if savefig:
            fig.savefig('./figures/exo_doppler_effect_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot.png', dpi=600)
