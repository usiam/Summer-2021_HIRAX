##############################################################
# Makes neccessary plots
##############################################################
import os
os.sys.path.append('./utils/')
from utils import functions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


class MakePlots():
    def __init__(self, so):
        # at class initialization, the storage_object is called and the matplotlib necessary functions are called
        font = {'size': 8}
        matplotlib.rc('font', **font)
        self.info = so
        plt.ion()

    def plot_example_spectra(self, savefig=False):
        '''

        '''
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
            fig.savefig('figures/example_spectra.png')

    def plot_no_doppler(self, savefig=False):
        '''
        :param savefig:
        :return:
        '''
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
        ax.set_title(
            f"Ratio of flux with and without exoplanet vs Center wavelengths\n(Throughput:{self.info.hirax.throughput[0]})")
        ax.grid()
        if savefig:
            fig.savefig(
                f'figures/no_doppler_effect_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot_{self.info.const.l0}-{self.info.const.l1}_tput-{self.info.hirax.throughput[0]}.png',
                dpi=600)

    def plot_exo_doppler(self, savefig=False):
        '''

        :param savefig:
        :return:
        '''
        xaxis = self.info.hirax.center_lam
        exo_speeds = self.info.exo.speed
        flux = [functions.calc_flux(self.info, exo_speed=speed, stel_speed=0) for speed in exo_speeds]
        no_exo_flux = [functions.calc_flux(self.info, exoplanet=False, exo_speed=speed, stel_speed=0) for
                       speed in exo_speeds]
        ratios = [fl[0] / no_exo_fl[0] for fl, no_exo_fl in zip(flux, no_exo_flux)]
        error_in_ratio = [r * np.sqrt((f[1] / f[0]) ** 2 + (ef[1] / ef[0]) ** 2) for r, f, ef in
                          zip(ratios, flux, no_exo_flux)]
        labels = [f"{str(speed / 1000)} km/s" for speed in exo_speeds]
        fig, ax = plt.subplots()
        for ratio, yerr, label in zip(ratios, error_in_ratio, labels):
            ax.errorbar(xaxis, ratio, yerr=yerr, fmt='.', color="black", alpha=0.3)
            ax.plot(xaxis, ratio, '-', label=f"{label}")
            ax.grid()
        plt.xlabel("Center wavelengths, nm")
        plt.ylabel("Flux ratio")
        plt.title(
            "Ratio of flux with and without exoplanet vs Center wavelengths\n(including effects of doppler shift)")
        plt.legend()
        if savefig:
            fig.savefig(
                f'figures/exo_doppler_effect_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot_{self.info.const.l0}-{self.info.const.l1}.png',
                dpi=600)

    ########## THIS CODE FOR ONE STEL SPEED AT A TIME ################
    def plot_stel_exo_doppler_one(self, savefig=False):
        '''

        :param savefig:
        :return:
        '''
        xaxis = self.info.hirax.center_lam
        stel_speed = self.info.stel.speed
        exo_speeds = self.info.exo.speed + stel_speed
        flux = [functions.calc_flux(self.info, exo_speed=speed, stel_speed=stel_speed) for speed in exo_speeds]
        no_exo_flux= [functions.calc_flux(self.info, exoplanet=False, exo_speed=speed, stel_speed=stel_speed) for speed in exo_speeds]
        ratios = [fl[0] / no_exo_fl[0] for fl, no_exo_fl in zip(flux, no_exo_flux)]
        error_in_ratio = [r * np.sqrt((f[1] / f[0]) ** 2 + (ef[1] / ef[0]) ** 2) for r, f, ef in
                          zip(ratios, flux, no_exo_flux)]

        fig, ax = plt.subplots()
        labels = [f"{str(speed/1000)} km/s" for speed in exo_speeds]
        for ratio, yerr, label in zip(ratios, error_in_ratio, labels):
            ax.errorbar(xaxis, ratio, yerr=yerr, fmt='.', color="black", alpha=0.3)
            ax.plot(xaxis, ratio, '-', label=f"{label}")
        plt.xlabel("Center wavelengths, nm")
        plt.ylabel("Flux ratio")
        plt.title(
            f"Ratio of flux with and without exoplanet vs Center wavelengths\nStellar speed = {stel_speed/1000} km/s")
        plt.legend()
        plt.savefig(f'./figures/ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot-exospeeds_{exo_speeds[0]/100}-{exo_speeds[-1]/1000}kms_stelspeed_{stel_speed/1000}kms.png', dpi=600)

    def plot_stel_exo_doppler(self, savefig=False):
        '''

        :param savefig:
        :return:
        '''
        xaxis = self.info.hirax.center_lam
        stel_speeds = np.arange(-100, 150, 50) * 10 ** 3

        fig, ax = plt.subplots(len(stel_speeds), 1, sharex=True, constrained_layout=True, figsize=(6, 20))
        for stel_speed in stel_speeds:
            exo_speeds = self.info.exo.speed + stel_speed
            flux = [functions.calc_flux(self.info, exo_speed=speed, stel_speed=stel_speed) for speed in exo_speeds]
            no_exo_flux = [functions.calc_flux(self.info, exoplanet=False, exo_speed=speed, stel_speed=stel_speed) for
                           speed in exo_speeds]
            ratios = [fl[0] / no_exo_fl[0] for fl, no_exo_fl in zip(flux, no_exo_flux)]
            error_in_ratio = [r * np.sqrt((f[1] / f[0]) ** 2 + (ef[1] / ef[0]) ** 2) for r, f, ef in
                              zip(ratios, flux, no_exo_flux)]

            labels = [f"{str(speed / 1000)} km/s" for speed in exo_speeds]
            num = np.arange(0, len(ratios), 1)
            for ratio, yerr, label, i in zip(ratios, error_in_ratio, labels, num):
                ax[i].errorbar(xaxis, ratio, yerr=yerr, fmt='.', color="black", alpha=0.3)
                ax[i].plot(xaxis, ratio, '-', label=f"{label}")
                ax[i].legend(bbox_to_anchor=(1.01, 1.0), loc='upper left')
                ax[i].grid()
        fig.suptitle(f'Ratio of flux with and without exoplanet vs Center wavelengths')
        fig.supxlabel('Center wavelengths, nm')
        fig.supylabel('Flux ratio')
        if savefig:
            plt.savefig(
                f'figures/exo_stel_doppler_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot.png',
                dpi=600)

    def plot_hirax_over_spectra(self, savefig=False):
        '''

        :param savefig:
        :return:
        '''
        fig, ax1 = plt.subplots()

        lns1 = ax1.plot(self.info.xgrid, self.info.tel.s, color='blue', label="telluric")

        for ind, profile in enumerate(self.info.hirax.hfp):
            if ind == 0:
                lns2 = ax1.plot(self.info.xgrid, profile, color='black', alpha=0.3, label="hirax")
            else:
                ax1.plot(self.info.xgrid, profile, color='black', alpha=0.3)
        lns3 = ax1.plot(self.info.xgrid, (self.info.oh.s) ** 100, color = "green", label = "OH")

        ax2 = ax1.twinx()
        lns4 = ax2.plot(self.info.xgrid, self.info.exo.depth, color="brown", label = "exoplanet")

        lns = lns1 + lns2 + lns3 + lns4
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc=0)

        ax1.set(xlabel="Center wavelengths, nm", ylabel='Hirax throughput',
                title=f"Telluric, Exoplanet, OH lines, and Hirax model\nCenter wavelengths: {self.info.hirax.center_lam}")
        ax2.set(ylabel="Depth")
        plt.tight_layout()
        if savefig:
            fig.savefig(
                f'figures/hirax_telluric_exoplanet_OH.png', dpi=600)
