##############################################################
# Makes neccessary plots
##############################################################
import os

from scipy.optimize import curve_fit
os.sys.path.append('./utils/')
from utils import functions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

class MakePlots():
    def __init__(self, so):
        font = {'size': 12}
        matplotlib.rc('font', **font)
        matplotlib.style.use('seaborn-pastel')
        self.info = so
        plt.ion()
        self.center_wavelengths = self.info.hirax.center_lam

    def plot_example_spectra(self, savefig=False):
        '''
        Plots the telluric, exoplanet, stellar and hirax models that are loaded into the storage object from
        the config file
        :param savefig: Saves the generated figure
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

    def plot_hirax_profile(self, savefig=False):
        fig, ax = plt.subplots()
        for profile in self.info.hirax.hfp:
            ax.fill(self.info.hirax.wavegrid, profile, color='black', alpha=0.3)
        ax.set_ylabel('Throughput')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_xlim((min(self.info.hirax.center_lam)-2, max(self.info.hirax.center_lam)+2))
        ax.set_ylim((0.4, 1.05))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xticks(ticks=np.arange(min(self.info.hirax.center_lam), max(self.info.hirax.center_lam)+2, 2))
        if savefig:
            fig.savefig(f"figures/hirax_profile_{self.info.hirax.hirax_file.split('_')[2].split('.')[0]}.png", dpi=600)

    def plot_no_doppler(self, savefig=False):
        '''
        Plots flux vs center wavelengths of the hirax bandpass when exoplanet and stellar speeds are 0 km/s
        :param savefig: Saves the generated figure
        '''
        # xaxis = self.info.hirax.center_lam
        flux, flux_err = functions.calc_flux(self.info)
        no_exo_flux, no_exo_flux_err = functions.calc_flux(self.info, exoplanet=False)
        ratio = flux / no_exo_flux
        error_in_ratio = ratio * np.sqrt((flux_err / flux) ** 2 + (no_exo_flux_err / no_exo_flux) ** 2)
        fig, ax = plt.subplots()
        ax.errorbar(self.center_wavelengths, ratio, xerr=self.info.hirax.width, yerr=error_in_ratio, fmt='.', color='black')
        ax.plot(self.center_wavelengths, ratio, '-', alpha=0.3, color='blue')
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel(r'$1-\left(\frac{R_{p}}{R_{s}}\right)^{2}$')
        ax.set_title(
            f"Ratio of flux with and without exoplanet vs Wavelength\nThroughput:{self.info.hirax.throughput[0]} VMag:{self.info.var.vmag}")
        ax.grid()
        if savefig:
            fig.savefig(
                f'figures/no_doppler_effect_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot_{self.info.const.l0}-{self.info.const.l1}_tput-{self.info.hirax.throughput[0]}.png',
                dpi=600)

    def plot_exo_doppler(self, savefig=False):
        '''
        Plots flux vs center wavelengths of the hirax bandpass when exoplanet speeds are nonzero and stellar speeds are
        0 km/s i.e. it plots the flux due to doppler shift caused by exoplanet's motion.
        :param savefig: Saves the generated figure
        '''
        # xaxis = self.info.hirax.center_lam
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
            ax.errorbar(self.center_wavelengths, ratio, xerr=self.info.hirax.width, yerr=yerr, fmt='.', color="black", alpha=0.3)
            ax.plot(self.center_wavelengths, ratio, '-', label=f"{label}")
            ax.grid()
        plt.xlabel("Wavelength (nm)")
        plt.ylabel(r'$1-\left(\frac{R_{p}}{R_{s}}\right)^{2}$')
        plt.title(
            f"Flux Ratio vs Wavelength with effect of doppler shift\nThroughput:{self.info.hirax.throughput[0]} VMag:{self.info.var.vmag}")
        plt.legend()
        if savefig:
            fig.savefig(
                f'figures/exo_doppler_effect_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot_{self.info.const.l0}-{self.info.const.l1}.png',
                dpi=600)

    ########## THIS CODE FOR ONE STEL SPEED AT A TIME ################
    def plot_stel_exo_doppler_one(self, savefig=False):
        '''
        Plots flux vs center wavelengths of the hirax bandpass when exoplanet speeds are nonzero and stellar speeds are
        nonzero i.e. it plots the flux due to doppler shift caused by exoplanet's AND star's motion
        :param savefig: Saves the generated figure
        '''
        # xaxis = self.info.hirax.center_lam
        stel_speed = self.info.stel.speed
        exo_speeds = self.info.exo.speed + stel_speed
        flux = [functions.calc_flux(self.info, exo_speed=speed, stel_speed=stel_speed) for speed in exo_speeds]
        no_exo_flux = [functions.calc_flux(self.info, exoplanet=False, exo_speed=speed, stel_speed=stel_speed) for speed
                       in exo_speeds]
        ratios = [fl[0] / no_exo_fl[0] for fl, no_exo_fl in zip(flux, no_exo_flux)]
        error_in_ratio = [r * np.sqrt((f[1] / f[0]) ** 2 + (ef[1] / ef[0]) ** 2) for r, f, ef in
                          zip(ratios, flux, no_exo_flux)]

        fig, ax = plt.subplots()
        labels = [f"{str(speed / 1000)} km/s" for speed in exo_speeds]
        for ratio, yerr, label in zip(ratios, error_in_ratio, labels):
            ax.errorbar(self.center_wavelengths, ratio, yerr=yerr, fmt='.', color="black", alpha=0.3)
            ax.plot(self.center_wavelengths, ratio, '-', label=f"{label}")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel(r'$1-\left(\frac{R_{p}}{R_{s}}\right)^{2}$')
        plt.title(
            f"Ratio of flux with and without exoplanet vs Wavelength\nStellar speed = {stel_speed / 1000} km/s")
        plt.legend()
        if savefig:
            plt.savefig(
                f'./figures/ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot-exospeeds_{exo_speeds[0] / 100}-{exo_speeds[-1] / 1000}kms_stelspeed_{stel_speed / 1000}kms.png',
                dpi=600)

    # need to change this a little to make it look nicer
    def plot_stel_exo_doppler(self, savefig=False):
        '''
        Plots n x 1 flux vs center wavelengths of the hirax bandpass when exoplanet speeds are nonzero and stellar speeds
        are nonzero i.e. it plots the flux due to doppler shift caused by exoplanet's AND star's motion. In this function
        the stellar speed is in the form of an array.
        :param savefig: Saves the generated figure
        '''
        # xaxis = self.info.hirax.center_lam
        stel_speeds = self.info.stel.speed

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
                ax[i].errorbar(self.center_wavelengths, ratio, xerr=self.info.hirax.width, yerr=yerr, fmt='.', color="black", alpha=0.3)
                ax[i].plot(self.center_wavelengths, ratio, '-', label=f"{label}")
                ax[i].legend(bbox_to_anchor=(1.01, 1.0), loc='upper left')
                ax[i].grid()
        fig.suptitle(f"Flux Ratio vs Wavelength with effect of doppler shift\nThroughput:{self.info.hirax.throughput[0]} VMag:{self.info.var.vmag}")
        fig.supxlabel('Wavelength (nm)')
        fig.supylabel(r'$1-\left(\frac{R_{p}}{R_{s}}\right)^{2}$')
        if savefig:
            plt.savefig(
                f'figures/exo_stel_doppler_ratio_of_fluxes_w_and_wo_exoplanet_errorbar_plot.png',
                dpi=600)

    def plot_hirax_over_spectra(self, savefig=False):
        '''
        Plots the hirax profiles over the telluric, exoplanet, and OH spectra, with the exoplnaet spectra on a twinx
        axis due to it very different y-scale
        :param savefig: Saves the generated figure
        '''
        fig, ax1 = plt.subplots()

        lns1 = ax1.plot(self.info.xgrid, self.info.tel.s, color='blue', label="Telluric")

        for ind, profile in enumerate(self.info.hirax.hfp):
            if ind == 0:
                lns2 = ax1.plot(self.info.xgrid, profile, color='black', alpha=0.3, label="HIRAX")
                ax1.fill_between(self.info.xgrid, profile, 0, alpha=0.4,color='gray')
            else:
                ax1.plot(self.info.xgrid, profile, color='black', alpha=0.3)
                ax1.fill_between(self.info.xgrid, profile, 0, alpha=0.4,color='gray')
        lns3 = ax1.plot(self.info.xgrid, (self.info.oh.s) ** 100, color="green", label="OH")

        ax2 = ax1.twinx()
        lns4 = ax2.plot(self.info.xgrid, self.info.exo.depth, color="brown", label="exoplanet")
        
        ratio, error_in_ratio = functions.calc_flux_ratio(self.info)
        ax2.errorbar(self.info.hirax.center_lam, ratio, xerr=self.info.hirax.width,\
            yerr=error_in_ratio, fmt='o',\
            label="Simulated Data",color='r')

        lns = lns1 + lns2 + lns3 + lns4
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc=0)

        ax1.set(xlabel="Wavelength (nm)", ylabel='Hirax throughput',
                title=f"No. Transits: {self.info.var.num_transits}"
                      f"\nThroughput:{self.info.hirax.throughput[0]} VMag:{self.info.var.vmag}")
        ax2.set(ylabel=r'$1-\left(\frac{R_{p}}{R_{s}}\right)^{2}$')
        plt.tight_layout()
        if savefig:
            fig.savefig(
                f'figures/hirax_telluric_exoplanet_OH.png', dpi=600)

    def plot_hirax_over_shifted_exo(self, savefig=False):
        '''
        Plots the hirax profiles over the shifted exoplanet models due to doppler effect with the hirax profile and
        the exoplanet spectra being on different y-scales
        :param savefig: Saves the generated figure
        '''
        exo_speeds = np.arange(-200, 300, 100) * 10 ** 3
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), sharey=True, constrained_layout=True)
        labels = [f"{str(speed / 1000)} km/s" for speed in exo_speeds]
        for i in range(len(axs)):
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)
            for ind, speed in enumerate(exo_speeds):
                shifted_exo_spec = functions.shift_spectra(self.info, xgrid=self.info.xgrid, speed=speed, choice = 'exoplanet')
                non_zero_shifted_exo_spec = shifted_exo_spec[np.nonzero(shifted_exo_spec)]
                x = self.info.exo.v[np.nonzero(shifted_exo_spec)]
                axs[i].plot(x, non_zero_shifted_exo_spec, label=labels[ind])
                axs[i].set_xlim((580, 600))
                axs[i].grid()
                if i == 1:
                    ax2.set_ylabel("Throughput")
                    axs[i].set_xlim((587, 591))
                    axs[i].legend(bbox_to_anchor=(1.01, 1.0), loc='upper left')

            ax2 = axs[i].twinx()
            for profile in self.info.hirax.hfp:
                ax2.fill(self.info.hirax.wavegrid, profile, color='black', alpha=0.3)

        ax2.legend(['Hirax'], loc='upper left')
        fig.suptitle(f'Exoplanet spectra vs Wavelength'
                     f'\nThroughput:{self.info.hirax.throughput[0]}')
        fig.supxlabel('Wavelengths, nm')
        fig.supylabel(r'$1-\left(\frac{R_{p}}{R_{s}}\right)^{2}$')
        if savefig:
            fig.savefig(f"figures/hirax_over_shifted_exo.png", dpi=600)

    def plot_lorentz_fit(self, savefig=False, feature = 'sodium'):
        '''
        :param savefig - saves the generated figure if True
        :param feature - takes a string input of what kind of feature is being targeted
        '''

        def lorentz_sodium(x, a, w):
            '''
            :param a: amplitude
            :param w: width
            :param c: center
            :param k: constant
            :return: function output
            '''
            c = 5.89610638 * 10 ** 2
            k = 9.80597083 * 10 ** (-1)
            return ((a * w) / (w ** 2 + (x - c) ** 2)) + k

        def lorentz_potassium(x, a, w):
            '''
            :param a: amplitude
            :param w: width
            :param c: center
            :param k: constant
            :return: function output
            '''
            c = 7.66907675 * 10 ** 2
            k = 9.80710064 * 10 ** (-1)
            return ((a * w) / (w ** 2 + (x - c) ** 2)) + k

        def fit_lorentz_aw(x_data, y_data, y_err, init_guess, choice=feature):
            '''
            :param x_data: x data to be fitted
            :param y_data: y data to be fitted
            :param y_err: error in y data to be fitted
            :param init_guess: initial guess as a list
            :return: fitted parameters in the initial guess
            '''
            if choice == 'sodium':
                popt, pcov = curve_fit(lorentz_sodium, x_data, y_data, sigma=y_err, p0=init_guess,
                                       absolute_sigma=True)
            elif choice == 'potassium':
                popt, pcov = curve_fit(lorentz_potassium, x_data, y_data, sigma=y_err, p0=init_guess,
                                       absolute_sigma=True)
            return -popt, pcov

        ratio, error_in_ratio = functions.calc_flux_ratio(self.info)
        popt, pcov = fit_lorentz_aw(self.info.hirax.center_lam, ratio, y_err=error_in_ratio, init_guess=[-0.00035, 3], choice=feature)

        # plotting
        fig, ax = plt.subplots()
        ax.errorbar(self.info.hirax.center_lam, ratio, xerr=self.info.hirax.width, yerr=error_in_ratio, fmt='-', label="Simulated Data")
        # if feature == 'sodium':
        #     ax.plot(self.info.xgrid, lorentz_sodium(self.info.xgrid, *popt), '-', label='fit')
        #     ax.fill_between(self.info.xgrid, lorentz_sodium(self.info.xgrid, popt[0] + np.sqrt(pcov[0][0]), popt[1]),
        #                     lorentz_sodium(self.info.xgrid, popt[0] - np.sqrt(pcov[0][0]), popt[1]), alpha=0.2)
        #     ax.plot(self.info.hirax.center_lam, lorentz_sodium(self.info.hirax.center_lam, *popt), '.')
        #     ax.set_xlim((575, 600))
        # if feature == 'potassium':
        #     ax.plot(self.info.xgrid, lorentz_potassium(self.info.xgrid, *popt), '-', label='fit')
        #     ax.fill_between(self.info.xgrid, lorentz_potassium(self.info.xgrid, popt[0] + np.sqrt(pcov[0][0]), popt[1]),
        #                     lorentz_potassium(self.info.xgrid, popt[0] - np.sqrt(pcov[0][0]), popt[1]), alpha=0.2)
        #     ax.plot(self.info.hirax.center_lam, lorentz_potassium(self.info.hirax.center_lam, *popt), '.')
        #     ax.set_xlim((750, 780))

        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Ratio of flux')
        #ax.set_title(f'Fitting the data to an inverse lorentz profile with center and constant fixed'
        #             f'\nThroughput:{self.info.hirax.throughput[0]} VMag:{self.info.var.vmag}')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
        if savefig:
            fig.savefig(f'figures/lorentz_fit_{feature}', dpi=600)

    def plot_noise_mag(self, cload, magnitudes, savefig=False):
        '''
        :param cload - loads the new info into the code
        :param magnitudes - an array of magnitudes used as the independent variable in the plot
        :param savefig - saves the generated figure if True
        '''
        photon_noise, read_noise, dark_noise = [], [], []
        default = self.info.var.vmag
        for mag in magnitudes:
            self.info.var.vmag = mag
            cload.stellar(self.info)
            cload.hirax(self.info, tfac=1)
            flux, _ = functions.calc_flux(self.info)  # only one band
            w_prime = (self.info.const.theta_s * self.info.const.focal_length * self.info.var.magnification) / 206265
            psf = w_prime / (self.info.const.pixel_size * 10 ** -6)
            n_pix = np.pi * psf ** 2
            ph_noise = np.sqrt(flux[0])
            rd_noise = self.info.var.read_n * np.sqrt((flux[0] * n_pix / (self.info.var.saturation * n_pix)))
            drk_noise = np.sqrt(self.info.var.dark_n * self.info.var.transit_duration * n_pix)
            photon_noise.append(ph_noise)
            read_noise.append(rd_noise)
            dark_noise.append(drk_noise)
        fig, ax = plt.subplots()
        ax.plot(magnitudes, photon_noise, '-', label='photon noise')
        ax.plot(magnitudes, read_noise, '-', label='read noise')
        ax.plot(magnitudes, dark_noise, '-', label='dark noise')
        ax.set_xlabel('v mag')
        ax.set_ylabel('noise')
        ax.set_title(f'Noise vs Magnitude\nThroughput:{self.info.hirax.throughput[0]}')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(loc='upper right')
        ax.legend()
        self.info.var.vmag = default
        if savefig:
            fig.savefig(f"figures/noise_mag_plot.png", dpi=600)
        np.savetxt(f"data/output/noise_vmag_{self.info.hirax.hirax_file.split('_')[2].split('.')[0]}.txt",
                   np.transpose([magnitudes, photon_noise, read_noise, dark_noise]), fmt='%1.4e', delimiter=',',
                   header='magnitude, photon_noise, read_noise, dark_noise')

    def plot_SNR_mag_tput(self, magnitudes, snr, tfac = np.arange(0.2 , 1, 0.2), savefig=False):
        '''
        :param magnitudes - an array of magnitudes used as the independent variable in the plot
        :param amplitudes - an array of SNR used as the dependent variable in the plot
        :param tfac - an array used to change the throughput of the hirax profile
        :param savefig - saves the generated figure if True
        '''

        fig, ax = plt.subplots()
        for i in range(len(tfac)):
            ax.plot(magnitudes, np.array(snr[i]), '-', label = f'tput_factor = {round(tfac[i], 3)}')
        ax.set_xlabel('V Magnitude')
        ax.set_ylabel('SNR')
        ax.set_title(f"SNR vs Magnitude - {self.info.hirax.hirax_file.split('_')[2].split('.')[0]}")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(loc='upper right')
        if savefig:
            fig.savefig(f"figures/SNR_mag_tput_plot_{self.info.hirax.hirax_file.split('_')[2].split('.')[0]}", dpi=600)

    def plot_vmag_signal_percent_noise_tfac(self, vmag_signal_data_file_name=\
        'data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv', \
        config='config1', tfactor=np.arange(0.2 , 1, 0.2), savefig = False,\
        sigma_thresh=3):
            """
            """
            vmag_signal_data = pd.read_csv(vmag_signal_data_file_name)
            vmag = vmag_signal_data['sy_vmag']
            signal_percent = vmag_signal_data['signal_percent']
            labels = vmag_signal_data['pl_name']
            path = 'data/output/'
            file = 'amplitude_vmag_'
            extension = '.csv'
            noise = pd.read_csv(path + file + config + extension)
            fig, ax = plt.subplots(figsize=(12, 10))
            ax.semilogx(signal_percent, vmag, '+k')
            ax.set_ylim((16, 7))
            ax.set_xlim(0.01, 0.1)
            ax.set_xlabel('Transit signal of 2 Atmospheric Scale Height (%)')
            ax.set_ylabel('V magnitude')
            ax.set_xticks([0.01, 0.1], minor=True)

            for i, txt in enumerate(labels):
                ax.annotate(txt, (signal_percent[i], vmag[i]), fontsize=5)

            for j,tfac in enumerate(tfactor):
                ax.semilogx(100 * (sigma_thresh)\
                    * noise[f"amplitude_err_tfac{j + 1}"],noise['magnitude'],\
                    '-',label=f"tfac = {round(tfac, 3)}", alpha=0.8) 

            fig.suptitle(file.split('_')[2])
            if savefig:
                fig.savefig(f"figures/vmag_signal_percent_{config}.png", dpi=600)


    def plot_vmag_signal_percent_noise_tfac_ashley(self, vmag_signal_data_file_name=\
        'data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv', \
        config='config1', tfactor=np.arange(0.2 , 1, 0.2), savefig = False,\
        sigma_thresh=3):
            """
            for palomar proposal - highlight observable planets with red star, compare to smaller aperture telescope
            assume 40%
            """
            def load_noise_data(config):
                path = 'data/output/'
                file = 'amplitude_vmag_'
                extension = '.csv'
                return pd.read_csv(path + file + config + extension)

            vmag_signal_data = pd.read_csv(vmag_signal_data_file_name)
            vmag = vmag_signal_data['sy_vmag']
            signal_percent = vmag_signal_data['signal_percent']* 2.5/2
            labels = vmag_signal_data['pl_name']
            noise = load_noise_data(config)
            noise2 = load_noise_data(config + '_60m_3transits')

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.semilogx(signal_percent , vmag, '+k')
            ax.set_ylim((12, 7))
            ax.set_xlim(0.01, 0.1)
            ax.set_xlabel('Differential Transit Signal (%)')
            ax.set_ylabel('V magnitude')
            ax.set_xticks([0.01, 0.1], minor=True)

            for i, txt in enumerate(labels):
                if txt in ['HD 209458 b', 'WASP-69 b', 'KELT-9 b', 'MASCARA-2 b/KELT-20 b']:
                    ax.plot(signal_percent[i]-0.0001, vmag[i]+.05, '*r',ms=10)
                    ax.annotate(txt, (signal_percent[i], vmag[i]), fontsize=9)
                else:
                    ax.annotate(txt, (signal_percent[i], vmag[i]), fontsize=5)

            tfacs = np.arange(0.2 , 1, 0.2)
            itfac = 1 # just plot 0.4 throughput line
            ax.semilogx(100 * sigma_thresh\
                    * noise[f"amplitude_err_tfac{itfac + 1}"],noise['magnitude'],\
                    'goldenrod',label=f"tfac = {round(tfacs[itfac], 3)}", alpha=0.8) 
            ax.semilogx(100 * sigma_thresh/np.sqrt(2)\
                    * noise[f"amplitude_err_tfac{itfac + 1}"],noise['magnitude'],\
                    'goldenrod',ls='--',label=f"tfac = {round(tfacs[itfac], 3)}", alpha=0.8) 
            ax.semilogx(100 * sigma_thresh\
                    * noise2[f"amplitude_err_tfac{itfac + 1}"],noise2['magnitude'],\
                    'c-',label=f"tfac = {round(tfacs[itfac], 3)}", alpha=0.8) 
            itfac = 3 # 0.8 throughput
            ax.semilogx(100 * sigma_thresh/np.sqrt(3)\
                    * noise[f"amplitude_err_tfac{itfac + 1}"],noise['magnitude'],\
                    'm',ls='-.',label=f"tfac = {round(tfacs[itfac], 3)}", alpha=0.8) 

            if savefig:
                plt.savefig(f"figures/vmag_signal_percent_{config}_ashley.png", dpi=600)





