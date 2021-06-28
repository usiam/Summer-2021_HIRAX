# calc signal to noise
# max for the calcium H&K
import sys
import matplotlib
import os

import numpy as np
import matplotlib.pylab as plt

from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.table import Table

from utils import functions

font = {'size': 8}
matplotlib.rc('font', **font)

os.sys.path.append('./utils/')
from objects import load_object
from load_inputs import fill_data
from functions import *

# from specutils import Spectrum1D
# from specutils.manipulation import FluxConservingResampler

plt.ion()

if __name__ == '__main__':
    # load inputs
    configfile = 'calc_snr_max.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    flux = functions.calc_flux(so)
    no_exo_flux = functions.calc_flux(so, exoplanet=False)
   
    plt.plot(so.hirax.center_lam, flux, '.', label="With exoplanet")
    plt.plot(so.hirax.center_lam, no_exo_flux, '.', label="No exoplanet")
    plt.savefig('./figures/no_exoplanet_plot.png')
    plt.xlabel("Center wavelengths, nm")
    plt.ylabel("Flux/sec")
    plt.title("Flux/sec vs Center wavelengths with exoplanet and without exoplanet")
    plt.legend()

    # fig, axs = plt.subplots(2, 2, sharex=True)
    # axs[0][0].plot(so.tel.v, so.tel.s, color='blue')
    # axs[0][0].set(ylabel='Spectrum',
    #               title="Telluric model")
    #
    # axs[0][1].plot(so.exo.v, so.exo.depth, color='maroon')
    # axs[0][1].set(ylabel='Depth',
    #               title="Exoplanet model")
    #
    # axs[1][0].plot(so.stel.vraw, so.stel.sraw, color='yellow')
    # axs[1][0].set(ylabel='Spectrum',
    #               title="Stellar model")
    #
    # axs[1][1].plot(so.hirax.wavegrid, so.hirax.hfp[0], color='black')
    # axs[1][1].set(ylabel='Bandpass',
    #               title="HIRAX model")
    # fig.text(0.525, 0.02, 'wavelength, nm', ha='center')
    # plt.tight_layout()
    # plt.savefig('./figures/presentation.png')
    
    # plt.plot(so.hirax.wavegrid, so.hirax.hfp[0], color='black')
    # plt.xlabel("wavelength, nm")
    # plt.ylabel("Bandpass")
    # plt.title('Gaussian filter profile with FWHM = 0.3 nm, center = 589.3 nm and throughput = 0.95')
    # plt.xlim([585, 593])
    # plt.savefig('./figures/hirax.png')


    # # scale spectra
    # vkpf, skpf = scale_to_kpf(so)
    # vhk, shk = scale_to_cahk(so)
    # # #sigma_rv_val, wvl_arr, snr_ord, dv_ord = kpf_photon_noise_estimate(so.var.teff, so.var.vmag, so.var.exp_time)
    #
    # plt.figure()
    # plt.plot(vkpf, skpf, 'g', label='KPF')
    # plt.plot(vhk, shk, 'orange', label='H&K')
    # # plt.plot(wvl_arr, snr_ord,'k',label='KPF-ETC')
    # plt.xlabel('Wavelength (nm)')
    # plt.ylabel('SNR')
    # plt.legend()
    # plt.title('$t_{exp}$=%s  Teff=%s  vmag=%s' % (so.var.exp_time, so.var.teff, so.var.vmag))
    # plt.savefig('fig1.png')
    #
    # # Make grids **to do: download phoenix files
    # # make_hk_grids(so)
    #
    # # run for solar data, convert to per pixel instead of per pixel column
    # so.var.vmag = 0
    # so.var.teff = 5800
    # so.var.exp_time = 1
    # cload.stellar(so)  # reload star
    #
    # vkpf, skpf = scale_to_kpf(so)
    # vhk, shk = scale_to_cahk(so)
    #
    # # plot pixel column version
    # plt.figure()
    # plt.plot(vkpf, skpf, 'g', label='KPF')
    # plt.plot(vhk, shk, 'orange', label='H&K')
    # plt.xlabel('Wavelength (nm)')
    # plt.ylabel('SNR')
    # plt.legend()
    # plt.title('$t_{exp}$=%s  Teff=%s  vmag=%s' % (so.var.exp_time, so.var.teff, so.var.vmag))
    # plt.savefig('fig2.png')
    #
    # # plot pixel version
    # plt.figure()
    # kpf_pixels = 45  # pixels in a column
    # hk_pixels = 8  # pixels in a column
    # plt.plot(vkpf, (0.8 / 0.65) * skpf ** 2 / kpf_pixels, 'g',
    #          label='KPF')  # 0.8/0.65 scales to actual fiber transmission
    # plt.plot(vhk, (0.4 / 0.65) * shk ** 2 / hk_pixels, 'orange',
    #          label='H&K')  # 0.4/0.65 scales to actual fiber transmission
    # plt.xlabel('Wavelength (nm)')
    # plt.ylabel('Photons per pixel')
    # plt.legend()
    # plt.title('$t_{exp}$=%s  Teff=%s  vmag=%s' % (so.var.exp_time, so.var.teff, so.var.vmag))
    # plt.savefig('fig3.png')
