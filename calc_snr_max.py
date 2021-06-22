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

font = {'size': 14}
matplotlib.rc('font', **font)

sys.path.append('./utils/')
from objects import load_object
from load_inputs import fill_data
from functions import *
# from kpf_etc.etc import kpf_photon_noise_estimate, kpf_etc_rv, kpf_etc_snr


from specutils import Spectrum1D
from specutils.manipulation import FluxConservingResampler

plt.ion()


def scale_to_kpf(so):
    """
    scale phoenix spectrum to be the right SNR across the KPF bandpass, save scaling factor to extend to Ca H&K

    inputs
    ------
    snr: float, signal to noise ratio of a resolution element in KPF

    outputs:
    --------
    kpf_spec: 1D array, phoenix spectrum scaled to input snr
    """
    # times exp time, telescope aperture, transmission
    s_ccd_hires = so.stel.s * so.var.exp_time * so.const.tel_area * so.kpf.ytransmit * np.abs(so.tel.s) ** 1.1547

    # convolve to lower res
    s_ccd_lores = degrade_spec(so.stel.v, s_ccd_hires, so.const.res_kpf)

    s_order_avg = np.zeros_like(so.kpf.order_wavelengths)
    for i, wl in enumerate(so.kpf.order_wavelengths):
        # find order subset
        fsr = 1.61e-5 * wl ** 2  # empirical stab at fsr, too lazy for constants
        order_sub = np.where(np.abs(so.stel.v - wl) < fsr / 2)[0]
        # resample for that order
        sig = wl / so.const.res_kpf / 3.5  # lambda/res = dlambda, 5 pixel sampling
        v_resamp, s_resamp = resample(so.stel.v[order_sub], s_ccd_lores[order_sub], sig=sig, dx=0, eta=1, mode='fast')

        # plt.plot(v_resamp,np.sqrt(s_resamp))
        # average the order and store
        s_order_avg[i] = np.sqrt(np.median(s_resamp))

    return so.kpf.order_wavelengths, s_order_avg


def scale_to_cahk(so):
    """
    scale output from kpf scaling to Ca H&K throughput
    """
    s_ccd_hires = so.stel.s * so.var.exp_time * so.const.tel_area * so.hk.ytransmit * np.abs(so.tel.s) ** 1.2

    # convolve to lower res
    s_ccd_lores = degrade_spec(so.stel.v, s_ccd_hires, so.const.res_hk)

    # resample onto res element grid
    sig = 395.0 / so.const.res_hk / 4  # lambda/res = dlambda
    v_resamp, s_resamp = resample(so.stel.v, s_ccd_lores, sig=sig, dx=0, eta=1, mode='fast')

    return v_resamp, np.sqrt(s_resamp)


def make_hk_grids(so):
    """
    loop through kpt_etc grids and save snr in same way

    saves to 'snr_master_hk.fits'

    """
    # Taken from Sam's KPF_ETC grids
    order_wls = np.array([401.29, 398.68, 396.11, 393.57, 391.07, 388.59])
    orders = np.array([153, 154, 155, 156, 157, 158])
    teffs = np.array([2700., 3000., 3300., 3600., 3900., 4200., 4500., 4800., 5100.,
                      5400., 5700., 6000., 6300., 6600.])
    vmags = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.,
                      14., 15., 16., 17., 18., 19.])
    exp_times = np.array([10., 50., 200., 500., 800., 1200., 2400., 3600.])

    # initiate snr_master_hk array
    snr_master_hk = np.zeros([len(orders), len(exp_times), len(vmags), len(teffs)])

    # Loop through
    for k, teff in enumerate(teffs):
        print(k)
        so.var.vmag = 0
        so.var.teff = teff
        cload.stellar(so)  # reload star
        for j, vmag in enumerate(vmags):
            so.stel.s *= 10 ** (-0.4 * vmag) / 10 ** (-0.4 * so.var.vmag)  # apply vmag scaling outside cload for speed
            so.var.vmag = vmag
            for i, exp_time in enumerate(exp_times):
                v, s = scale_to_cahk(so)  # full spectrum scaled to exp_time, vmag, teff
                for h, order in enumerate(orders):
                    isub = np.where(np.abs(v - order_wls[h]) < 1.23)[0]
                    snr_master_hk[h, i, j, k] = np.median(s[isub])

    # save to fits
    filename = 'snr_master_hk.fits'
    if os.path.isfile(filename): os.system('rm ' + filename)
    fits.writeto(filename, snr_master_hk)


if __name__ == '__main__':
    # load inputs
    configfile = 'calc_snr_max.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    # scale spectra
    vkpf, skpf = scale_to_kpf(so)
    vhk, shk = scale_to_cahk(so)
    # #sigma_rv_val, wvl_arr, snr_ord, dv_ord = kpf_photon_noise_estimate(so.var.teff, so.var.vmag, so.var.exp_time)

    plt.figure()
    plt.plot(vkpf, skpf, 'g', label='KPF')
    plt.plot(vhk, shk, 'orange', label='H&K')
    # plt.plot(wvl_arr, snr_ord,'k',label='KPF-ETC')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('SNR')
    plt.legend()
    plt.title('$t_{exp}$=%s  Teff=%s  vmag=%s' % (so.var.exp_time, so.var.teff, so.var.vmag))
    plt.savefig('fig1.png')

    # Make grids **to do: download phoenix files
    # make_hk_grids(so)

    # run for solar data, convert to per pixel instead of per pixel column
    so.var.vmag = 0
    so.var.teff = 5800
    so.var.exp_time = 1
    cload.stellar(so)  # reload star

    vkpf, skpf = scale_to_kpf(so)
    vhk, shk = scale_to_cahk(so)

    # plot pixel column version
    plt.figure()
    plt.plot(vkpf, skpf, 'g', label='KPF')
    plt.plot(vhk, shk, 'orange', label='H&K')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('SNR')
    plt.legend()
    plt.title('$t_{exp}$=%s  Teff=%s  vmag=%s' % (so.var.exp_time, so.var.teff, so.var.vmag))
    plt.savefig('fig2.png')

    # plot pixel version
    plt.figure()
    kpf_pixels = 45  # pixels in a column
    hk_pixels = 8  # pixels in a column
    plt.plot(vkpf, (0.8 / 0.65) * skpf ** 2 / kpf_pixels, 'g',
             label='KPF')  # 0.8/0.65 scales to actual fiber transmission
    plt.plot(vhk, (0.4 / 0.65) * shk ** 2 / hk_pixels, 'orange',
             label='H&K')  # 0.4/0.65 scales to actual fiber transmission
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Photons per pixel')
    plt.legend()
    plt.title('$t_{exp}$=%s  Teff=%s  vmag=%s' % (so.var.exp_time, so.var.teff, so.var.vmag))
    plt.savefig('fig3.png')
