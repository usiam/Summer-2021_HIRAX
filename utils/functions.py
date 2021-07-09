##############################################################
# General functions for calc_snr_max
##############################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from scipy import signal
from scipy import interpolate, constants
from utils.objects import storage_object

all = {'integrate', 'gaussian', 'define_lsf', 'vac_to_stand', 'setup_band', 'resample'}


def integrate(x, y):
    """
    Integrate y wrt x
    """
    return trapz(y, x=x)


def gaussian(x, shift, sig):
    ' Return normalized gaussian with mean shift and var = sig^2 '
    return np.exp(-.5 * ((x - shift) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))


def define_lsf(v, res):
    """
    define gaussian in pixel elements to convolve resolved spectrum with to get rightish resolution
    """
    dlam = np.median(v) / res
    fwhm = dlam / np.mean(
        np.diff(v))  # desired lambda spacing over current lambda spacing resolved to give sigma in array elements
    sigma = fwhm / 2.634  # FWHM is dl/l but feed sigma
    x = np.arange(sigma * 10)
    gaussian = (1. / sigma / np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - 0.5 * len(x)) / sigma) ** 2)

    return gaussian


def degrade_spec(x, y, res):
    """
    given wavelength, flux array, and resolving power R, return  spectrum at that R
    """
    lsf = define_lsf(x, res=res)
    y_lowres = np.convolve(y, lsf, mode='same')

    return y_lowres


def vac_to_stand(wave_vac):
    """Convert vacuum wavelength (Ang) to standard wavelength in air since we're
    doing ground based stuff. 

	https://idlastro.gsfc.nasa.gov/ftp/pro/astro/vactoair.pro
    Equation from Prieto 2011 Apogee technical note
    and equation and parametersfrom Cidor 1996
    
    inputs: 
    -------
    wave_fac: 1D array, wavelength [A]

    outputs:
    -------
    
    """
    # eqn
    sigma2 = (1e4 / wave_vac) ** 2.
    fact = 1. + 5.792105e-2 / (238.0185 - sigma2) + \
           1.67917e-3 / (57.362 - sigma2)

    # return l0 / n which equals lamda
    return wave_vac / fact


def setup_band(x, x0=0, sig=0.3, eta=1):
    """
    give step function

    inputs:
    ------
    x0
    sig
    eta
    """
    y = np.zeros_like(x)

    ifill = np.where((x > x0 - sig / 2) & (x < x0 + sig / 2))[0]
    y[ifill] = eta

    return y


def gen_filter_profile(wavegrid, lam, width, throughput, mode='tophat', plot=False, savefig=False):
    '''
    Generates a distribution of throughputs as a function of wavelength centered at a specific wavelength

    :parameter
    wavegrid: Array of wavelengths (in nm) that the plot covers
    lam: The wavelength (in nm) that the distribution is centered around
    width: the FWHM of the distribution
    throughput: Maximum throughput
    modes: tophat, gaussian
    plot: whether a plot should be created or not
    savefig: whether the figure should be saved or not

    :returns
    filter_profile - returns the arrafilter_profile of distributed throughput as a function of the wavelength grid
    '''
    if mode == 'tophat':
        filter_profile = (np.where(abs(wavegrid - lam) <= width / 2, throughput, 0))

    elif mode == 'gaussian':  # 2.355 factor converts FWHM to stdev
        filter_profile = np.exp(-.5 * ((wavegrid - lam) / (width / 2.355)) ** 2) / (
                (width / 2.355) * np.sqrt(2 * np.pi))
        filter_profile = filter_profile / max(filter_profile) * throughput

    if plot:
        plt.plot(wavegrid, filter_profile)
        plt.ylim(0, 1)
        plt.ylabel('throughput')
        plt.xlabel('wavelength/nm')
        plt.xticks(np.arange(min(wavegrid), max(wavegrid) + 1, 1.0))
        plt.yticks(np.arange(min(filter_profile), max(filter_profile) + 0.1, 0.1))
        plt.axvline(lam)
        if savefig:
            plt.savefig(f'./figures/{mode}_plot_center{lam}_wid_{width}.png')
        plt.show()

    return filter_profile


def calc_flux(so: storage_object, exoplanet=True, exo_speed=0,
              stel_speed=0) -> list:
    '''
    This function applies the relevant doppler shifts to the relevant spectra (so due to motion of exoplanet
    there is a shift in the exoplanet spectra), re-interpolates the shifted spectra back into the original wave-grid
    and then calculates the total flux that is received on the earth by a telescope of a given area given a certain
    exposure time.
    Note: This function works for one speed. To use it for an array of speeds, enclose the function call in a loop.
    :parameter
    so: storage_object - storage_object class that stores all the data
    exoplanet: boolean - True if exoplanet spectra is being used, False otherwise
    :returns
    flux : the result of the integration
    flux_err: the photon noise for each flux calculated using sqrt(flux)
    '''
    area, exposure_time = so.const.hale_area, 2 * 3600  # 2 hour transit period
    original_xgrid = so.xgrid

    # accounts for stellar proper motion
    shift_stel = np.array(stel_speed * original_xgrid / constants.c)
    shifted_xgrid_stel = np.array(original_xgrid + shift_stel)
    tck_stel = interpolate.splrep(shifted_xgrid_stel, so.stel.s, k=2, s=0)
    stel_spec = interpolate.splev(original_xgrid, tck_stel, der=0, ext=1)

    # accounts for exoplanet motion
    shift_exo = np.array(exo_speed * original_xgrid / constants.c)
    shifted_xgrid_exo = np.array(
        original_xgrid + shift_exo)  # new wavelength grid that is shifted due to doppler effect
    if exoplanet == False:
        depth = [1] * len(original_xgrid)
    else:
        tck_exo = interpolate.splrep(shifted_xgrid_exo, so.exo.depth, k=2, s=0)
        depth = interpolate.splev(original_xgrid, tck_exo, der=0, ext=1)

    # integrates to give flux
    flux_arr = np.array(
        [np.trapz((stel_spec * depth * so.tel.s * profile), original_xgrid) * area * exposure_time for profile in
         so.hirax.hfp])

    # flux_err (photon noise) calculations
    flux_err_arr = np.array([np.sqrt(flux) for flux in flux_arr])
    return flux_arr, flux_err_arr

def resample(x, y, sig=0.3, dx=0, eta=1, mode='slow'):
    """
    resample using convolution

    x: wavelength array in nm
    y_in/y_out: two y arrays (evaluated at x) to resample, units in spectral density (e.g. photons/nm)

    sig in nanometers - width of bin, default 0.3nm
    dx - offset for taking first bin, defaul 0
    eta 0-1 for efficiency (amplitude of bin) default 1
    
    modes: slow, fast
    slow more accurate (maybe?), fast uses fft

    slow method uses trapz so slightly more accurate, i think? both return similar flux values

    """
    if mode == 'fast':
        dlam = np.median(np.diff(x))  # nm per pixel, most accurate if x is uniformly sampled in wavelength
        if sig <= dlam: raise ValueError('Sigma value is smaller than the sampling of the provided wavelength array')
        nsamp = int(sig / dlam)  # width of tophat
        tophat = eta * np.ones(nsamp)  # do i need to pad this?

        int_spec_oversample = dlam * signal.fftconvolve(y, tophat, mode='same')  # dlam integrating factor

        int_lam = x[int(nsamp / 2 + dx / dlam):][
                  ::nsamp]  # shift over by dx/dlam (npoints) before taking every nsamp point
        int_spec = int_spec_oversample[int(nsamp / 2 + dx / dlam):][::nsamp]

    elif mode == 'slow':
        i = 0
        int_lam, int_spec = [], []
        # step through and integrate each segment
        while i * sig / 2 + dx < np.max(x) - sig / 2 - np.min(x):  # check
            xcent = np.min(x) + dx + i * sig / 2
            tophat = setup_band(x, x0=xcent, sig=sig, eta=eta)  # eta throughput of whole system
            int_spec.append(integrate(x, tophat * y))
            int_lam.append(xcent)
            i += 1

    return int_lam, int_spec
