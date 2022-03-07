##############################################################
# Load variables into objects object
###############################################################

import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from astropy.io import fits
from scipy import interpolate
import sys

from functions import *

__all__ = ['fill_data', 'load_phoenix']


def load_phoenix(stelname, wav_start=750, wav_end=780):
    """
    load fits file stelname with stellar spectrum from phoenix
    http://phoenix.astro.physik.uni-goettingen.de/?page_id=15

    return subarray

    wav_start, wav_end specified in nm

    convert s from egs/s/cm2/cm to phot/cm2/s/nm using
    https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf
    """

    # conversion factor

    f = fits.open(stelname)
    # ergs/s/cm2/cm to ergs/s/cm2/Angstrom for conversion
    spec = f[0].data / (1e8)
    f.close()

    path = stelname.split('/')
    f = fits.open(path[0] + '/' + path[1] + '/' + path[2] + '/' +
                  'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
    lam = f[0].data  # angstroms
    f.close()

    # Convert
    conversion_factor = 5.03 * 10 ** 7 * lam  # lam in angstrom here
    spec *= conversion_factor  # phot/cm2/s/angstrom

    # Take subarray requested
    isub = np.where((lam > wav_start * 10.0) & (lam < wav_end * 10.0))[0]

    # Convert
    return lam[isub] / 10.0, spec[isub] * 10 * 100 ** 2  # nm, phot/m2/s/nm


def calc_nphot(dl_l, zp, mag):
    """
    http://astroweb.case.edu/ssm/ASTR620/mags.html

    Values are all for a specific bandpass, can refer to table at link ^ for values
    for some bands. Function will return the photons per second per meter squared
    at the top of Earth atmosphere for an object of specified magnitude

    inputs:
    -------
    dl_l: float, delta lambda over lambda for the passband
    zp: float, flux at m=0 in Jansky
    mag: stellar magnitude

    outputs:
    --------
    photon flux
    """
    phot_per_s_m2_per_Jy = 1.51 * 10 ** 7  # convert to phot/s/m2 from Jansky

    return dl_l * zp * 10 ** (-0.4 * mag) * phot_per_s_m2_per_Jy


def scale_stellar(so, vmag):
    """
    scale spectrum by Vmag
    """
    so.filt.v, so.filt.s = load_phoenix(
        so.stel.phoenix_file, wav_start=470.0, wav_end=700.0)  # nm, phot/m2/s/nm
    # nm, transmission out of 1
    xtemp, ytemp = np.loadtxt(so.filt.filter_file).T
    f = interpolate.interp1d(
        xtemp / 10, ytemp, bounds_error=False, fill_value=0)
    # filter profile sampled at stellar
    so.filt.x, so.filt.y = so.filt.v, f(so.filt.v)
    # filter profile resampled to phoenix times phoenix flux density
    filtered_stellar = so.filt.s * so.filt.y

    so.filt.dl_l = np.mean(integrate(so.filt.x, so.filt.y) / so.filt.x)
    nphot_expected_0 = calc_nphot(so.filt.dl_l, so.filt.zp_v,
                                  vmag)  # what's the integrated V flux supposed to be in photons/m2/s?
    # what's the integrated V flux now? in same units as ^
    nphot_phoenix = integrate(so.filt.v, filtered_stellar)

    return nphot_expected_0 / nphot_phoenix


class fill_data():
    """
    Load variables into storage object

    Inputs: so (storage object with user defined things loaded)
    Outputs: so (storage object with data and other stuff loaded)

    Edits
    -----
    Ashley - initial implementation Oct 26, 2018
    """

    def __init__(self, so):
        # define
        so.xgrid = np.arange(so.const.l0, so.const.l1, 0.00005)
        self.x = so.xgrid
        self.stellar(so)
        self.exoplanet(so)
        self.hirax(so)
        self.telluric(so)
        self.oh(so)

    def stellar(self, so):
        """
        loads stellar spectrum
        returns spectrum scaled to input V band mag

        everything in nm
        """
        # Part 1: load raw spectrum
        #
        teff = str(int(so.var.teff)).zfill(5)
        so.stel.phoenix_file = so.stel.phoenix_folder + \
            'lte%s-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' % (teff)
        so.stel.vraw, so.stel.sraw = load_phoenix(so.stel.phoenix_file, wav_start=so.const.l0,
                                                  wav_end=so.const.l1)  # phot/m2/s/nm

        # find scaling factor
        so.stel.factor_0 = scale_stellar(so, so.var.vmag)

        so.stel.speed = (so.var.stel_speed +
                         so.var.barycentric_speed) * 10 ** 3

        # Scale stellar spectrum to v
        tck_stel = interpolate.splrep(so.stel.vraw, so.stel.sraw, k=2, s=0)

        so.stel.s = so.stel.factor_0 * \
            interpolate.splev(self.x, tck_stel, der=0, ext=1)

        # adding vel shift
        shift_spectra = np.array(so.stel.speed * so.xgrid / constants.c)
        # new wavelength grid that is shifted due to doppler effect
        shifted_xgrid = np.array(so.xgrid + shift_spectra)
        tck_spec = interpolate.splrep(shifted_xgrid, so.stel.s, k=2, s=0)
        so.stel.s = interpolate.splev(
            so.xgrid, tck_spec, der=0, ext=1)

        so.stel.v = self.x

        so.stel.units = 'photons/s/m2/nm'  # stellar spec is in photons/s/m2/nm

    def telluric(self, so):
        """
        load tapas telluric file
        """
        data = fits.getdata(so.tel.telluric_file)
        tck_tel = interpolate.splrep(data['Wave/freq'], data['Total'], k=2,
                                     s=0)  # the documentation says even k should be avoided for small 0
        so.tel.v, so.tel.s = self.x, interpolate.splev(
            self.x, tck_tel, der=0, ext=1)

    def exoplanet(self, so):
        '''
        loads the goyal exoplanet file
        '''

        goyal_file = np.loadtxt(so.exo.exoplanet_file)
        v_temp, depth_temp = (
            goyal_file[~np.isnan(goyal_file).any(axis=1), :]).T
        tck_exo = interpolate.splrep(v_temp * 1000, 1 - depth_temp, k=2, s=0)
        so.exo.v, so.exo.depth = self.x, interpolate.splev(
            self.x, tck_exo, der=0, ext=1)
        so.exo.speed = np.arange(so.var.min_exo_speed,
                                 so.var.max_exo_speed + 50, 50) * 10 ** 3

    def hirax(self, so, tfac=1):  # do you want me to modify this function a bit? I could do something like
        '''
        loads the hirax file
        '''
        hirax_file = np.loadtxt(so.hirax.hirax_file)
        lam, width, throughput = hirax_file.T
        so.hirax.wavegrid = self.x
        so.hirax.width = width
        so.hirax.throughput = throughput * tfac
        so.hirax.center_lam = lam
        so.hirax.hfp = [gen_filter_profile(self.x, l, w, t) for l, w, t in
                        zip(lam, width, so.hirax.throughput)]

    def oh(self, so):
        '''
        loads the oh file
        '''
        oh_file = np.loadtxt(so.oh.oh_file)
        vraw, sraw = (1*10**7)/(oh_file.T[0])[::-1], ((oh_file.T[1])[::-1])
        tck_oh = interpolate.splrep(vraw, sraw, k=2, s=0)
        so.oh.v, so.oh.s = self.x, interpolate.splev(
            self.x, tck_oh, der=0, ext=1)
