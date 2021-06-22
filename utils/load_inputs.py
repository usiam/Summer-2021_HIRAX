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
    spec = f[0].data / (1e8)  # ergs/s/cm2/cm to ergs/s/cm2/Angstrom for conversion
    f.close()

    path = stelname.split('/')
    f = fits.open(path[0] + '/' + path[1] + '/' + path[2] + '/' + \
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
    so.filt.v, so.filt.s = load_phoenix(so.stel.phoenix_file, wav_start=470.0, wav_end=700.0)  # nm, phot/m2/s/nm
    xtemp, ytemp = np.loadtxt(so.filt.filter_file).T  # nm, transmission out of 1
    f = interpolate.interp1d(xtemp / 10, ytemp, bounds_error=False, fill_value=0)
    so.filt.x, so.filt.y = so.filt.v, f(so.filt.v)  # filter profile sampled at stellar
    filtered_stellar = so.filt.s * so.filt.y  # filter profile resampled to phoenix times phoenix flux density

    so.filt.dl_l = np.mean(integrate(so.filt.x, so.filt.y) / so.filt.x)
    nphot_expected_0 = calc_nphot(so.filt.dl_l, so.filt.zp_v,
                                  vmag)  # what's the integrated V flux supposed to be in photons/m2/s?
    nphot_phoenix = integrate(so.filt.v, filtered_stellar)  # what's the integrated V flux now? in same units as ^

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
        self.x = np.arange(so.const.l0, so.const.l1, 0.00005)
        self.hk(so)
        self.stellar(so)
        self.kpf(so)
        self.telluric(so)

    def hk(self, so):
        ###########
        # load hk transmission
        xtemp, ytemp = np.loadtxt(so.hk.transmission_file).T
        f = interp1d(xtemp, ytemp, kind='linear', bounds_error=False, fill_value=0)

        so.hk.xtransmit, so.hk.ytransmit = self.x, f(self.x)

    def stellar(self, so):
        """
        loads stellar spectrum
        returns spectrum scaled to input V band mag

        everything in nm
        """
        # Part 1: load raw spectrum
        #
        teff = str(int(so.var.teff)).zfill(5)
        so.stel.phoenix_file = so.stel.phoenix_folder + 'lte%s-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits' % (teff)
        so.stel.vraw, so.stel.sraw = load_phoenix(so.stel.phoenix_file, wav_start=so.const.l0,
                                                  wav_end=so.const.l1)  # phot/m2/s/nm

        # find scaling factor
        so.stel.factor_0 = scale_stellar(so, so.var.vmag)

        # Scale stellar spectrum to v
        tck_stel = interpolate.splrep(so.stel.vraw, so.stel.sraw, k=2, s=0)
        so.stel.s = so.stel.factor_0 * interpolate.splev(self.x, tck_stel, der=0, ext=1)
        so.stel.v = self.x
        so.stel.units = 'photons/s/m2/nm'  # stellar spec is in photons/s/m2/nm

    def kpf(self, so):
        """
        load kpf things
        """
        # kpf transmission
        xtemp, ytemp = np.loadtxt(so.kpf.transmission_file, delimiter=',').T
        f = interp1d(xtemp, ytemp, kind='linear', bounds_error=False, fill_value=0)

        so.kpf.xtransmit, so.kpf.ytransmit = self.x, f(self.x)

    def telluric(self, so):
        """
        load tapas telluric file
        """
        data = fits.getdata(so.tel.telluric_file)
        tck_tel = interpolate.splrep(data['Wave/freq'], data['Total'], k=2, s=0) # the documentation says even k should be avoided for small 0
        so.tel.v, so.tel.s = self.x, interpolate.splev(self.x, tck_tel, der=0, ext=1)

    def exoplanet(self, so):
        '''
        loads the goyal exoplanet file
        '''

        goyal_file = np.loadtxt('./data/goyal/trans-eqpt_WASP-019_0.25_+1.7_0.75_model.txt')
        v_temp, depth_temp = (goyal_file[~np.isnan(goyal_file).any(axis=1), :]).T
        tck_exo = interpolate.splrep(v_temp, depth_temp, k=2, s=0)
        so.exo.v, so.exo.depth = self.x, interpolate.splev(self.x, tck_exo, der=0, ext=1)

    # look into website, and paper
    # look into interpolate and set it up similar to telluric -- DONE??

    def hirax(self, so): # we have to call gen_filter_profile in this function possibly?
        '''
        loads the hirax file
        '''
        hirax_file = np.loadtxt('./data/hirax/hirax_bandpass.txt')
        lam, width, throughput = hirax_file.T
        so.hirax.wavetime, so.hirax.hfp = self.x, gen_filter_profile(self.x, lam, width, throughput, savefig=True)
