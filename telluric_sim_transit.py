import os
os.sys.path.append('./utils/')
from utils.objects import load_object
from utils.load_inputs import fill_data
from utils.make_plots import MakePlots
import numpy as np
from astropy.time import Time
import matplotlib.pylab as plt
from scipy import signal
#from pytransit import QuadraticModel
import astropy.units as u
import astropy.constants as c
from scipy import interpolate
import matplotlib

font = {'size': 12}
matplotlib.rc('font', **font)
matplotlib.style.use('seaborn-pastel')



def define_lsf(v,res):
    """
    define gaussian in pixel elements to convolve resolved spectrum with to get rightish resolution
    """
    dlam  = np.median(v)/res
    fwhm = dlam/np.mean(np.diff(v)) # desired lambda spacing over current lambda spacing resolved to give sigma in array elements
    sigma = fwhm/2.634 # FWHM is dl/l but feed sigma to gaussian equation 
    x = np.arange(sigma*10)
    gaussian = (1./sigma/np.sqrt(2*np.pi)) * np.exp(-0.5*( (x - 0.5*len(x))/sigma)**2 )

    return gaussian

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

    ifill = np.where((x > x0-sig/2) & (x < x0 + sig/2))[0]
    y[ifill] = eta

    return y

def resample_two(x,y_in,y_out,sig=0.3, dx=0, eta=1,mode='slow'):
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
    if mode=='fast':
        dlam= np.median(np.diff(x)) # nm per pixel
        nsamp = int(sig / dlam)     # width of tophat
        temp_band   = eta * np.ones(nsamp)

        int_spec_in_oversample         = dlam * signal.fftconvolve(y_in,temp_band,mode='same') # dlam integrating factor
        int_spec_out_oversample        = dlam * signal.fftconvolve(y_out,temp_band,mode='same') # dlam integrating factor
        int_lam = x[int(nsamp/2 + dx/dlam):][::nsamp] # shift over by dx/dlam (npoints) before taking every nsamp point

        int_spec_in , int_spec_out = int_spec_in_oversample[int(nsamp/2 + dx/dlam):][::nsamp], int_spec_out_oversample[int(nsamp/2 + dx/dlam):][::nsamp]

    elif mode=='slow':
        i=0
        int_lam, int_spec_in, int_spec_out = [], [], []
        # step through and integrate each segment
        while i*sig/2 + dx< np.max(x)-sig/2 - np.min(x): # check
            xcent = np.min(x) + dx + i*sig/2
            temp_band   = setup_band(x, x0=xcent, sig=sig, eta=eta) # eta throughput of whole system
            int_spec_in.append(integrate(x,temp_band * y_in))
            int_spec_out.append(integrate(x,temp_band * y_out))
            int_lam.append(xcent)
            i += 1

    return int_lam, int_spec_in, int_spec_out

def calc_planet_vperp1(t, t_mid, ms, mp, period,v_star):
    """
    take star RV equation and scale by sqrt(ms/mp) to get the planet's velocity

    Inputs:
    -------
    t: time array (hours)
    t_mid: time of transit midpoint (hrs)
    v_star - stellar RV orbital max amplitude , meters/second
    ms - stars mass (solar masses)
    mp - planets mass (jupiter masses)
    period - orbital period (days)

    assumes circular orbit
    """
    period, v_star = period * u.day, v_star * u.m/u.s# v_star is orbital max amplitude only
    mp, ms     =  mp * u.jupiterMass, ms * u.solMass
    v_planet   = v_star * u.m/u.s  *  ms/mp
    phase      = u.rad * (t - t_mid) * 2 * np.pi /period
    v_pl       = v_planet * np.sin(phase.decompose()) # eh i think this is right

    return v_pl.decompose()

def calc_planet_vperp2(t, t_mid, ms,a, period):
    """
    take star RV equation and scale by sqrt(ms/mp) to get the planet's velocity

    Inputs:
    -------
    t: time array (hours)
    t_mid: time of transit midpoint (hrs)
    a - semi major axis (AU)
    ms - stars mass (solar masses)
    period - orbital period (days)

    assumes circular orbit
    """
    v_planet  =  np.sqrt(c.G * ms / a)
    period    = period * u.day 
    phase     = u.rad * (t - t_mid) * 2 * np.pi /period
    v_pl      = v_planet * np.sin(phase.decompose()) # eh i think this is right

    return v_pl.decompose()

def sim_bary_vel():
    """
    simulate barycentric velocity shift during transit for getting tellurics right
    using ephemeris file from IAG data for sun just to get the velocity change roughly right
    later can add offset to take different parsts of the curve

    returns: interpolation function, takes time in hours between 0 and 24
    """
    f = open('data/ephem/2015-03-23_SunI2.txt','r')
    lines = f.readlines()
    f.close()

    k = []
    ts = []
    read = False
    for l in lines:
        if l.startswith('$$EOE'):
            read = False
        if read:
            k.append(float(l.split('  ')[-3]))
            ts.append(l.split('  ')[0][1:18])
        if l.startswith('$$SOE'):
            read = True

    months = {'Jan': '01', 'Feb': '02', 'Mar': '03', \
              'Apr': '04', 'May': '05', 'Jun': '06', \
              'Jul': '07', 'Aug': '08', 'Sep': '09', \
              'Oct': '10', 'Nov': '11', 'Dec': '12'}
    
    for i,t in enumerate(ts):
        if t[5:8] in months.keys():
            ts[i] = ts[i].replace(t[5:8],months[t[5:8]])
            ts[i] = ts[i] + ':00.00'
    
    tarr = Time(ts, format='iso', scale='utc').jd

    # interpolate k on tarr array
    int_ks = interpolate.interp1d((tarr - tarr[0])*24,k,kind='quadratic')

    return int_ks

def calc_flux_ref(so,sig,eta,amass1=1.0, amass2=1.0, res=4000, grism=False,sample_mode='fast'):
    """
    calc flux for reference star
    """
    source = so.ref.s # for using different spectral type

    # shifted stel spectrum
    #vel = 20
    #x_temp        = so.exo.v * (1 +(1.0 * v_star/300000.0))
    #source        = np.interp(so.exo.v,x_temp,so.stel.s,left=1,right=1)

    at_scope1     = source     * (so.tel.h2o* so.tel.o2*so.tel.rayleigh)**amass1
    at_ccd1       = at_scope1   * so.inst.exp_time   * so.inst.tel_area

    at_scope2     = source     * (so.tel.h2o* so.tel.o2*so.tel.rayleigh)**amass2
    at_ccd2       = at_scope2   * so.inst.exp_time   * so.inst.tel_area

    if grism == True:
        # convolve w/ PSF first - can make PFS variable
        lsf         = define_lsf(so.exo.v,res=res)
        at_ccd1     = np.convolve(at_ccd1,lsf,mode='same')
        at_ccd2     = np.convolve(at_ccd2,lsf,mode='same')  

    # resample data
    int_lam, int_spec1, int_spec2 = resample(so.exo.v,at_ccd1,at_ccd2,sig=sig,dx=0,eta=eta,mode=sample_mode)

    return int_spec1, int_spec2

def make_transit(k=0.1, a=12, p=3.86):
    """
    Inputs:
    k (float or array)
        - planet to star radius ratio (or average over bandpass considering)
    p (float):
         - period of exoplanet orbit in days
    a (float)
         - scaled semi-major axis (a/R*), unitless

    default values are for WASP-69 

    Assumptions:
    zero-epoch t0, orbital period p, scaled semi-major axis a (a/R*), orbital inclination i, eccentricity e, and argument of periastron w), ldc limb darkening coefficients

    To DO:
    make a, p free params based on exoplanet catalog
    """
    tm = QuadraticModel()
    t_dur   = 3/24. # apprxo duration, days, for time range
    times   = np.arange(-2*t_dur,2*t_dur,t_dur/100)
    tm.set_data(times)#,exptimes=3/3600/24, nsamples=10) # use exptimes and nsamples if running long exp times
    transit = tm.evaluate(k=k, ldc=[0.2, 0.1], t0=0.0, p=p, a=a, i=0.5*np.pi, w=0)
    if plot: plt.figure(-104); plt.plot(times,transit); plt.xlabel('Time (days)'); plt.ylabel('Flux')

    return transit


def plot_vel_transit(times, vpl, transit):
    """
    plot transit and velocity of planet during transit
    """
    fig, ax = plt.subplots(1,1,figsize=(6,4))#, sharex=True, sharey=False)
    ax.plot(times,transit)
    ax1 = ax.twinx()
    ax1.plot(times,v_pl.decompose(),c='orange')

def resample(x,y,sig=0.3, dx=0, eta=1,mode='slow'):
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
    if mode=='fast':
        dlam    = np.median(np.diff(x)) # nm per pixel, most accurate if x is uniformly sampled in wavelength
        if sig <= dlam: raise ValueError('Sigma value is smaller than the sampling of the provided wavelength array')
        nsamp   = int(sig / dlam)     # width of tophat
        tophat  = eta * np.ones(nsamp) # do i need to pad this?

        int_spec_oversample    = dlam * signal.fftconvolve(y,tophat,mode='same') # dlam integrating factor
        
        int_lam  = x[int(nsamp/2 + dx/dlam):][::nsamp] # shift over by dx/dlam (npoints) before taking every nsamp point
        int_spec =  int_spec_oversample[int(nsamp/2 + dx/dlam):][::nsamp]

    elif mode=='slow':
        i=0
        int_lam, int_spec  = [], []
        # step through and integrate each segment
        while i*sig/2 + dx< np.max(x)-sig/2 - np.min(x): # check
            xcent    = np.min(x) + dx + i*sig/2
            tophat   = setup_band(x, x0=xcent, sig=sig, eta=eta) # eta throughput of whole system
            int_spec.append(integrate(x,tophat * y))
            int_lam.append(xcent)
            i += 1

    return int_lam, int_spec


def calc_flux(so,sig, eta, vel_p=0, amass1=1.0, amass2=1.0, res=4000, grism=False,sample_mode='slow'):
    """
    
    """
    #hack, let's add these to config
    so.stel.vel = 0
    so.exo.vel  = 0

    # exoplanet transmission is 1-(rp/rs)^2
    transmission  = so.exo.depth

    # shift star to user-defined velocity

    x_temp        = so.exo.v * (1 +(1.0 * so.stel.vel/300000.0))
    stel          = np.interp(so.exo.v,x_temp,so.stel.s,left=1,right=1)

    x_temp        = so.exo.v * (1 +(1.0 * (so.exo.vel + so.stel.vel)/300000.0))
    transmission  = np.interp(so.exo.v,x_temp,transmission,left=1,right=1)

    # compute source
    source     = transmission * stel  #
    source_out = stel

    # multiply by telluric
    at_scope     = source     * (np.abs(so.tel.h2o* so.tel.s))**amass1
    #at_scope     = source     * so.tel.h2o * so.tel.o2 * so.tel.rayleigh
    at_scope_out = source_out * (np.abs(so.tel.h2o * so.tel.s))**amass2

    # Instrument
    at_ccd       = at_scope     * so.var.transit_duration   * so.const.tel_area
    at_ccd_out   = at_scope_out * so.var.transit_duration   * so.const.tel_area

    if grism == True:
        # convolve w/ PSF first - can make PFS variable
        lsf    = define_lsf(so.exo.v,res=res)
        at_ccd      = np.convolve(at_ccd,lsf,mode='same')
        at_ccd_out  = np.convolve(at_ccd_out,lsf,mode='same')   

    # resample data
    int_lam, int_spec_in, int_spec_out = resample(so.exo.v,at_ccd,at_ccd_out,sig=sig,dx=0,eta=eta,mode=sample_mode)

    # compute errors
    int_lam, int_spec_in, int_spec_out = np.array(int_lam), np.array(int_spec_in), np.array(int_spec_out)
    err_in = np.sqrt(int_spec_in)
    err_out = np.sqrt(int_spec_out)

    in_over_out     = int_spec_in/int_spec_out
    in_over_out_err = in_over_out * np.sqrt(1/int_spec_in + 1/int_spec_out)

    # put things into so.out. - make this a thing instead of afterthought*
    so.out.s_in  = at_ccd
    so.out.s_out = at_ccd_out
    so.out.s_in_obvs = int_spec_in
    so.out.s_in_obvs = int_spec_out 

    return int_lam, in_over_out, in_over_out_err

def get_differential_telluric(sig=0.3, eta=0.4, a1=1.0, a2=1.1, res=None, grism=False, sample_mode='fast'):
    """
    
    """
    atm1 = (np.abs(so.tel.h2o * so.tel.s))**a1
    atm2 = (np.abs(so.tel.h2o * so.tel.s))**a2   

    if grism:
        # convolve w/ PSF first - can make PFS variable
        lsf    = define_lsf(so.exo.v,res=res)
        atm1   = np.convolve(atm1,lsf,mode='same')
        atm2  = np.convolve(atm2,lsf,mode='same')   

    # resample data
    int_lam, spec1, spec2 = resample(so.exo.v,atm1,atm2,sig=sig,dx=0,eta=eta,mode=sample_mode)

    return np.array(spec1), np.array(spec2)

def plot_telluric_contamination(so):
    """
    plot exoplanet signal with no tellurics, one with tellurics but corrected to see differential error
    no velocity shift in tellurics assumed 
    """
    fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True,figsize=(10,6),gridspec_kw={'height_ratios': [2, 1]})

    # plot with raw signal
    lam, frat, frat_err = calc_flux(so, .01, eta, vel_p=0,\
                     amass1=a1, amass2=a1, res=None, grism=False,
                     sample_mode=mode)
    ax1.plot(lam, frat ,c='k',ls='-',label='Signal')

    lam, frat, frat_err = calc_flux(so, sig, eta, vel_p=0,\
                     amass1=a1, amass2=a1, res=res, grism=True,
                     sample_mode=mode)
    #ax1.errorbar(lam, frat , frat_err, c='r',fmt='o',label='Signal Recovered')

    # plot with telluric correction
    lam, frat, frat_err = calc_flux(so, sig, eta, vel_p=0,\
                     amass1=a1, amass2=a2, res=res, grism=grism,
                     sample_mode=mode)

    # again but with slightly wrong airmass
    tel1, tel2 = get_differential_telluric(sig=sig, eta=eta, a1=a1, a2=a2+0.02,\
                                            res=res, grism=grism, sample_mode=mode)    
    eb = ax1.errorbar(lam,frat * (tel2/tel1), frat_err, alpha=0.3)
    #ax1.plot(lam,frat,c=eb[0].get_color(),label=str(a2[i])) 
    ax1.plot(lam, frat * (tel2/tel1),c=eb[0].get_color(),ls='--',label='Observed')

    # divide out low resolution spectrum
    ax1.set_ylabel('1 - $(R_p/R_s)^2$')
    ax1.set_title('R~%s, Standard Star' %(res))

    ax1.set_ylabel('$F_{in}/F_{out}$')
    ax1.legend(fontsize=12)

    # plot solar + telluric below
    ax2.plot(so.exo.v,so.stel.s,c='orange')#,label='Stellar (%s)' %so.stel.type)
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel('Stellar Flux Density \nArea (photons/s/nm/m$^2$)',fontsize=14)
    ax3 = ax2.twinx()
    ax3.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2 
    ax3.set_ylabel('Sky Transmission')

    ax3.plot(so.exo.v,so.tel.s * so.tel.h2o,c='steelblue',label='Telluric')

    ax2.legend(loc=4,fontsize=12)
    ax3.legend(loc=3,fontsize=12)

    fig.subplots_adjust(bottom=0.13,left=0.22,hspace=0,right=0.88,top=0.9)


def doppler(v_over_c):
    """
    """
    return np.sqrt( (1 + v_over_c)/ (1 - v_over_c))

def shift_spectrum(x,y,v,xnew=None):
    """
    shift spectrum by v
    v must be in m/s
    return y
    type: type of interpolation, nearest works best for retaining noise structure

    if x in observer frame and you want to shift to star's frame, and star moving w.r.t. earth by -v
    then you want to shift x by +v to undo shift
    """
    if np.all(xnew == None):
        xnew = x

    v_over_c = v * (u.m/u.s) / c.c
    doppler_factor = doppler(v_over_c)

    x_shift = doppler_factor.value * x

    # interp onto new wavelength grid
    tck        = interpolate.splrep(x_shift,y, s=0, k=3)
    y_int      = interpolate.splev(xnew, tck, der=0,ext=1)

    return y_int


def gen_velocities():
    """
    get velocity shifts
    """
    pass

def get_transits(so):
    """
    integrate stellar + spectrum
    """
    pass


def scintillation():
    """
    make scintillation noise
    """
    pass

# notes:
# make stellar(v_bc) * telluric(a,pwv) spectra
# integrate over hirax bandpasses for 4 hours
# multiply by transit w/ depth for exoplanet radius
# can add stellar activity effects later (by hack of changing presumed stellar disk integrated spectrum?)
# or can simulate stellar activity during transit by changing core depth of sodium line

def degrade_spec(x, y, res):
    """
    given wavelength, flux array, and resolving power R, return  spectrum at that R
    """
    lsf = define_lsf(x, res=res)
    y_lowres = np.convolve(y, lsf, mode='same')

    return y_lowres


if __name__ == '__main__':
    # load inputs
    configfile = 'hirax_snr_sodium.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    # define time array, pwv array, airmass array
    ts    = np.arange(0,4,0.01) # 36 sec intervals
    pwvs  = 1.2 + 0.5 *np.sin(2*np.pi * (ts/ts[-1])) 
    amass = 1/np.cos(np.arange(0,np.pi/3,.01))

    telluric_notch = so.tel.s**2 * (1 - np.sum(so.hirax.hfp,axis=0))
    tel_lowres= degrade_spec(so.tel.v,telluric_notch,630/0.14)

    plt.figure()
    plt.plot(so.exo.v,tel_lowres,'b')
    plt.ylim(0.1,1.3)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.xlim(580,960)
    plt.show()

    # plot telluric notch with noise
    amass=3
    telluric_notch = so.tel.s**amass * (1 - np.sum(so.hirax.hfp,axis=0))

    exptime = 60
    specthru = 0.01 # spectrometer throughput
    spec = so.stel.s * telluric_notch * so.const.hale_area * exptime * specthru
    spec_lowres = degrade_spec(so.tel.v,spec,630/0.08)
    xnew, spec_lowres_resamp = resample(so.exo.v,spec_lowres,sig=0.01, dx=0, eta=1,mode='fast')

    plt.figure()
    #plt.plot(so.exo.v, spec_lowres) 
    plt.errorbar(xnew, spec_lowres_resamp, yerr=np.sqrt(spec_lowres_resamp))





