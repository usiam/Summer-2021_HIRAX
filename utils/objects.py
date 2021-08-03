import configparser
import numpy as np
from distutils.util import strtobool

all = {'storage_object', 'load_object'}


class storage_object():
    """
    Main storage object for organization
    """

    def __init__(self):
        # Classes
        self.xgrid = None
        self.run = RUN()
        self.const = CONSTANT()
        self.var = VARIABLE()
        self.filt = FILTER()
        self.stel = STELLAR()
        self.tel = TELLURIC()
        self.out = OUTPUT()
        self.hirax = HIRAX()
        self.exo = EXOPLANET()
        self.oh = OH()
        # non class things
        self.info = "see objects.py in utils/ for info"


class RUN():
    "star info and spectrum"

    def __init__(self):
        self.plot_prefix = None  # stellar spec file name
        self.savename = None  # wavelength like normal (should match exoplanet and be in standard wavelength)


class CONSTANT():
    "float values"

    def __init__(self):
        self.l0 = None  # starting wavelength
        self.l1 = None  # ending wavelength
        self.snr_kpf = None  # kpf snr
        self.res_hk = None
        self.tel_area = None
        self.hale_area = None
        self.focal_length = None # in meter
        self.pixel_size = None # in micrometer
        self.theta_s = None # in arcsec
        self.psf = None
        # wheres res_hk and tel_area


class VARIABLE():
    "float values, variables"

    def __init__(self):
        self.vmag = None  # v band magntidue
        self.teff = None  # K, 100K steps
        self.transit_duration = None # in seconds
        self.min_exo_speed = None # in km/s
        self.max_exo_speed = None # in km/s
        self.stel_speed = None # in km/s
        self.dark_n = None
        self.read_n = None
        self.saturation = None
        self.magnification = None

class FILTER():
    "float values"

    def __init__(self):
        self.x = None
        self.y = None
        self.zp_v = 3640  # Jy from http://astroweb.case.edu/ssm/ASTR620/mags.html


class STELLAR():
    "star info and spectrum"

    def __init__(self):
        # User optional define:
        self.phoenix_file = None  # stellar spec file name, **make this take temp value in future
        # Filled in by code:
        self.speed = None
        self.vraw = None  # wavelength like normal (should match exoplanet and be in standard wavelength)
        self.sraw = None  # spectrum


class TELLURIC():
    "telluric transmission file, static"

    def __init__(self):
        # User optional define:
        self.telluric_file = None  # spec file name
        # Filled in by code:
        self.v = None  # wavelength
        self.s = None  # spectrum


class OUTPUT():
    "output file info and storage"

    def __init__(self):
        # User optional defined
        self.savename = None  # output name
        # Filled in by code
        self.spectrum = None  # ca H&k spectrum


class EXOPLANET():
    "exoplanet transmission file, static"
    def __init__(self):
        self.exoplanet_file = None
        self.v = None  # wavelength
        self.depth = None  # transit depth
        self.speed = None


class HIRAX():
    "Hirax filter profile file, static"
    def __init__(self):
        self.hirax_file = None
        self.wavegrid = None
        self.center_lam = None
        self.throughput = None
        self.width = None
        self.hfp = None  # hfp : hirax_filter_profile


class OH():
    "OH transmission file, static"
    def __init__(self):
        self.oh_file = None
        self.v = None
        self.s = None

def LoadConfig(configfile, config={}):
    """
    Reads configuration file 'XXX.cfg'
    returns a dictionary with keys of the form
    <section>.<option> and the corresponding values
    """
    config = config.copy()
    cp = configparser.ConfigParser()
    cp.read(configfile)
    for sec in cp.sections():
        name = str(sec)
        for opt in cp.options(sec):
            config[name + "." + str(opt)] = str(
                cp.get(sec, opt)).strip()
    return config


def load_object(configfile):
    """
    Loads config file as dictionary using LoadConfig function
    Then loads stoar_object and fills in user-defined
    quantities
    """
    config = LoadConfig(configfile)
    so = storage_object()

    for key in config:
        s1, s2 = key.split('.')
        if s1 == 'const' or s1 == 'var':
            setattr(getattr(so, s1), s2, float(config[key]))
        else:
            setattr(getattr(so, s1), s2, config[key])
    return so

    # setattr(getattr(so,s1), s2, config[key]) general syntax: setattr(object, attribute, value)
    # what this does is - getattr(so,s1) returns an attribute (s1) of the storage_object so like run, constant etc
    # then s2 is the attribute of s1 that we will set a value to. So if s1 is run then s2 can be plot_prefix and savename
    # finally config[key] is the value that s2 of object s1 is assigned
