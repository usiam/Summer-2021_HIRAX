import os

os.sys.path.append('./utils/')
from utils.objects import load_object, storage_object
from utils.load_inputs import fill_data
from utils.make_plots import MakePlots
from utils import functions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit

font = {'size': 8}
matplotlib.rc('font', **font)
matplotlib.style.use('seaborn-pastel')
plt.ion()

def lorentz(x, a, w, c, k):
    '''
    :param a: amplitude
    :param w: width
    :param c: center
    :param k: constant
    :return: function output
    '''
    return ((a * w)/(w**2 + (x - c)** 2)) + k

def fit_lorentz(x_data, y_data, y_err, init_guess):
    '''
    :param x_data: x data to be fitted
    :param y_data: y data to be fitted
    :param y_err: error in y data to be fitted
    :param init_guess: initial guess as a list
    :return: fitted parameters in the initial guess
    '''
    popt, pcov = curve_fit(lorentz, x_data, y_data, sigma=y_err, p0=init_guess,
                           absolute_sigma=True)
    return popt, pcov


def lorentz_sodium(x, a, w):
    '''
    :param a: amplitude
    :param w: width
    :param c: center
    :param k: constant
    :return: function output
    '''
    c = 5.89610638 * 10**2
    k = 9.80597083 * 10**(-1)
    return ((a * w)/(w**2 + (x - c)** 2)) + k

def lorentz_potassium(x, a, w):
    '''
    :param a: amplitude
    :param w: width
    :param c: center
    :param k: constant
    :return: function output
    '''
    c = 7.66907675e+02
    k = 9.80710064e-01
    return ((a * w)/(w**2 + (x - c)** 2)) + k

def fit_lorentz_aw(x_data, y_data, y_err, init_guess, feature = 'sodium'):
    '''
    :param x_data: x data to be fitted
    :param y_data: y data to be fitted
    :param y_err: error in y data to be fitted
    :param init_guess: initial guess as a list
    :return: fitted parameters in the initial guess
    '''
    if feature == 'sodium':
        popt, pcov = curve_fit(lorentz_sodium, x_data, y_data, sigma=y_err, p0=init_guess,
                           absolute_sigma=True)
    elif feature == 'potassium':
        popt, pcov = curve_fit(lorentz_potassium, x_data, y_data, sigma=y_err, p0=init_guess,
                           absolute_sigma=True)
    return -popt, pcov

if __name__ == '__main__':
    # load inputs
    configfile = 'hirax_snr_sodium.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    magnitudes = np.arange(7, 17, 1)
    tfac = np.arange(0.2 , 1, 0.2)
    a = [[]*len(magnitudes) for i in range(len(tfac))]
    a_err = [[]*len(magnitudes) for i in range(len(tfac))]
    w = [[]*len(magnitudes) for i in range(len(tfac))]
    w_err = [[]*len(magnitudes) for i in range(len(tfac))]
    default = so.var.vmag
    for ind,fac in enumerate(tfac):
        for mag in magnitudes:
            so.var.vmag = mag
            cload.stellar(so)
            cload.hirax(so, fac)
            ratio, error_in_ratio = functions.calc_flux_ratio(so)
            popt, pcov = fit_lorentz_aw(so.hirax.center_lam, ratio, y_err=error_in_ratio, init_guess=[-0.00035, 3], feature=configfile.split('_')[2].split('.')[0])
            a[ind].append(popt[0]/np.sqrt(pcov[0][0]))
            w[ind].append(popt[1] / np.sqrt(pcov[1][1]))
            a_err[ind].append(np.sqrt(pcov[0][0]))
            w_err[ind].append(np.sqrt(pcov[1][1]))
    so.var.vmag = default

    make_plots = MakePlots(so)
    make_plots.plot_SNR_mag_tput(magnitudes=magnitudes, amplitudes=a, tfac = np.arange(0.2 , 1, 0.2), savefig=True)

