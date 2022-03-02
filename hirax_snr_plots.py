import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from utils import functions
from utils.make_plots import MakePlots
from utils.load_inputs import fill_data
from utils.objects import load_object, storage_object
import os

os.sys.path.append('./utils/')

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
    return ((a * w)/(w**2 + (x - c) ** 2)) + k


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
    return ((a * w)/(w**2 + (x - c) ** 2)) + k


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
    return ((a * w)/(w**2 + (x - c) ** 2)) + k


def fit_lorentz_aw(x_data, y_data, y_err, init_guess, feature='sodium'):
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

# Check out scipy.minimize to automate the code


if __name__ == '__main__':
    # load inputs
    configfile = 'hirax_snr_potassium.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    mag = np.arange(7, 17, 1)
    tfactor = np.arange(0.2, 1, 0.2)

    def calc_amplitude_width(magnitudes, tfac):
        a = [[]*len(magnitudes) for i in range(len(tfac))]
        a_err = [[]*len(magnitudes) for i in range(len(tfac))]
        w = [[]*len(magnitudes) for i in range(len(tfac))]
        w_err = [[]*len(magnitudes) for i in range(len(tfac))]
        default = so.var.vmag
        for ind, fac in enumerate(tfac):
            for mag in magnitudes:
                so.var.vmag = mag
                cload.stellar(so)
                cload.hirax(so, fac)
                ratio, error_in_ratio = functions.calc_flux_ratio(so)
                popt, pcov = fit_lorentz_aw(so.hirax.center_lam, ratio, y_err=error_in_ratio, init_guess=[
                                            -0.00035, 3], feature=configfile.split('_')[2].split('.')[0])
                a[ind].append(popt[0])
                w[ind].append(popt[1])
                a_err[ind].append(np.sqrt(pcov[0][0]))
                w_err[ind].append(np.sqrt(pcov[1][1]))

        so.var.vmag = default
        return np.array(a), np.array(a_err), abs(np.array(w)), np.array(w_err)

    def write_snr_data(amplitudes, amplitude_errs, widths, width_errs, magnitudes, tfac):
        mydict = {
            'magnitude': magnitudes,
        }
        for i in range(len(tfac)):
            mydict[f'amplitude_tfac{i+1}'] = amplitudes[i]
            mydict[f'amplitude_err_tfac{i+1}'] = amplitude_errs[i]
            mydict[f'width_tfac{i+1}'] = widths[i]
            mydict[f'width_err_tfac{i+1}'] = width_errs[i]
        df = pd.DataFrame(mydict)
        df.to_csv(
            f"data/output/amplitude_vmag_{so.hirax.hirax_file.split('_')[2].split('.')[0]}.csv", sep=',', index=False)

    a, a_err, w, w_err = calc_amplitude_width(magnitudes=mag, tfac=tfactor)
    snr = a/a_err
    write_snr_data(amplitudes=a, amplitude_errs=a_err, widths=w,
                   width_errs=w_err, magnitudes=mag, tfac=tfactor)
