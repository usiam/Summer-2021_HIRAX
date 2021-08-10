import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib

matplotlib.style.use('seaborn-pastel')

def translate_eu(d_eu):
    """
    make eu dataset match columns of arxv
    """
    k_eu = np.array(['name', 'radius', 'semi_major_axis', 'orbital_period', 'mass', 'mag_v', 'star_teff', 'star_mass',
                     'star_radius', 'ra', 'dec'])
    k_arxv = np.array(
        ['pl_name', 'pl_radj', 'pl_orbsmax', 'pl_orbper', 'pl_massj', 'sy_vmag', 'st_teff', 'st_mass', 'st_rad', 'ra',
         'dec'])

    mapper = {}
    for i, key in enumerate(k_eu):
        mapper[key] = k_arxv[i]

    return d_eu.rename(columns=mapper)


def sync_catalogs():
    """
    http://exoplanet.eu/catalog/
    and exoplanet archive have different formats

    load both and sync into one dictionary
    """
    d_eu = pd.read_csv('data/exoplanet_catalogs/exoplanet.eu_catalog_hotjupiter_072621.csv', delimiter=',')
    d_arxv = pd.read_csv('data/exoplanet_catalogs/hotjupiters_northerhemisphere_072621.csv', delimiter=',', header=212)

    d_eu_new = translate_eu(d_eu)
    df = pd.merge(d_eu_new, d_arxv, 'outer',
                  on=['pl_name', 'pl_massj', 'pl_orbper', 'pl_radj', 'pl_orbsmax', 'st_teff', 'st_mass', 'st_rad', 'ra',
                      'dec', 'sy_vmag'], copy=False, indicator=True)
    hj = df.drop_duplicates('pl_name')

    return hj.reset_index()

def write_dec_20_data(path):
    d_eu = pd.read_csv('data/exoplanet_catalogs/exoplanet.eu_catalog_hotjupiter_072621.csv', delimiter=',')
    translate_eu(d_eu)
    cmplt_catalog = sync_catalogs()
    filtered_df = cmplt_catalog[cmplt_catalog['dec'] > -20]
    df_tosave = filtered_df.reset_index()
    df_tosave.to_csv(path)

def calc_absorbtion_signal(data):
    kb = 1.38064852e-23 #m2 kg s-2 K-1
    mu = 2.3 * 1.6605390e-27 # kg
    G = 6.67408e-11 #m3 kg-1 s-2

    # estimate Temp planet in kelvin
    Teq = (1/4.)**(1/4.) * data['st_teff'] * np.sqrt(0.00465047 * data['st_rad']/data['pl_orbsmax']) # 0.00465047 AU per Rsun
    gravity = G * data['pl_massj'] * 1.898e27 / (data['pl_radj'] * 69911000.0)**2 #kg from m_jup

    # get H -> kb * T/ mu / g
    H = kb * Teq / mu / gravity # meters

    # calculate A like Sing did
    A = 2 * data['pl_radj']*0.10049 * H*1.4374e-9 / data['st_rad']**2
    return A

def write_file_with_no_nan_signal(data):
    signal = 2 * calc_absorbtion_signal(data)
    signal_percent = signal * 100
    data['signal_percent'] = signal_percent
    exoplanets = data[~np.isnan(data['signal_percent'])]
    exoplanets.to_csv('data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv')


def quad_fit(params, x):
    return params[0] * x ** 2 + params[1] * x + params[2]


if __name__=='__main__':

    def plot_vmag_signal_percent_noise_tfac(
            vmag_signal_data_file_name='data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv',
            config='config1', tfactor=np.arange(0.2, 1, 0.2), sigma = 3, savefig=False):

        # plots the exoplanets
        vmag_signal_data = pd.read_csv(vmag_signal_data_file_name)
        vmag = vmag_signal_data['sy_vmag']
        signal_percent = vmag_signal_data['signal_percent']
        labels = vmag_signal_data['pl_name']

        fig, ax = plt.subplots(figsize=(12, 10))
        ax.semilogx(signal_percent, vmag, '+k')
        ax.set_ylim((16, 7))
        ax.set_xlim(0.01, 0.1)
        ax.set_xlabel('Transit signal of 2 Atmospheric Scale Height (%)')
        ax.set_ylabel('V magnitude')
        ax.set_xticks([0.01, 0.1], minor=True)

        # writes the exoplanet names
        for i, txt in enumerate(labels):
            ax.annotate(txt, (signal_percent[i], vmag[i]), fontsize=5)

        # plots the noise lines for the throughput factors
        ax2 = ax.twiny()

        path = 'data/output/'
        file = 'amplitude_vmag_'
        extension = '.csv'
        noise_data = pd.read_csv(path + file + config + extension)

        for j, tfac in enumerate(tfactor):
            if j % 2:
                ax2.semilogx(sigma * noise_data[f"amplitude_err_tfac{j + 1}"], noise_data['magnitude'], '-',
                             label=f'tfac = {round(tfac, 3)}', alpha=0.8)
            else:
                ax2.semilogx(sigma * noise_data[f"amplitude_err_tfac{j + 1}"], noise_data['magnitude'], '--',
                             label=f'tfac = {round(tfac, 3)}', alpha=0.8)

            # finding functional form of the noise lines
            quad_params = np.polyfit(x=3 * noise_data[f"amplitude_err_tfac{j + 1}"], y=noise_data['magnitude'], deg=2)
            print(quad_params)
            exo_mags = quad_fit(quad_params, 3 * noise_data[f"amplitude_err_tfac{j + 1}"])
            ax2.semilogx(3 * noise_data[f"amplitude_err_tfac{j + 1}"], exo_mags, '-k')

        ax2.axis('off')
        ax2.legend()
        fig.suptitle(file.split('_')[2])
        if savefig:
            fig.savefig(f"figures/vmag_signal_percent_{config}.png", dpi=600)


    plot_vmag_signal_percent_noise_tfac(config='config3', savefig=False)
    # exoplanets = pd.read_csv('data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv')
    # signal = exoplanets['signal_percent']
    #
    # path = 'data/output/'
    # file = 'amplitude_vmag_config1'
    # extension = '.csv'
    # noise_data = pd.read_csv(path + file + extension)
    #
    # tfactor = np.arange(0.2, 1, 0.2)
    # x = []
    # for ind, tfac in enumerate(tfactor):
    #     noise_3s = 3 * noise_data[f'amplitude_err_tfac{ind + 1}']
    #     mag = noise_data['magnitude']
    #     quad_params = np.polyfit(x=noise_3s, y=mag, deg=2)
    #     exo_mags = quad_fit(quad_params, signal)
    #     x.append(exo_mags)









