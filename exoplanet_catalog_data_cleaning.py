import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def make_one_catalog(d_eu, d_arxv):
    """
    make translation b/n one header to another
    for exoplanet archive and exoplanet eu catalogs
    inputs:
    d_eu - exoplanet.eu catalog dictionary
    d_arxv - exoplanet archive download file

    return:
    dic - pandas dataframe, combined dictionary of inputs with only needed parameters out

    Note:
    - no RV amplitude or stellar velocity in exoplanet.eu dataset
    - for description of variables see exoplanet archive csv file header
    """
    dic = {}

    # define things we need
    dtype = np.array(
        ['|S15', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float',
         'float', 'float'])
    keys = np.array(['NAME', 'R', 'MSINI', 'A', 'PER', 'K', 'MASS', 'V', 'TEFF', 'MSTAR', 'RSTAR', 'RA', 'DEC', 'RV'])

    k_eu = np.array(
        ['name', 'radius', 'mass_sini', 'semi_major_axis', 'orbital_period', np.nan, 'mass', 'mag_v', 'star_teff',
         'star_mass', 'star_radius', 'ra', 'dec', np.nan])
    k_arxv = np.array(
        ['pl_name', 'pl_radj', 'pl_massj', 'pl_orbsmax', 'pl_orbper', 'pl_rvamp', 'pl_massj', 'sy_vmag', 'st_teff',
         'st_mass', 'st_rad', 'ra', 'dec', 'st_radv'])

    for i, key in enumerate(keys):
        if k_eu[i] == 'nan':
            dic[key] = np.concatenate((np.nan * np.ones_like(d_eu['name']), d_arxv[k_arxv[i]]))
        else:
            dic[key] = np.concatenate((d_eu[k_eu[i]], d_arxv[k_arxv[i]]))

    return pd.DataFrame.from_dict(dic)


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
    # f1 = 'exoplanet.eu_catalog_hotjupiter_072621.csv'
    # f2 = 'hotjupiters_northerhemisphere_072621.csv'

    d_eu = pd.read_csv('data/exoplanet_catalogs/exoplanet.eu_catalog_hotjupiter_072621.csv', delimiter=',')
    d_arxv = pd.read_csv('data/exoplanet_catalogs/hotjupiters_northerhemisphere_072621.csv', delimiter=',', header=212)

    # dic_all = make_one_catalog(d_eu, d_arxv)
    # dic = dic_all.drop_duplicates('NAME')
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

if __name__=='__main__':

    exoplanets = pd.read_csv('data/exoplanet_catalogs/exoplanets_dec_over_20.csv')
    signal = 2 * calc_absorbtion_signal(exoplanets)
    vmag = exoplanets['sy_vmag']

    fig, ax = plt.subplots()
    ax.semilogx(signal * 100, vmag, '+k')
    ax.set_ylim((15, 7))
    ax.set_xlim(0.01, 1)
    ax.set_xlabel('Transit signal of 2 Atmospheric Scale Height (%)')
    ax.set_ylabel('V magnitude')

    ax2 = ax.twiny()
    path = 'data/output/'
    file = 'amplitude_vmag_config4'
    extension = '.txt'
    noise = pd.read_csv(path + file + extension)
    tfactor = np.arange(0.2 , 1, 0.2)

    for j,tfac in enumerate(tfactor):
        if j % 2:
            ax2.semilogx(noise[f"amplitude_err_tfac{j + 1}"], noise['magnitude'], '-k', label=f'tfac = {round(tfac, 3)}', alpha=0.5)
        else:
            ax2.semilogx(noise[f"amplitude_err_tfac{j + 1}"], noise['magnitude'], '--k', label=f'tfac = {round(tfac, 3)}', alpha=0.5)

    ax2.legend()
    fig.suptitle(file.split('_')[2])
    fig.savefig(f"vmag_signal_percent_{file.split('_')[2]}")







