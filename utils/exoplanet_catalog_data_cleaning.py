import numpy as np
import pandas as pd

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

def sync_catalogs(to_sync, sync_wth):
    """
    http://exoplanet.eu/catalog/
    and exoplanet archive have different formats

    load both and sync into one dictionary
    """
    d_eu = pd.read_csv(to_sync, delimiter=',')
    d_arxv = pd.read_csv(sync_wth, delimiter=',', header=212)

    d_eu_new = translate_eu(d_eu)
    df = pd.merge(d_eu_new, d_arxv, 'outer',
                  on=['pl_name', 'pl_massj', 'pl_orbper', 'pl_radj', 'pl_orbsmax', 'st_teff', 'st_mass', 'st_rad', 'ra',
                      'dec', 'sy_vmag'], copy=False, indicator=True)
    hj = df.drop_duplicates('pl_name')

    return hj.reset_index()


def filter_over_neg20_dec(cmplt_catalog_file):
    '''
    Reads in the eu catalog, translates all the columns to match with the arxv catalog and then syncs the two catalogs.
    Then the complete catalog is filtered to removed anything below a declination of -20 degrees and saved.
    :param catalog_eu: EU catalog
    :param catalog_arxv: ARXV catalog
    :param file_to_save: file name that is saved after filtering
    '''
    cmplt_catalog = pd.read_csv(cmplt_catalog_file)
    filtered_df = cmplt_catalog[cmplt_catalog['dec'] > -20].reset_index(drop=True)
    return filtered_df


def calc_absorbtion_signal(data):
    kb = 1.38064852e-23  # m2 kg s-2 K-1
    mu = 2.3 * 1.6605390e-27  # kg
    G = 6.67408e-11  # m3 kg-1 s-2

    # estimate Temp planet in kelvin
    Teq = (1 / 4.) ** (1 / 4.) * data['st_teff'] * np.sqrt(
        0.00465047 * data['st_rad'] / data['pl_orbsmax'])  # 0.00465047 AU per Rsun
    gravity = G * data['pl_massj'] * 1.898e27 / (data['pl_radj'] * 69911000.0) ** 2  # kg from m_jup

    # get H -> kb * T/ mu / g
    H = kb * Teq / mu / gravity  # meters

    # calculate A like Sing did
    A = 2 * data['pl_radj'] * 0.10049 * H * 1.4374e-9 / data['st_rad'] ** 2
    return A


def filter_catalog_no_nan_signal(data_file):
    data = pd.read_csv(data_file)
    signal = 2 * calc_absorbtion_signal(data)
    signal_percent = signal * 100
    data['signal_percent'] = signal_percent
    exoplanets_no_nan_signal = data[~np.isnan(data['signal_percent'])].reset_index(drop=True)
    return exoplanets_no_nan_signal


if __name__ == '__main__':

    eu_catalog = 'data/exoplanet_catalogs/exoplanet.eu_catalog_hotjupiter_072621.csv'
    arxv_catalog = 'data/exoplanet_catalogs/hotjupiters_northerhemisphere_072621.csv'

    cmplt_exoplanet_catalog = sync_catalogs(to_sync=eu_catalog, sync_wth=arxv_catalog)
    cmplt_exoplanet_catalog.to_csv('data/exoplanet_catalogs/complete_catalog.csv')

    complete_catalog = 'data/exoplanet_catalogs/complete_catalog.csv'
    exoplanets_over_neg20_dec = filter_over_neg20_dec(complete_catalog)
    exoplanets_over_neg20_dec.to_csv('data/exoplanet_catalogs/exoplanets_dec_over_20.csv')

    exoplanets_no_nan_signal = filter_catalog_no_nan_signal('data/exoplanet_catalogs/exoplanets_dec_over_20.csv')
    exoplanets_no_nan_signal.to_csv('data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv')

