import numpy as np
import pandas as pd

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


if __name__=='__main__':
    d_eu = pd.read_csv('data/exoplanet_catalogs/exoplanet.eu_catalog_hotjupiter_072621.csv', delimiter=',')
    translate_eu(d_eu)
    x = sync_catalogs()
    rslt_df = x[x['dec'] > -20]
    rslt_df = rslt_df.reset_index()
    rslt_df.to_csv('data/exoplanet_catalogs/exoplanets.csv')



