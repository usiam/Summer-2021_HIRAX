Anaconda environment: py3


Purpose: 
-------
Optimizing bandpass configuration for a novel instrument - HIRAX - to allow observations of exoplanet (Hot Jupiter) atmosphere from the Earth.

Usage:
----------
edit [var] parameters (T_eff (K), exp_time (s), vmag (generic Bessell V)) in calc_snr_max.cfg then run calc_snr_max.py to plot SNR of H&K and KPF spectrometers. Currently only certain T_eff values are allowed because this relies on which stellar Phoenix models are downloaded and placed in ./data/phoenix/

This code also has a function to generate  snr_master_hk.fits which is interpolated to match the KPF-ETC grids

Before running, will need to unzip telluric file, check paths.

To Do:
---------
* Clean up plotting sequences, put into functions
* Incorporate this into the KPF-ETC interpolator
* Calculate errors on S_HK index and add saturation flag
* Adapt this code to be back end of online GUI ETC


Notes:
---------
Apply velocity shifts on load
Test new bandpass configurations out
