Anaconda environment: py3


Purpose: 
-------
Calculate the SNR of H&K and KPF spectrometers give exposure time, stellar V magnitude, and stellar temperature

Usage:
----------
edit [var] parameters (T_eff (K), exp_time (s), vmag (generic Bessell V)) in calc_snr_max.cfg then run calc_snr_max.py to plot SNR of H&K and KPF spectrometers. Currently only certain T_eff values are allowed because this relies on which stellar Phoenix models are downloaded and placed in ./data/phoenix/

This code also has a function to generate  snr_master_hk.fits which is interpolated to match the KPF-ETC grids

Before running, will need to unzip telluric file, check paths

To Do:
---------
* Clean up plotting sequences, put into functions
* Incorporate this into the KPF-ETC interpolator
* Calculate errors on S_HK index and add saturation flag
* Adapt this code to be back end of online GUI ETC


Notes:
---------
This is a work in progress, dumping code here for future edits

Can check number of pixels KPF res element spans in following fits file: 
KPF_1000000rays_orders103-138_-_cal_111.11W_lfc_-_sci_500.00W_lfc_-_sky_500.00W_lfc_-_normalized_145_Green.FITS
