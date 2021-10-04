import os
os.sys.path.append('./utils/')
from utils.objects import load_object
from utils.load_inputs import fill_data
from utils.make_plots import MakePlots
import numpy as np
import matplotlib

font = {'size': 12}
matplotlib.rc('font', **font)
matplotlib.style.use('seaborn-pastel')


if __name__ == '__main__':
    # load inputs
    configfile = 'hirax_snr_sodium.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    make_plots = MakePlots(so)
    # make_plots.plot_example_spectra(savefig=True)
    # make_plots.plot_hirax_profile(savefig=True)
    # make_plots.plot_no_doppler(savefig=True)
    # make_plots.plot_exo_doppler(savefig=True)
    # make_plots.plot_stel_exo_doppler(savefig=True)
    #cload.hirax(so,0.4)
    #make_plots.plot_hirax_over_spectra(savefig=True)
    # make_plots.plot_hirax_over_shifted_exo(savefig=True)
    # make_plots.plot_lorentz_fit(savefig=True, feature=configfile.split('_')[2].split('.')[0])
    # make_plots.plot_noise_mag(cload, magnitudes=np.arange(7, 17, 1), savefig=True)
    make_plots.plot_vmag_signal_percent_noise_tfac_ashley(config='config6',  sigma_thresh=3)


