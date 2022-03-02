from utils.make_plots import MakePlots
from utils.load_inputs import fill_data
from utils.objects import load_object
import os
os.sys.path.append('./utils/')


if __name__ == '__main__':
    # load inputs
    configfile = 'hirax_snr_sodium.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    make_plots = MakePlots(so)
    # make_plots.plot_example_spectra(savefig=False)
    # make_plots.plot_hirax_profile(savefig=False)
    # make_plots.plot_no_doppler(savefig=False)
    # make_plots.plot_exo_doppler(savefig=False)
    make_plots.plot_stel_exo_doppler_one(savefig=False)
    # make_plots.plot_hirax_over_spectra(savefig=False)
    # make_plots.plot_hirax_over_shifted_exo(savefig=False)
    # make_plots.plot_lorentz_fit(savefig=False, feature=configfile.split('_')[2].split('.')[0])
    # make_plots.plot_noise_mag(cload, magnitudes=np.arange(7, 17, 1), savefig=False
