import os
os.sys.path.append('./utils/')
from utils.objects import load_object
from utils.load_inputs import fill_data
from utils.make_plots import MakePlots

if __name__ == '__main__':
    # load inputs
    configfile = 'calc_snr_max.cfg'
    so = load_object(configfile)
    cload = fill_data(so)

    make_plots = MakePlots(so)
    make_plots.plot_no_doppler(savefig=False)
    make_plots.plot_exo_doppler(savefig=False)
