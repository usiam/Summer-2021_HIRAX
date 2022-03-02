from matplotlib import pyplot as plt
import matplotlib
from scipy import interpolate
import numpy as np
import pandas as pd

matplotlib.style.use('seaborn-pastel')


def plot_vmag_signal_percent_noise_tfac(
        exoplanet_data_file_name='data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv',
        config='config1', tfactor=np.arange(0.2, 1, 0.2), nsigma=3, savefig=False):

    # plots the exoplanets
    exoplanet_data = pd.read_csv(exoplanet_data_file_name)
    vmag = exoplanet_data['sy_vmag']
    signal_percent = exoplanet_data['signal_percent']
    labels = exoplanet_data['pl_name']

    fig, ax = plt.subplots(figsize=(17, 9))
    ax.semilogx(signal_percent, vmag, '+k', alpha=1)
    ax.set_ylim((16, 7))
    ax.set_xlim(min(signal_percent), 1)
    ax.set_xlabel('Transit signal of 2 Atmospheric Scale Height (%)')
    ax.set_ylabel('V magnitude')
    ax.set_xticks([0.01, 0.025, 0.050, 0.075, 0.1], minor=True)

    # writes the exoplanet names
    for i, txt in enumerate(labels):
        ax.annotate(txt, (signal_percent[i], vmag[i]), fontsize=5, alpha=0.5)

    # plots the nsigma noise lines for the throughput factors
    path = 'data/output/'
    file = 'amplitude_vmag_'
    extension = '.csv'
    model_noise_data = pd.read_csv(path + file + config + extension)
    for j, tfac in enumerate(tfactor):
        model_3sig_noise_percent = nsigma * \
            model_noise_data[f"amplitude_err_tfac{j + 1}"] * 100
        if j % 2:
            marker = '-'
        else:
            marker = '--'
        ax.semilogx(model_3sig_noise_percent, model_noise_data['magnitude'], marker,
                    label=f'Throughput={round(tfac, 3)}', alpha=1)
    ax.legend()
    if config == 'potassium':
        cfg = 'potassium'
    else:
        cfg = config[-1]
    fig.suptitle(f'Configuration {cfg}')
    plt.tight_layout()
    if savefig:
        fig.savefig(f"figures/vmag_signal_percent_{config}.png", dpi=600)


def find_detectable_exoplanet_counts(tfactor_ind, exoplanets_data, model_noise_data):
    signal = exoplanets_data['signal_percent']
    # add the signal percent to detectables csv
    model_3sig_percent = 3 * \
        model_noise_data[f'amplitude_err_tfac{tfactor_ind + 1}'] * 100

    tck = interpolate.splrep(
        model_3sig_percent, model_noise_data['magnitude'], k=3, s=0)
    exoplanets_in_limit = exoplanets_data[(signal > min(model_3sig_percent)) & (signal < max(model_3sig_percent))
                                          & (exoplanets_data['sy_vmag'] > min(model_noise_data['magnitude']))].reset_index(drop=True)
    target_mags = interpolate.splev(
        exoplanets_in_limit['signal_percent'], tck, der=0)
    detectables = exoplanets_in_limit[(
        exoplanets_in_limit['sy_vmag'] < target_mags)].reset_index(drop=True)

    return len(detectables)


def display_counts(counts, tfactor, savefig=False):
    throughput_labels = [f'Throughput={round(tfac,2)}' for tfac in tfactor]
    config_labels = [f'C{i + 1}' for i in range(len(counts))]

    data_lists = []
    for i in range(len(config_labels)):
        group = [config_labels[i]] + counts[i]
        data_lists.append(group)

    df = pd.DataFrame(data_lists, columns=[
                      'Configurations'] + throughput_labels)

    ax = df.plot(x='Configurations',
                 kind='bar',
                 width=0.8,
                 stacked=False,
                 figsize=(10, 8),
                 title='Exoplanet counts')
    for container in ax.containers:
        ax.bar_label(container)
    ax.set_ylabel('Counts')
    ax.legend(loc='upper left')
    plt.tight_layout()
    plt.show()
    if savefig:
        plt.savefig('figures/exoplanet_counts.png', dpi=600)


def detectables_data(model_noise_data, config_name):
    exoplanets_data = pd.read_csv(
        'Summer-2021_HIRAX/data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv')
    signal = exoplanets_data['signal_percent']
    # add the signal percent to detectables csv
    model_3sig_percent = 3 * \
        model_noise_data[f'amplitude_err_tfac3'] * 100
    tck = interpolate.splrep(
        model_3sig_percent, model_noise_data['magnitude'], k=3, s=0)
    exoplanets_in_limit = exoplanets_data[(signal > min(model_3sig_percent)) & (signal < max(model_3sig_percent))
                                          & (exoplanets_data['sy_vmag'] > min(model_noise_data['magnitude']))].reset_index(drop=True)
    target_mags = interpolate.splev(
        exoplanets_in_limit['signal_percent'], tck, der=0)
    detectables = exoplanets_in_limit[(
        exoplanets_in_limit['sy_vmag'] < target_mags)].reset_index(drop=True)

    detectables.to_csv(f'{config_name}.csv')


if __name__ == '__main__':
    num_configs = 5
    # for i in range(num_configs):
    #     plot_vmag_signal_percent_noise_tfac(
    #         config=f'config{i+1}', savefig=True)

    exoplanets = pd.read_csv(
        'Summer-2021_HIRAX/data/exoplanet_catalogs/exoplanets_dec_over_20_no_nan_signal.csv')
    tfactor = [0.5]
    counts = [[] for i in range(num_configs)]
    for i in range(len(counts)):
        path = 'Summer-2021_HIRAX/data/output/'
        file = f'amplitude_vmag_config{i+1}'
        extension = '.csv'
        model_noise_data = pd.read_csv(path + file + extension)
        detectables_data(model_noise_data=model_noise_data,
                         config_name=f'cofig{i+1}')
        # for j in range(len(tfactor)):
        #     counts[i].append(find_detectable_exoplanet_counts(
        #         tfactor_ind=j, exoplanets_data=exoplanets, model_noise_data=model_noise_data, config_name=f'cofig{i+1}'))
    # display_counts(counts, tfactor, savefig=False)

    # need to fix it for potassium
    # path = 'data/output/'
    # file = f'amplitude_vmag_potassium'
    # extension = '.csv'
    # model_noise_data = pd.read_csv(path + file + extension)
    # counts_potassium = []
    # for j in range(len(tfactor)):
    #     counts_potassium.append(find_detectable_exoplanet_counts(tfactor_ind=j, exoplanets_data=exoplanets,
    #                                                       model_noise_data=model_noise_data))
    # display_counts(counts_potassium, tfactor)
