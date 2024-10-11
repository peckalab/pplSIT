import h5py, os, sys, json
import numpy as np
import itertools
from scipy import stats, signal
from sklearn import decomposition


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.population import unit_activity_matrix
from utils.neurosuite import get_unit_names_sorted
from utils.psth import staple_spike_times
from utils.spiketrain import get_shuffled


def make_smooth(data, width_in_bins):
    kernel = np.ones(width_in_bins)
    return np.convolve(data, kernel, 'same') / kernel.sum()


# read datasets
s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
animal  = session.split('_')[0]
electrodes = snakemake.config['nMAP_electrodes'][animal]  # electrodes in A1
win_sizes = snakemake.config['coherence']['win_sizes']
comp_count = snakemake.config['coherence']['comp_count']

meta_file = snakemake.input[0]
unit_file = snakemake.input[1]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(unit_file, 'r') as f:
    unit_names = get_unit_names_sorted([name for name in f if int(name.split('-')[0]) in electrodes])
    for unit_name in f:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

assert len(unit_names) >= comp_count

# periods of no sound stimulation
#periods = get_event_periods(tl, 0)  # old fashion

# no stimulus periods defined by fake sound pulses
SIL_times_beg = sound_events[sound_events[:, 1] == 0][:, 0]
SIL_times_end = SIL_times_beg + cfg['sound']['latency']
periods_SIL = np.vstack([SIL_times_beg, SIL_times_end]).T

# BGR periods - sustained part only: 0.120 to 0.25 for each pulse
period_dur = cfg['sound']['latency'] - snakemake.config['coherence']['sustained_start']
BGR_times_beg = sound_events[sound_events[:, 1] == 0][:, 0] + snakemake.config['coherence']['sustained_start']
BGR_times_end = BGR_times_beg + period_dur
periods_BGR = np.vstack([BGR_times_beg, BGR_times_end]).T

prefixes = ['SIL', 'BGR']
for i, periods in enumerate([periods_SIL, periods_BGR]):

    # just use equally-spaced binning for spiketrains
    bins = np.arange(0, np.diff(periods, axis=1).sum(), snakemake.config['coherence']['bin_size'])  # ignore last uneven bin
    unit_mx_real = np.zeros([len(unit_names), len(bins)-1])
    unit_mx_shuf = np.zeros([len(unit_names), len(bins)-1])

    for k, unit_name in enumerate(unit_names):
        # shrink all spikes as if there is no other periods.
        # allows to do shuffling without effects from periods with 
        # sounds where mean firing rates can be different
        s_times = spike_times[unit_name]
        strain = staple_spike_times(s_times, periods, mode='sequence')  # result is in periods!
        strain_real = np.array([item for sublist in strain for item in sublist])  # flatten to one array
        strain_shuf = get_shuffled(strain_real)

        unit_mx_real[k] = np.histogram(strain_real, bins=bins)[0]
        unit_mx_shuf[k] = np.histogram(strain_shuf, bins=bins)[0]

    # for pairwise correlations
    unit_pair_idxs = [x for x in itertools.combinations(range(len(unit_names)), 2)]
    corrs_real = np.zeros([len(win_sizes), len(unit_pair_idxs)])
    corrs_shuf = np.zeros([len(win_sizes), len(unit_pair_idxs)])
    var_real   = np.zeros([len(win_sizes), comp_count])
    var_shuf   = np.zeros([len(win_sizes), comp_count])
    comp_real  = np.zeros([len(win_sizes), comp_count, len(unit_names)])

    for w, win_size in enumerate(win_sizes):
        # z-score and smooth
        for j in range(len(unit_names)):
            unit_mx_real[j] = stats.zscore(unit_mx_real[j])
            unit_mx_real[j] = make_smooth(unit_mx_real[j], win_size)
            unit_mx_shuf[j] = stats.zscore(unit_mx_shuf[j])
            unit_mx_shuf[j] = make_smooth(unit_mx_shuf[j], win_size)

        # pairwise correlations
        for j, pair in enumerate(unit_pair_idxs):
            corr, pval = stats.pearsonr(unit_mx_real[pair[0]], unit_mx_real[pair[1]])
            corrs_real[w][j] = corr

            corr, pval = stats.pearsonr(unit_mx_shuf[pair[0]], unit_mx_shuf[pair[1]])
            corrs_shuf[w][j] = corr

        # PCA explained ratios and components
        pca_real = decomposition.PCA(n_components=comp_count)
        X_real   = pca_real.fit_transform(unit_mx_real.T)

        pca_shuf = decomposition.PCA(n_components=comp_count)
        X_shuf   = pca_shuf.fit_transform(unit_mx_shuf.T)

        var_real[w] = pca_real.explained_variance_ratio_
        var_shuf[w] = pca_shuf.explained_variance_ratio_

        # principal components
        comp_real[w] = pca_real.components_

    # dump to H5, overwrite
    prefix = prefixes[i]
    to_save = [corrs_real, corrs_shuf, var_real, var_shuf, comp_real]
    names = ['corrs_real', 'corrs_shuf', 'var_real', 'var_shuf', 'comp_real']
    with h5py.File(snakemake.output[0], 'a') as f:
        for n, name in enumerate(names):
            name_pfx = name + '_' + prefix
            if name_pfx in f:
                del f[name_pfx]

            f.create_dataset(name_pfx, data=to_save[n])

            