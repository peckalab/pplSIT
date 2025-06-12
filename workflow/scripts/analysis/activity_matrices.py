import h5py, os, sys, json
import numpy as np
import scipy.stats as scipystats


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import staple_spike_times


s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
animal  = session.split('_')[0]
#electrodes = snakemake.config['nMAP_electrodes'][animal]  # electrodes in A1

# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_ids = [x for x in f]
with h5py.File(snakemake.input[1], 'r') as f:
    for unit_name in unit_ids:
        spike_times[unit_name]  = np.array(f[unit_name]['spike_times'])


# matrix is always [units x time]
unit_mx_database = {}

# 10, 50, 250 ms binning
for bin_size in [0.01, 0.05, 0.25]:
    t_bins = sound_events[:, 0]

    bins_to_add = int(0.25/bin_size) - 1
    if bins_to_add > 0:
        res = t_bins.copy()
        for i in range(bins_to_add):
            res = np.vstack([res, t_bins + (i+1) * bin_size])
        t_bins = res.T.flatten()
        t_bins.sort()

    unit_mx = np.zeros([len(spike_times), len(t_bins)-1])
    for i, (unit_name, s_times) in enumerate(spike_times.items()):
        unit_mx[i], _ = np.histogram(s_times, bins=t_bins)

    unit_mx_database[int(bin_size) * 100] = {
        'mx': unit_mx.copy(),
        'mx_z': scipystats.zscore(unit_mx, axis=1),
        'bins': t_bins.copy(),
    }


# next two are interpolated to sound event binnning
t_bins = sound_events[:, 0]

# for 500ms
t_bins_even = t_bins[0::2]
t_bins_odd  = t_bins[1::2]

unit_mx = np.zeros([len(spike_times), len(t_bins)-1])
for i, (unit_name, s_times) in enumerate(spike_times.items()):
    spike_count_even, _ = np.histogram(s_times, bins=t_bins_even)
    spike_count_odd , _ = np.histogram(s_times, bins=t_bins_odd)

    indices_even = np.arange(0, unit_mx.shape[1], 2)[:len(spike_count_even)]
    indices_odd  = np.arange(0, unit_mx.shape[1], 2)[:len(spike_count_odd)]

    unit_mx[i][indices_even] = spike_count_even
    unit_mx[i][indices_odd]  = spike_count_odd

unit_mx_database[500] = {
    'mx': unit_mx.copy(),
    'mx_z': scipystats.zscore(unit_mx, axis=1),
    'bins': t_bins.copy(),
}

# for 1000ms
t_bins_1 = t_bins[0::4]
t_bins_2 = t_bins[1::4]
t_bins_3 = t_bins[2::4]
t_bins_4 = t_bins[3::4]

unit_mx = np.zeros([len(spike_times), len(t_bins)-1])
for i, (unit_name, s_times) in enumerate(spike_times.items()):
    for j, c_bins in enumerate([t_bins_1, t_bins_2, t_bins_3, t_bins_4]):
        spike_count, _ = np.histogram(s_times, bins=c_bins)
        indices  = np.arange(0, unit_mx.shape[1], 4)[:len(spike_count)]
        unit_mx[i][indices] = spike_count

unit_mx_database[1000] = {
    'mx': unit_mx.copy(),
    'mx_z': scipystats.zscore(unit_mx, axis=1),
    'bins': t_bins.copy(),
}

with h5py.File(snakemake.output[0], 'w') as f:
    for key, mx_dict in unit_mx_database.items():
        grp = f.create_group(f"mx_{key}ms")

        for name, ds in mx_dict.items():
            grp.create_dataset(name, data=ds)