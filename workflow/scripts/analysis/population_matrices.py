import h5py, os, sys, json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import staple_spike_times


s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
animal  = session.split('_')[0]
electrodes = snakemake.config['nMAP_electrodes'][animal]  # electrodes in A1

# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = [x for x in f if int(x.split('-')[0]) in electrodes]
with h5py.File(snakemake.input[1], 'r') as f:
    for unit_name in unit_names:
        spike_times[unit_name]  = np.array(f[unit_name]['spike_times'])

# compute population responses
bin_size = snakemake.config['nMAP_EV_SU']['bin_size']
ev_su_lag = snakemake.config['nMAP_EV_SU']['ev_su_lag']

ev_bin_count = int(ev_su_lag/bin_size)
ev_times_beg = sound_events[:, 0]
ev_times_end = ev_times_beg + ev_su_lag
ev_periods = np.vstack([ev_times_beg, ev_times_end]).T

su_bin_count = int((cfg['sound']['latency'] - ev_su_lag)/bin_size)
su_times_beg = sound_events[:, 0] + ev_su_lag
su_times_end = su_times_beg + cfg['sound']['latency'] - ev_su_lag
su_periods = np.vstack([su_times_beg, su_times_end]).T

ev_bins = np.arange(0, np.diff(ev_periods, axis=1).sum(), bin_size)  # ignore last uneven bin
ev_unit_mx = np.zeros([len(unit_names), len(ev_bins)-1])
su_bins = np.arange(0, np.diff(su_periods, axis=1).sum(), bin_size)  # ignore last uneven bin
su_unit_mx = np.zeros([len(unit_names), len(su_bins)-1])

for k, unit_name in enumerate(unit_names):
    s_times = spike_times[unit_name]

    # shrink all spikes as if there is no other periods.
    ev_strain = staple_spike_times(s_times, ev_periods, mode='sequence')  # result is in periods!
    ev_strain = np.array([item for sublist in ev_strain for item in sublist])  # flatten to one array
    ev_unit_mx[k] = np.histogram(ev_strain, bins=ev_bins)[0]

    su_strain = staple_spike_times(s_times, su_periods, mode='sequence')  # result is in periods!
    su_strain = np.array([item for sublist in su_strain for item in sublist])  # flatten to one array
    su_unit_mx[k] = np.histogram(su_strain, bins=su_bins)[0]

# finally dump everything
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('ev_bin_count', data=ev_bin_count)
    f.create_dataset('ev_periods', data=ev_periods)
    f.create_dataset('ev_bins', data=ev_bins)
    f.create_dataset('ev_unit_mx', data=ev_unit_mx)  # original activity matrix, 10 ms bins

    f.create_dataset('su_bin_count', data=su_bin_count)
    f.create_dataset('su_periods', data=su_periods)
    f.create_dataset('su_bins', data=su_bins)
    f.create_dataset('su_unit_mx', data=su_unit_mx)  # original activity matrix, 10 ms bins