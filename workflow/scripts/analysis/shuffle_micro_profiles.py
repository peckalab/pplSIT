import os, sys
import h5py
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import compute_shuffled_metrics, staple_pulsetrain, staple_spike_times
from utils.events import get_event_periods
from utils.neurosuite import get_unit_names_sorted

# some configs
hw = snakemake.config['psth']['micro']['latency']
bc = snakemake.config['psth']['micro']['bin_count']
iter_count = snakemake.config['psth']['micro']['shuf_iter_count']
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
event_type_names = ['SIL', 'BGR', 'TGT', 'NOI']  # SIL, BGR, TGT, NOI - order matters

# reading timeline and spike times
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    
spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in f:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

# compute shuffled PSTHs
shuffled = np.zeros([len(event_types), len(unit_names), 5, bc-1])
for j, event_id in enumerate(event_types):  # sound events
    pulses  = sound_events[sound_events[:, 1] == event_id][:, 0]
    periods = get_event_periods(tl, event_id)

    # shrink all selected pulses as if there is no other periods.
    # allows to do shuffling without effects from other periods
    # where mean firing rates can be different
    adjusted_pulses = staple_pulsetrain(pulses, periods)
    
    for k, unit_name in enumerate(unit_names):
        s_times = spike_times[unit_name]
        strain = staple_spike_times(s_times, periods, mode='sequence')  # result is in periods!
        strain = np.array([item for sublist in strain for item in sublist])  # flatten to one array
        shuffled[j][k] = compute_shuffled_metrics(strain, adjusted_pulses, offset=hw, bin_count=bc, iter_count=iter_count)

# save to H5
with h5py.File(snakemake.output[0], 'w') as f:
    for i, event_name in enumerate(event_type_names):
        grp_ev = f.create_group(event_name)

        for j, unit_name in enumerate(unit_names):
            grp_unit = grp_ev.create_group(unit_name)
            grp_unit.create_dataset('shuffled', data=shuffled[i][j])




# def compute_shuffle_psth_micro(meta_file, units_file, dst_unit_file, event_type, latency, bin_count):
#     unit_name = os.path.basename(dst_unit_file).split('_')[0]
#     with h5py.File(units_file, 'r') as f:
#         s_times = np.array(f[unit_name]['spike_times'])

#     with h5py.File(meta_file, 'r') as f:
#         tl = np.array(f['processed']['timeline'])
#         sound_events = np.array(f['processed']['sound_events'])

#     pulses  = sound_events[sound_events[:, 1] == event_type][:, 0]
#     periods = get_event_periods(tl, event_type)

#     # first shrink all selected pulses as if there is no other periods
#     adjusted_pulses = staple_pulsetrain(pulses, periods)

#     # spikes during certain event type only
#     strain = staple_spike_times(s_times, periods, mode='sequence')  # result is in periods!
#     strain = np.array([item for sublist in strain for item in sublist])  # flatten to one array
#     data = compute_shuffled_metrics(strain, adjusted_pulses, offset=latency, bin_count=bin_count)

#     np.save(dst_unit_file, data)


# # compute for silence
# compute_shuffle_psth_micro(
#     snakemake.input[0], 
#     snakemake.input[1], 
#     snakemake.output[0],
#     0,  # silence event type
#     snakemake.config['psth']['micro']['latency'],
#     snakemake.config['psth']['micro']['bin_count']
# )
