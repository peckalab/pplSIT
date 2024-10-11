import os, sys
import h5py
import json
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import compute_shuffled_metrics, staple_pulsetrain, staple_spike_times, get_spike_counts
from utils.neurosuite import get_unit_names_sorted
from utils.events import get_sound_event_periods


# configuration
hw = snakemake.config['psth']['micro']['latency']
bc = snakemake.config['psth']['micro']['bin_count']
iter_count = snakemake.config['psth']['micro']['boot_iter_count']
bin_size   = hw/((bc-1)/2)
split_ratio = snakemake.config['metronome']['split_ratio']


# reading events and spiking data
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in f:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])


# 1 - split SIL condition into profile and metronome parts
event_type = 0  # no stimulus
t_periods = np.array(get_sound_event_periods(sound_events, event_type))  # these are absolute times

prof_count = int(len(t_periods)*split_ratio)
test_count = len(t_periods) - prof_count
prof_idxs = np.random.choice(np.arange(len(t_periods)), prof_count, replace=False)
test_idxs = np.array([x for x in range(len(t_periods)) if x not in prof_idxs])
prof_idxs.sort()
test_idxs.sort()

# save periods and profile / test indices
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('periods', data=t_periods)
    f.create_dataset('prof_idxs', data=prof_idxs)
    f.create_dataset('test_idxs', data=test_idxs)
print('Periods created')


# 2 - create PSTH profiles for all units based on profile periods
idxs_sil_prof_events = []
for p in t_periods[prof_idxs]:
    idxs_period = np.where( (sound_events[:, 0] > p[0]-0.00001) & (sound_events[:, 0] < p[1]+0.00001) )[0]
    idxs_period[0] = idxs_period[1]  # replace first response with the second to avoid previous stimulus effect
    idxs_sil_prof_events.append(idxs_period)
idxs_sil_prof_events = np.array([item for sublist in idxs_sil_prof_events for item in sublist])  # flatten

# do bootstrapping
profiles = np.zeros([len(unit_names), iter_count, bc-1])
profile_stats = np.zeros([len(unit_names), 5, bc-1])

for u, unit_name in enumerate(unit_names):
    
    # bootstrap PSTHs
    for k in range(iter_count):
        # get bootstrapped events from the pool
        idxs_rand = np.random.choice(idxs_sil_prof_events, len(idxs_sil_prof_events), replace=True)
        times_rand = sound_events[idxs_rand][:, 0]
        
        # random jitter spike times for smoothing
        strain = spike_times[unit_name]
        strain = strain + ((np.random.rand(len(strain)) - 0.5) * bin_size)
                    
        bins, psth = get_spike_counts(strain, times_rand, hw=hw, bin_count=bc)
        profiles[u][k] = psth
        
    # compute stats right away
    confidence_low  = np.zeros(bc-1)
    confidence_high = np.zeros(bc-1)
    for k, col in enumerate(profiles[u].T):
        confidence_low[k]  = np.percentile(col, 5)
        confidence_high[k] = np.percentile(col, 95)

    profile_stats[u] = np.vstack([
        bins[:-1],  
        np.median(profiles[u], axis=0),
        profiles[u].std(axis=0),
        confidence_low, 
        confidence_high
    ])
        
# save to H5
with h5py.File(snakemake.output[1], 'w') as f:
    for j, unit_name in enumerate(unit_names):
        grp_unit = f.create_group(unit_name)
        grp_unit.create_dataset('profile_stats', data=profile_stats[j])
print('PSTH complete')


# 3 - compute shuffled PSTHs
shuffled = np.zeros([len(unit_names), 5, bc-1])
pulses   = sound_events[sound_events[:, 1] == 0][:, 0]

# shrink all selected pulses as if there is no other periods.
# allows to do shuffling without effects from other periods
# where mean firing rates can be different
adjusted_pulses = staple_pulsetrain(pulses, t_periods[prof_idxs])  # important - here only profile periods

for k, unit_name in enumerate(unit_names):
    s_times = spike_times[unit_name]
    strain = staple_spike_times(s_times, t_periods[prof_idxs], mode='sequence')  # result is in periods!
    strain = np.array([item for sublist in strain for item in sublist])  # flatten to one array
    shuffled[k] = compute_shuffled_metrics(strain, adjusted_pulses, offset=hw, bin_count=bc, iter_count=iter_count)

# save to H5
with h5py.File(snakemake.output[2], 'w') as f:
    for j, unit_name in enumerate(unit_names):
        grp_unit = f.create_group(unit_name)
        grp_unit.create_dataset('shuffled', data=shuffled[j])
print('Shuffle complete')