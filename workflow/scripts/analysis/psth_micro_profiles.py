import os, sys
import h5py
import json
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import get_spike_counts
from utils.neurosuite import get_unit_names_sorted


# configuration
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
event_type_names = ['SIL', 'BGR', 'TGT', 'NOI']  # SIL, BGR, TGT, NOI - order matters
hw = snakemake.config['psth']['micro']['latency']
bc = snakemake.config['psth']['micro']['bin_count']
iter_count = snakemake.config['psth']['micro']['boot_iter_count']
bin_size   = hw/((bc-1)/2)


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

# first do bootstrapping
profiles = np.zeros([len(event_types), len(unit_names), iter_count, bc-1])
profile_stats = np.zeros([len(event_types), len(unit_names), 5, bc-1])

for j, event_id in enumerate(event_types):  # sound events
    # absolute indices to sound events
    idxs_pool = []  # some events here will be added twice to avoid bias from previous event at the moment of switch
    for k in range(len(sound_events)):
        # first event if special as there is no preceeding event
        if k == 0 and sound_events[k][1] == event_id:
            idxs_pool.append(k)
            
        elif sound_events[k][1] == event_id:
            if sound_events[k-1][1] != event_id:
                if not k+1 >= len(sound_events):
                    idxs_pool.append(k+1)  # add next event instead of the first one in a sequence to avoid bias
            else:
                idxs_pool.append(k)
    idxs_pool = np.array(idxs_pool)
    
    for u, unit_name in enumerate(unit_names):
        
        # bootstrap PSTHs
        for k in range(iter_count):
            # get bootstrapped events from the pool
            idxs_rand = np.random.choice(idxs_pool, len(idxs_pool), replace=True)
            times_rand = sound_events[idxs_rand][:, 0]
            
            # random jitter spike times for smoothing
            strain = spike_times[unit_name]
            strain = strain + ((np.random.rand(len(strain)) - 0.5) * bin_size)
                        
            bins, psth = get_spike_counts(strain, times_rand, hw=hw, bin_count=bc)
            profiles[j][u][k] = psth
            
        # compute stats right away
        confidence_low  = np.zeros(bc-1)
        confidence_high = np.zeros(bc-1)
        for k, col in enumerate(profiles[j][u].T):
            confidence_low[k]  = np.percentile(col, 5)
            confidence_high[k] = np.percentile(col, 95)

        profile_stats[j][u] = np.vstack([
            bins[:-1],  
            np.median(profiles[j][u], axis=0),
            profiles[j][u].std(axis=0),
            confidence_low, 
            confidence_high
        ])
        
# save to H5
with h5py.File(snakemake.output[0], 'w') as f:
    for i, event_name in enumerate(event_type_names):
        grp_ev = f.create_group(event_name)

        for j, unit_name in enumerate(unit_names):
            grp_unit = grp_ev.create_group(unit_name)
            grp_unit.create_dataset('profiles', data=profiles[i][j])
            grp_unit.create_dataset('profile_stats', data=profile_stats[i][j])
