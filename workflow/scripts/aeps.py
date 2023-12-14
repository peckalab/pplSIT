import os
import h5py
import json
import numpy as np


aep_dur = snakemake.config['lfp']['aep_dur']  # AEP duration in sec
s_rate  = snakemake.config['lfp']['s_rate']   # LFP samp rate Hz
dur = int(aep_dur*s_rate)

with h5py.File(snakemake.input[1], 'r') as f:
    lfp_mx = np.array(f['lfp'])

with h5py.File(snakemake.input[2], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])

with open(snakemake.input[0]) as json_file:
    channels = json.load(json_file)['AEPs']

for area, channel in channels.items():
    lfp = lfp_mx[:, channel]

    # constructing AEP matrix
    aeps = np.zeros([len(sound_events), dur])
    for i, event in enumerate(sound_events):
        idx_event = int(event[0]*s_rate)
        aep = lfp[idx_event:idx_event + dur]
        if len(aep) == dur:
            aeps[i] = aep

    # save AEPs to file
    aeps_name = 'aeps_%d_%d' % (channel, dur)
    with h5py.File(snakemake.output[0], 'a') as f:
        if not area in f:
            f.create_group(area)
        if aeps_name in f[area]:
            del f[area][aeps_name]
        f[area].create_dataset(aeps_name, data=aeps)