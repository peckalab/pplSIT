import os
import h5py
import json
import numpy as np
from scipy import signal


aep_dur = snakemake.config['lfp']['aep_dur']  # AEP duration in sec
s_rate  = snakemake.config['lfp']['s_rate']   # LFP samp rate Hz
dur = int(aep_dur*s_rate)
nyquist = 0.5 * s_rate

with h5py.File(snakemake.input[1], 'r') as f:
    lfp_mx = np.array(f['lfp'])

with h5py.File(snakemake.input[2], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])

with open(snakemake.input[0]) as json_file:
    channels = json.load(json_file)['AEPs']

for area, channel in channels.items():
    lfp = lfp_mx[:, channel]  # raw LFP

    # first cut out outliers (is this correct?)
    lfp[lfp >  4*lfp.std()] =  4*lfp.std()
    lfp[lfp < -4*lfp.std()] = -4*lfp.std()

    # try to filter out 50Hz
    low  = 47 / nyquist
    high = 53 / nyquist
    sos = signal.butter(10, [low, high], btype='bandstop', output='sos')
    lfp_filt = signal.sosfiltfilt(sos, lfp)

    # try to filter < 1Hz slow staff
    h_cut = 1 / nyquist
    sos = signal.butter(10, h_cut, btype='highpass', output='sos')
    lfp_filt = signal.sosfiltfilt(sos, lfp_filt)

    # constructing AEP matrix
    aeps   = np.zeros([len(sound_events), dur])
    aeps_f = np.zeros([len(sound_events), dur])
    for i, event in enumerate(sound_events):
        idx_event = int(event[0]*s_rate)
        aep   = lfp[idx_event:idx_event + dur]
        aep_f = lfp_filt[idx_event:idx_event + dur]
        if len(aep) == dur:
            aeps[i] = aep
            aeps_f[i] = aep_f

    # save AEPs to file
    aeps_name   = 'aeps_%d_%d' % (channel, dur)
    aeps_name_f = 'aeps_%d_%d_filt' % (channel, dur)
    with h5py.File(snakemake.output[0], 'a') as f:
        if not area in f:
            f.create_group(area)
        if aeps_name in f[area]:
            del f[area][aeps_name]
        if aeps_name_f in f[area]:
            del f[area][aeps_name_f]
        f[area].create_dataset(aeps_name, data=aeps)
        f[area].create_dataset(aeps_name_f, data=aeps_f)