import UnitMatchPy.extract_raw_data as erd
import numpy as np 
# from pathlib import Path
from joblib import Parallel, delayed
# import matplotlib.pyplot as plt
import os
import json

#Set Up Parameters
sample_amount = snakemake.config['unit_match']['sample_amount'] # 1000 # for both CV, at least 500 per CV
spike_width = snakemake.config['unit_match']['spike_width'] # 61 # assuming 30khz sampling, 82 and 61 are common choices (KS1-3, KS4), covers the AP and space around needed for processing
half_width = np.floor(spike_width/2).astype(int)
max_width = np.floor(spike_width/2).astype(int) #Size of area at start and end of recording to ignore to get only full spikes
n_channels = snakemake.config['unit_match']['n_channels'] #384 #neuropixels default, the number of channels EXCLUDING sync channels
extract_good_units_only = snakemake.config['unit_match']['extract_good_units_only'] #True # bool, set to true if you want to only extract units marked as good 

n_jobs = -1

KS4_data = snakemake.config['unit_match']['KS4_data'] #True #bool, set to true if using Kilosort, as KS4 spike times refer to start of waveform not peak
if KS4_data:
    samples_before = snakemake.config['unit_match']['samples_before'] #20
    samples_after = spike_width - samples_before
    max_width = samples_after #Number of samples on either side of the 

# KS_dirs is the directory containing the spk_times file
KS_dirs = [os.path.dirname(snakemake.input.spk_times)]
# n_sessions = len(KS_dirs) #How many session are being extracted
spike_ids, spike_times, good_units, all_unit_ids = erd.extract_KS_data(KS_dirs, extract_good_units_only = extract_good_units_only)
# extract_KS_data returns lists of length n_sessions, but here we have only one session
spike_ids = spike_ids[0]
spike_times = spike_times[0]
good_units = good_units[0]
all_unit_ids = all_unit_ids[0]

#Extract the units 

if extract_good_units_only:
    #load metadata
    with open(snakemake.input.oebin_file, 'r') as file:
        meta = json.load(file)
    n_bytes = os.path.getsize(snakemake.input.dat_file)
    n_channels_tot = int(meta['continuous'][0]['num_channels'])
    n_samples = int(n_bytes / (2*n_channels_tot))

    #create memmap to raw data, for that session
    data = np.memmap(snakemake.input.dat_file, dtype = 'int16', shape =(n_samples, n_channels_tot))

    # Remove spike which won't have a full waveform recorded
    spike_ids_tmp = np.delete(spike_ids, np.logical_or( (spike_times < max_width), ( spike_times > (data.shape[0] - max_width))))
    spike_times_tmp = np.delete(spike_times, np.logical_or( (spike_times < max_width), ( spike_times > (data.shape[0] - max_width))))


    #might be slow extracting sample for good units only?
    sample_idx = erd.get_sample_idx(spike_times_tmp, spike_ids_tmp, sample_amount, units = good_units)

    if KS4_data:
        avg_waveforms = Parallel(n_jobs = n_jobs, verbose = 10, mmap_mode='r', max_nbytes=None )(delayed(erd.extract_a_unit_KS4)(sample_idx[uid], data, samples_before, samples_after, spike_width, n_channels, sample_amount)for uid in range(good_units.shape[0]))
        avg_waveforms = np.asarray(avg_waveforms)           
    else:
        avg_waveforms = Parallel(n_jobs = n_jobs, verbose = 10, mmap_mode='r', max_nbytes=None )(delayed(erd.extract_a_unit)(sample_idx[uid], data, half_width, spike_width, n_channels, sample_amount)for uid in range(good_units.shape[0]))
        avg_waveforms = np.asarray(avg_waveforms)

    #Save in file named 'RawWaveforms' in the KS Directory
    erd.save_avg_waveforms(avg_waveforms, KS_dirs[0], all_unit_ids, good_units = good_units, extract_good_units_only = extract_good_units_only)
else:
    
    #Extracting ALL the Units
    n_units = len(np.unique(spike_ids))
    #load metadata
    with open(snakemake.input.oebin_file, 'r') as file:
        meta = json.load(file)
    n_bytes = os.path.getsize(snakemake.input.dat_file)
    n_channels_tot = int(meta['continuous'][0]['num_channels'])
    n_samples = int(n_bytes / (2*n_channels_tot))

    #create memmap to raw data, for that session
    data = np.memmap(snakemake.input.dat_file, dtype = 'int16', shape =(n_samples, n_channels_tot))

    # Remove spikes which won't have a full waveform recorded
    spike_ids_tmp = np.delete(spike_ids, np.logical_or( (spike_times < max_width), ( spike_times > (data.shape[0] - max_width))))
    spike_times_tmp = np.delete(spike_times, np.logical_or( (spike_times < max_width), ( spike_times > (data.shape[0] - max_width))))


    sample_idx = erd.get_sample_idx(spike_times_tmp, spike_ids_tmp, sample_amount, units= np.unique(spike_ids))
    
    if KS4_data:
        avg_waveforms = Parallel(n_jobs = n_jobs, verbose = 10, mmap_mode='r', max_nbytes=None )(delayed(erd.extract_a_unit_KS4)(sample_idx[uid], data, samples_before, samples_after, spike_width, n_channels, sample_amount)for uid in range(n_units))
        avg_waveforms = np.asarray(avg_waveforms)           
    else:
        avg_waveforms = Parallel(n_jobs = n_jobs, verbose = 10, mmap_mode='r', max_nbytes=None )(delayed(erd.extract_a_unit)(sample_idx[uid], data, half_width, spike_width, n_channels, sample_amount)for uid in range(n_units))
        avg_waveforms = np.asarray(avg_waveforms)

    #Save in file named 'RawWaveforms' in the KS Directory
    erd.save_avg_waveforms(avg_waveforms, KS_dirs[0], all_unit_ids, good_units = good_units, extract_good_units_only = extract_good_units_only)

del data

# touch RawWaveforms.ready to signal completion
with open(os.path.join(KS_dirs[0], 'RawWaveforms', 'RawWaveforms.ready'), 'w') as f:
    f.write('')