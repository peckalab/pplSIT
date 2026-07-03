import json
import os
import fcntl

thread_budget = max(1, int(getattr(snakemake, "threads", 1)))
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"
os.environ["NUMBA_NUM_THREADS"] = "1"

import UnitMatchPy.extract_raw_data as erd
import numpy as np 
# from pathlib import Path
from joblib import Parallel, delayed
# import matplotlib.pyplot as plt


def infer_n_channels_total(dat_file, base_n_channels, dtype_bytes=2):
    n_bytes = os.path.getsize(dat_file)
    candidates = [base_n_channels, base_n_channels + 1]
    for n_channels_total in candidates:
        if n_bytes % (dtype_bytes * n_channels_total) == 0:
            return n_channels_total
    raise ValueError(
        f"Cannot infer channel count for {dat_file}. File size {n_bytes} is not "
        f"divisible by int16 samples for {base_n_channels} or {base_n_channels + 1} channels."
    )


def acquire_io_heavy_lock():
    lock_config = snakemake.config.get("io_heavy_lock", {})
    if not lock_config.get("enabled", False):
        return None

    lock_path = str(lock_config.get("path", "/tmp/pplSIT_io_heavy.lock"))
    os.makedirs(os.path.dirname(lock_path), exist_ok=True)
    lock_handle = open(lock_path, "w")
    print(f"Waiting for shared heavy-I/O lock: {lock_path}", flush=True)
    fcntl.flock(lock_handle, fcntl.LOCK_EX)
    print(f"Acquired shared heavy-I/O lock: {lock_path}", flush=True)
    return lock_handle

#Set Up Parameters
sample_amount = snakemake.config['unit_match']['sample_amount'] # 1000 # for both CV, at least 500 per CV
spike_width = snakemake.config['unit_match']['spike_width'] # 61 # assuming 30khz sampling, 82 and 61 are common choices (KS1-3, KS4), covers the AP and space around needed for processing
half_width = np.floor(spike_width/2).astype(int)
max_width = np.floor(spike_width/2).astype(int) #Size of area at start and end of recording to ignore to get only full spikes
n_channels = snakemake.config['unit_match']['n_channels'] #384 #neuropixels default, the number of channels EXCLUDING sync channels
extract_good_units_only = snakemake.config['unit_match']['extract_good_units_only'] #True # bool, set to true if you want to only extract units marked as good 

configured_jobs = int(os.environ.get(
    "RAW_WAVEFORM_WORKERS",
    snakemake.config.get('unit_match', {}).get('raw_waveform_workers', thread_budget)
))
n_jobs = max(1, min(configured_jobs, thread_budget))
print(f"Raw waveform extraction: using n_jobs={n_jobs} with Snakemake threads={thread_budget}")
io_heavy_lock_handle = acquire_io_heavy_lock()

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
    n_bytes = os.path.getsize(snakemake.input.dat_file)
    n_channels_tot = infer_n_channels_total(snakemake.input.dat_file, n_channels)
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
    n_bytes = os.path.getsize(snakemake.input.dat_file)
    n_channels_tot = infer_n_channels_total(snakemake.input.dat_file, n_channels)
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

# Record the extraction settings that determine the saved waveform shape.
with open(os.path.join(KS_dirs[0], 'RawWaveforms', 'RawWaveforms.ready'), 'w') as f:
    json.dump({
        "spike_width": int(spike_width),
        "samples_before": int(samples_before) if KS4_data else None,
        "sample_amount": int(sample_amount),
        "extract_good_units_only": bool(extract_good_units_only),
        "KS4_data": bool(KS4_data),
    }, f, indent=2)
    f.write("\n")
