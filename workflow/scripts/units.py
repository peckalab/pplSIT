import os, sys
import h5py
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import load_clu_res, XMLHero
from utils.spiketrain import instantaneous_rate, spike_idxs
from utils.hdf import create_dataset, H5NAMES


session_path = os.path.dirname(snakemake.input[0])

# loading unit data
units = load_clu_res(session_path)  # spikes are in samples, not seconds
sampling_rate = XMLHero(snakemake.input[0]).get_sampling_rate()

# loading timeline
with h5py.File(snakemake.input[1], 'r') as f:
    tl = np.array(f['processed']['timeline'])  # time, X, Y, speed, HD, trials, sounds
    
for electrode_idx in units.keys():
    unit_idxs = units[electrode_idx]

    for unit_idx, spiketrain in unit_idxs.items():
        unit_name = '%s-%s' % (electrode_idx, unit_idx)

        # 3 main ways to store a spiketrain
        s_times = spiketrain/sampling_rate  # spike times in seconds
        i_rate  = instantaneous_rate(s_times, tl[:, 0])  # instantaneous rate, sampling to timeline
        s_idxs  = spike_idxs(s_times, tl[:, 0])  # timeline indices of individual spikes

        create_dataset(snakemake.output[0], unit_name, H5NAMES.spike_times, s_times)
        create_dataset(snakemake.output[0], unit_name, H5NAMES.inst_rate, i_rate)
        create_dataset(snakemake.output[0], unit_name, H5NAMES.spike_idxs, s_idxs)