import os, sys
import h5py
import json
import pywt
import numpy as np
from scipy import signal
#from scipy.ndimage import uniform_filter1d
#from scipy.ndimage import gaussian_filter1d

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.aeps import nanmean_filter1d, interp_nan_1d


area = snakemake.config['lfp']['aep_ITPC']['area']
fs = snakemake.config['lfp']['target_rate']
kernel_sizes = snakemake.config['lfp']['EV_SU']['kernel_sizes']

# AEPs
with h5py.File(snakemake.input[0], 'r') as f:
    cc_originals = {
        'avg_across_channels': np.array(f[area]['avg_across_channels']),
        'weighted_avg': np.array(f[area]['weighted_avg']),
        'zscored_avg': np.array(f[area]['zscored_avg']),
    }


cc_metrics = {}
for cc_orignal_name, cc_original in cc_originals.items():
    cc_m = {}
    for k_size in kernel_sizes:
        
        # 1. Smooth 0, 4, 8, 12 pulses - before metrics
        if k_size > 0:
            cc_smoothed = nanmean_filter1d(cc_original, size=k_size, axis=0)
        else:
            cc_smoothed = cc_original.copy()

        # 2. Compute EV / SU metrics
        EV = np.ptp(cc_smoothed[:, :125], axis=1)  # Evoked 0 - 125 ms
        SU = cc_smoothed[:, 180:250].mean(axis=1)  # Sustained 180 - 250 ms

        #mean_late = cc_smoothed[:, 0:250].mean(axis=1)  # Sustained 200 - 250 ms
        #mean_late = np.ptp(cc_smoothed[:, 170:250], axis=1)  # Sustained 130 - 250 ms

        sig = cc_smoothed[:, 0:125].std(axis=1)
        noise = cc_smoothed[:, 150:250].std(axis=1)
        #noise = np.abs(cc_smoothed[:, 0:250].mean(axis=1))
        
        cc_m[k_size] = np.stack([
            interp_nan_1d(EV), 
            interp_nan_1d(SU), 
            interp_nan_1d(sig), 
            interp_nan_1d(noise)
        ], axis=1)  # EV, SU, signal?, noise
    cc_metrics[cc_orignal_name] = dict(cc_m)

# write out
with h5py.File(snakemake.output[0], 'w') as f:
    for AEP_name, cc_m in cc_metrics.items():
        grp = f.create_group(AEP_name)
        for kernel_size, ds in cc_m.items():
            grp.create_dataset(str(kernel_size), data=ds)
