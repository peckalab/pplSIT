import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.decoder_state import *


#cfg = snakemake.config['trajectories']

STATE_KEYS = [
    "idxs_tgt_sta_succ",
    "idxs_bgr_sta",
    "idxs_sil_sta",
]

# -----------------------------
# Reading datasets
# -----------------------------
meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
actm_file = snakemake.input[2]
traj_file = snakemake.input[3]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

# fixed episode onset / length (22 bins)
tgt_succ_mx = tgt_mx[tgt_mx[:, 4] == 1]
targets_periods_bins = np.column_stack([tgt_succ_mx[:, 0] + 2, tgt_succ_mx[:, 0] + 24])

# 250 ms binning to match PCA representation
with h5py.File(actm_file, 'r') as f:
    X_counts = np.array(f['mx_250ms']['mx'])
    t_edges  = np.array(f['mx_250ms']['bins'])
    
#X_counts_50ms = gaussian_filter1d(X_counts_50ms, sigma=5, axis=0, mode="nearest")

with h5py.File(segm_file, 'r') as f:
    state_idxs = {}
    for idxs_name in STATE_KEYS:
        state_idxs[idxs_name] = np.array(f[idxs_name]).astype(np.int32)

