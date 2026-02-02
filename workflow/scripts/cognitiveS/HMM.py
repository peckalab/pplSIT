import os, sys
import h5py
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.cognitiveS.HMM import *


#cfg = snakemake.config['trajectories']

meta_file = snakemake.input[0]
actm_file = snakemake.input[1]
covm_file = snakemake.input[2]

s_path  = os.path.dirname(meta_file)
session = os.path.basename(s_path)
animal  = session.split('_')[0]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

# 50ms binning
with h5py.File(actm_file, 'r') as f:
    X_counts_50ms = np.array(f['mx_50ms']['mx'])
    t_edges_50ms  = np.array(f['mx_50ms']['bins'])
    
# same 50ms binning but 0.5sec smoothed covariance
grp_name = 'mx_50ms_mx_sm10'
with h5py.File(covm_file, 'r') as f:
    # take labels / sorting from the 'all' condition
    labels    = np.array(f[grp_name]['all']['labels'])
    idxs_sort = np.array(f[grp_name]['all']['idxs_sort'])


def get_filtered(labels, to_keep):
    labels_filt = []
    for l in labels:
        labels_filt.append(l if l in to_keep else 0)
    return np.array(labels_filt)

# manually select clusters
state_candidates = {
    '008229_hippoSIT_2022-05-16_20-36-44': [1, 3, 4, 6, 8],
    '008229_hippoSIT_2022-05-17_21-44-43': [1, 2, 6],
    '008229_hippoSIT_2022-05-18_14-36-18': [1, 2, 3, 5, 7],
    '008229_hippoSIT_2022-05-20_15-54-39': [1, 2, 3, 5, 8]
}
to_keep = state_candidates[session]

results = run_sticky_hmm_pipeline(
    spike_counts=X_counts_50ms,        # neurons x time
    ensemble_labels=get_filtered(labels, to_keep),        # neurons
    smooth_bins=10,                 # 250 ms
    drift_bins=2400,                # 800 = 40 s
    K=4,
    stickiness=15.0,
    output_h5=snakemake.output[0],
)