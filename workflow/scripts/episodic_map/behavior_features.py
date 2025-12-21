import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.behavior_features import *

#cfg = snakemake.config['trajectories']

# -----------------------------
# Reading datasets
# -----------------------------
meta_file = snakemake.input[0]
dlc_file  = snakemake.input[1]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

# fixed episode onset / length (22 sound bins)
tgt_succ_mx = tgt_mx[tgt_mx[:, 4] == 1]
targets_periods_tl_bins = np.column_stack([
    tgt_succ_mx[:, 2] + 2*25, 
    tgt_succ_mx[:, 2] + 24*25
]) 
targets_periods_times   = np.column_stack([
    tl[targets_periods_tl_bins[:, 0]][:, 0],
    tl[targets_periods_tl_bins[:, 1]][:, 0]
])

with h5py.File(dlc_file, 'r') as f:
    dlc_mat = np.array(f['df_DLC'])
    dlc_columns = f['df_DLC'].attrs['headers'].tolist()

# -----------------------------
# Compute behavior features
# -----------------------------
session_ts, ep_feats = compute_dlc_behavior_features(
    dlc_mat,
    dlc_columns,
    targets_periods_times,
    smooth_size=20,
    p_thr=0.8,
    body_points=("neck", "lower_spine", "tail_base"),
    use_ears_for_hd=True,
)
    
save_dlc_features_to_h5(
    h5_path=snakemake.output[0],
    session_ts=session_ts,
    episode_feats=ep_feats,
    group="dlc_features",
    mode="a",
    overwrite_group=True,
    extra_meta={"smooth_size": 20, "p_thr": 0.8},
)