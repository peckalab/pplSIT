import os, sys
import h5py
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.trajectories_tgt import *
from utils.episodic_map.trajectories_sta import *

import yaml
from pathlib import Path

cfg = snakemake.config['trajectories']

STATE_KEYS = [
    "tgt_sta_succ_mx",
    "bgr_sta_mx",
    "sil_sta_mx",
]

# -----------------------------
# Reading datasets
# -----------------------------
meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
actm_file = snakemake.input[2]
behf_file = snakemake.input[3]
ensm_file = snakemake.input[4]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

with h5py.File(actm_file, 'r') as f:
    X_counts_50ms = np.array(f['mx_50ms']['mx'])  # units x time
    t_edges_50ms  = np.array(f['mx_50ms']['bins'])
    unit_ids_all = eval(f['mx_50ms'].attrs['unit_ids'])
    
with h5py.File(segm_file, 'r') as f:
    state_periods_pulse = {}
    for idxs_name in STATE_KEYS:
        state_periods_pulse[idxs_name] = np.array(f[idxs_name])[:, :2].astype(np.int32)


ensm_file = Path(ensm_file)
with open(ensm_file, "r") as f:
    data = yaml.safe_load(f)

mode = "w"
for ens in data.get("ensembles", []):
    name = ens.get("name")
    unit_ids_to_filt = ens.get("unit_ids")

    idxs_filt_pos = [unit_ids_all.index(uid) for uid in unit_ids_to_filt if uid in unit_ids_all]
    idxs_filt_neg = [x for x in range(len(unit_ids_all)) if x not in idxs_filt_pos]

    out1 = build_core_trajectories_h5(
        out_h5_path=snakemake.output[0],
        subgroup=name + "_only",
        mode=mode,
        X_counts_50ms=X_counts_50ms[idxs_filt_pos],
        state_periods_pulse=state_periods_pulse,
        t_edges_50ms=t_edges_50ms,
        X_is_neurons_by_time=True,      # set correctly for your matrix
        episode_win_s=cfg['episode_win_s'],
        pca_n_components=cfg['pca_n_components'],
        transform=cfg['transform'],
        behavior_h5_path=behf_file,
    )

    if mode == "w":
        mode = "r+"

    out1 = build_core_trajectories_h5(
        out_h5_path=snakemake.output[0],
        subgroup=name + "_excluded",
        mode=mode,
        X_counts_50ms=X_counts_50ms[idxs_filt_neg],
        state_periods_pulse=state_periods_pulse,
        t_edges_50ms=t_edges_50ms,
        X_is_neurons_by_time=True,      # set correctly for your matrix
        episode_win_s=cfg['episode_win_s'],
        pca_n_components=cfg['pca_n_components'],
        transform=cfg['transform'],
        behavior_h5_path=behf_file,
    )

    
