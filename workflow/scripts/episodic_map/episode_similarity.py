import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.episode_similarity import *


#cfg = snakemake.config['trajectories']

# STATE_KEYS = [
#     "idxs_tgt_sta_succ",
#     "idxs_bgr_sta",
#     "idxs_sil_sta",
# ]

# -----------------------------
# Reading datasets
# -----------------------------
actm_file = snakemake.input[0]
traj_file = snakemake.input[1]
behf_file = snakemake.input[2]


#tgt_succ_mx = tgt_mx[tgt_mx[:, 4] == 1]
#targets_periods_bins = np.column_stack([tgt_succ_mx[:, 0] + 2, tgt_succ_mx[:, 0] + 24])

# 250 ms binning to match PCA representation
with h5py.File(actm_file, 'r') as f:
    X_counts = np.array(f['mx_50ms']['mx'])
    t_edges  = np.array(f['mx_50ms']['bins'])
    

res = compute_episode_similarity_and_drivers_to_h5(
    out_h5_path=snakemake.output[0],
    trajectories_h5=traj_file,
    behavior_h5=behf_file,

    # optional PV counts (for option B)
    X_counts_50ms=X_counts,
    X_is_neurons_by_time=True,

    # episode selection (indices into trajectories episodes)
    #ep_select_idx: Optional[np.ndarray] = None,

    # similarity choices
    latent_similarity="cosine",  # "cosine" or "neg_euclid"
    pv_similarity="cosine",      # "cosine" (recommended)

    # phase window
    #strip_first_bins: int = 10,         # e.g. 0.5 s at 50 ms
    #use_time_bins: Optional[np.ndarray] = None,  # overrides computed time bins
)