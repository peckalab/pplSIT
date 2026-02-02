import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.drift_PCA import *

# uncomment if needed
# meta_file = snakemake.input[0]
# with h5py.File(meta_file, 'r') as f:
#     tl = np.array(f['processed']['timeline'])
#     tgt_mx = np.array(f['processed']['target_matrix'])
#     sound_events = np.array(f['processed']['sound_events'])

actm_file = snakemake.input[0]
with h5py.File(actm_file, 'r') as f:
    spike_counts = np.array(f['mx_50ms']['mx'])
    #t_edges_50ms= np.array(f['mx_50ms']['bins'])
    
episode_sets = {
    "target":   "episodes/target_bin_windows",
    "all_sta":  "episodes_sta/all/raw/bin_windows",
    "bgr_sta":  "episodes_sta/bgr/raw/bin_windows",
    "sil_sta":  "episodes_sta/sil/raw/bin_windows",
}

compute_and_save_drift_pcs_single_session(
    traj_h5_path=snakemake.input[1], 
    out_h5_path=snakemake.output[0],
    spike_counts=spike_counts, 
    episode_sets=episode_sets
)
