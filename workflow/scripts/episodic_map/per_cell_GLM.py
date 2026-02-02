import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.per_cell_GLM import *


meta_file = snakemake.input[0]
traj_file = snakemake.input[1]
behf_file = snakemake.input[2]
dpca_file = snakemake.input[3]
actm_file = snakemake.input[4]

meta_file = snakemake.input[0]
with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])
stim_times_s = sound_events[:,0]

with h5py.File(actm_file, 'r') as f:
    Y_counts    = np.array(f['mx_50ms']['mx']).T  # (T, N)
    bin_edges_s = np.array(f['mx_50ms']['bins'])

spec = GLMSpec(
    model="poisson",
    alpha=1.0,
    n_splits=5,
    drift_pcs_k=3,
    stim_f_hz=4.0,
    stim_impulse_kernel_s=0.5,
    include_phi2=True,
)

compute_and_save_cell_glm_single_session(
    out_h5_path=snakemake.output[0],
    Y_counts=Y_counts,                       # (T, N)
    bin_edges_s=bin_edges_s,          # (T+1,)
    bin_size_s=0.05,
    behavior_h5_path=behf_file,
    stim_times_s=stim_times_s,        # (n_stim,)
    trajectories_h5_path=traj_file,
    drift_pcs_h5_path=dpca_file,
    drift_condition="target",
    fit_period="all",
    spec=spec,
    unit_meta=None,
)

compute_and_save_cell_glm_single_session(
    out_h5_path=snakemake.output[1],
    Y_counts=Y_counts,                       # (T, N)
    bin_edges_s=bin_edges_s,          # (T+1,)
    bin_size_s=0.05,
    behavior_h5_path=behf_file,
    stim_times_s=stim_times_s,        # (n_stim,)
    trajectories_h5_path=traj_file,
    drift_pcs_h5_path=dpca_file,
    drift_condition="target",
    fit_period="all_sta",
    spec=spec,
    unit_meta=None,
)