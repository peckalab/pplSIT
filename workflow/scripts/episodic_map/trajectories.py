import os, sys
import h5py
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.trajectories import *


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

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

with h5py.File(actm_file, 'r') as f:
    X_counts_50ms = np.array(f['mx_50ms']['mx'])
    t_edges_50ms  = np.array(f['mx_50ms']['bins'])
    
with h5py.File(segm_file, 'r') as f:
    state_periods_pulse = {}
    for idxs_name in STATE_KEYS:
        state_periods_pulse[idxs_name] = np.array(f[idxs_name])[:, :2].astype(np.int32)

out = build_core_trajectories_h5(
    out_h5_path=snakemake.output[0],
    X_counts_50ms=X_counts_50ms,
    state_periods_pulse=state_periods_pulse,
    t_edges_50ms=t_edges_50ms,
    X_is_neurons_by_time=True,      # set correctly for your matrix
    episode_win_s=cfg['episode_win_s'],
    pca_n_components=cfg['pca_n_components'],
    transform=cfg['transform'],
    behavior_h5_path=behf_file,
)

# def build_core_trajectories_h5(
#     out_h5_path: str,
#     X_counts_50ms: np.ndarray,
#     state_periods_pulse: dict,
#     t_edges_50ms: np.ndarray | None = None,
#     *,
#     X_is_neurons_by_time: bool = True,
#     bins_per_pulse: int = 5,           # 250ms / 50ms = 5
#     pulse_end_inclusive: bool = True,

#     # which states to use
#     key_target: str  = "tgt_sta_succ_mx",
#     key_bgr_sta: str = "bgr_sta_mx",
#     key_sil_sta: str = "sil_sta_mx",

#     # PCA fit choices
#     transform: str = "sqrt",
#     pca_n_components: int = 20,

#     # episode windowing
#     episode_win_s: float = 6.0,       # 6 s target episodes
#     bin_size_s: float = 0.05,         # 50 ms
#     align: str = "start",

#     # residualization inputs
#     behavior_h5_path: Optional[str] = None,
#     do_context_residual: bool = True,
#     do_stim_phase_residual: bool = True,
#     context_cov_spec: ContextCovariateSpec = ContextCovariateSpec(),
#     ridge_alpha_resid: float = 1.0,
# )



