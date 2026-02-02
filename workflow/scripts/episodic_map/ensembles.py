import os, sys
import h5py
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy import stats
from sklearn import decomposition

import yaml
from datetime import datetime
from pathlib import Path

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)


#cfg = snakemake.config['trajectories']

# -----------------------------
# Reading datasets
# -----------------------------
meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
actm_file = snakemake.input[2]
unit_file = snakemake.input[3]
sphl_file = snakemake.input[4]
cglm_file = snakemake.input[5]

# ------------- loading data -------------

with h5py.File(unit_file, 'r') as f:
    unit_ids = [x for x in f]

with h5py.File(cglm_file, 'r') as f:
    coef_full = np.array(f['fit/coef_full'])  # unit ids should match
    names = f['design/X_names'][...].astype(str)

# segments - TGt
with h5py.File(segm_file, 'r') as f:
    idxs_tgt_ev = np.array(f['idxs_tgt_sta_succ'])
    tgt_succ_start_mx = np.array(f['tgt_sta_succ_mx'])
idxs_tgt_for_episodes = []
counter = 0
for rec in tgt_succ_start_mx:
    diff = rec[1] - rec[0] + 1
    idxs_tgt_for_episodes.append([counter, counter + diff])
    counter += diff
idxs_tgt_for_episodes = np.array(idxs_tgt_for_episodes) * 5  # convert to 50ms-bin scale

# need to convert idices from pulse 250ms-bin scale to 50ms-bin scale
idxs_tgt_ev_50ms = []
offset = 0
for idx in idxs_tgt_ev:
    idxs_tgt_ev_50ms += list(range(idx*5 + offset, idx*5+5 + offset))
idxs_tgt_ev_50ms = np.array(idxs_tgt_ev_50ms)

# load condition-specific units
mx_type = 'mx'
with h5py.File(actm_file, 'r') as f:
    X_counts_50ms = np.array(f['mx_50ms'][mx_type])  # units x timebins
    t_edges_50ms  = np.array(f['mx_50ms']['bins'])
mx_z_cond  = stats.zscore(X_counts_50ms.T[idxs_tgt_ev_50ms].T, axis=1)
mx_z_clean = np.nan_to_num(mx_z_cond, nan=0.0)
unit_mx_state = mx_z_clean.T
    
# PCA and loadings
ev_pca = decomposition.PCA(n_components=5)
ev_X   = ev_pca.fit_transform(unit_mx_state)
loadings = ev_pca.components_[0]

# MRLs
MRLs_tgt, MRLs_bgr, p_vals = [], [], []
with h5py.File(sphl_file, 'r') as f:
    for unit_id in unit_ids:
        MRLs_tgt.append( np.array(f[f'tgt/{unit_id}/MRL_real']) )
        MRLs_bgr.append( np.array(f[f'bgr_sta/{unit_id}/MRL_real']) )
        p_vals.append( np.array(f[f'tgt/{unit_id}/p_value']) )
MRLs_tgt = np.array(MRLs_tgt)
MRLs_bgr = np.array(MRLs_bgr)
p_vals = np.array(p_vals)

# -----------------------------

def get_unit_ids_by_idxs(idxs):
    return [unit_ids[i] for i in idxs]

# 1. ensemble with high PC1 target loadings
idxs_low_PC  = np.where(loadings < -0.05)[0]  # check the threshold - maybe closer to 0?

# 2. ensemble with significant target MRLs
idxs_signif_tgt = np.where(p_vals < 0.05)[0]
idxs_ensemble1 = np.union1d(idxs_low_PC, idxs_signif_tgt)

# 3. All above + movement-related units (example)
idxs_mov = np.where(coef_full[:, 8] > 0.005)[0]  # body RMS coef
idxs_ensemble2 = np.union1d(idxs_ensemble1, idxs_mov)


ensembles = [
    {
        "name": "gates",
        "description": "All 'gating' units",
        "params": {
            "pc1_tgt_loading_less_than": -0.05,
        },
        "unit_ids": get_unit_ids_by_idxs(idxs_low_PC),
    },
    {
        "name": "sound_phase_lock",
        "description": "All significantly sound phase-locked units",
        "params": {
            "phase_lock_p": 0.001,
        },
        "unit_ids": get_unit_ids_by_idxs(idxs_signif_tgt),
    },
    {
        "name": "move",
        "description": "All significantly sound phase-locked units",
        "params": {
            "body_RMS_more_than": 0.005,
        },
        "unit_ids": get_unit_ids_by_idxs(idxs_mov),
    },
    {
        "name": "gate_lock",
        "description": "Negative PC1 loading + phase locking",
        "params": {
            "pc1_tgt_loading_less_than": -0.05,
            "phase_lock_p": 0.001,
        },
        "unit_ids": get_unit_ids_by_idxs(idxs_ensemble1),
    },
    {
        "name": "gate_lock_move",
        "description": "Negative PC1 loading + phase locking + movement related",
        "params": {
            "pc1_tgt_loading_less_than": -0.05,
            "phase_lock_p": 0.001,
            "body_RMS_more_than": 0.005,
        },
        "unit_ids": get_unit_ids_by_idxs(idxs_ensemble2),
    },
]

def save_ensembles_yaml(ensembles, out_path):
    """
    ensembles: list of dicts, each with keys:
        name, unit_ids, description (optional), params (optional), source (optional)
    out_path: str or Path
    """
    for ens in ensembles:
        ens.setdefault("description", "")
        ens.setdefault("params", {})
        ens.setdefault("created_at", datetime.now().isoformat(timespec="seconds"))

    data = {"ensembles": ensembles}

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as f:
        yaml.safe_dump(
            data,
            f,
            sort_keys=False,
            default_flow_style=False
        )

save_ensembles_yaml(ensembles, snakemake.output[0])