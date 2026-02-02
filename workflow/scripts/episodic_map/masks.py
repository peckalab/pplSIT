import os, sys
import h5py
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.masks import *

import yaml
from pathlib import Path

cfg = snakemake.config['trajectories']
ens_name = 'gates'

#meta_file = snakemake.input[0]
actm_file = snakemake.input[0]
ensm_file = snakemake.input[1]

# with h5py.File(meta_file, 'r') as f:
#     tl = np.array(f['processed']['timeline'])
#     tgt_mx = np.array(f['processed']['target_matrix'])
#     sound_events = np.array(f['processed']['sound_events'])

with h5py.File(actm_file, 'r') as f:
    X_counts_50ms = np.array(f['mx_50ms']['mx'])  # units x time
    t_edges_50ms  = np.array(f['mx_50ms']['bins'])
    unit_ids_all = eval(f['mx_50ms'].attrs['unit_ids'])
    
ensm_file = Path(ensm_file)
with open(ensm_file, "r") as f:
    ensembles = yaml.safe_load(f)['ensembles']

gate_ens = [e for e in ensembles if e['name'] == 'gates'][0]
units_mask = gate_ens['unit_ids']
idxs_filt_pos = np.array([unit_ids_all.index(uid) for uid in units_mask], dtype=np.int32)

# gate_counts: shape (n_gate_units, T)
gate_counts = X_counts_50ms[idxs_filt_pos]
out = build_gate_engagement_mask(gate_counts, smooth_sigma_bins=2, thr_quantile=0.8, min_on_bins=2)

save_gate_engagement_to_h5(
    h5_path=snakemake.output[0],
    h5_group=ens_name,
    out=out,
)