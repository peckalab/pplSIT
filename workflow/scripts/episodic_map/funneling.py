import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.funneling import *

#cfg = snakemake.config['trajectories']

# 6 seconds long targets
compute_and_save_funneling_from_core_h5(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],            
    group="funneling",
    early_window_s=(0.0, 1.0),
    late_window_s=(5.0, 6.0),
    representation="raw",              
    episode_kind="target",       
    use_resid=False,
)

compute_and_save_funneling_from_core_h5(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],            
    group="funneling",
    early_window_s=(0.0, 1.0),
    late_window_s=(5.0, 6.0),
    use_resid=True,
)

# 3 seconds long to compare with STAs
compute_and_save_funneling_from_core_h5(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],            
    group="funneling_early",
    early_window_s=(0.0, 1.0),
    late_window_s=(2.0, 3.0),
    use_resid=False,
)

compute_and_save_funneling_from_core_h5(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],            
    group="funneling_early",
    early_window_s=(0.0, 1.0),
    late_window_s=(2.0, 3.0),
    use_resid=True,
)

compute_and_save_funneling_from_core_h5(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],            
    group="funneling_late",
    early_window_s=(3.0, 4.0),
    late_window_s=(5.0, 6.0),
    use_resid=False,
)

compute_and_save_funneling_from_core_h5(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],            
    group="funneling_late",
    early_window_s=(3.0, 4.0),
    late_window_s=(5.0, 6.0),
    use_resid=True,
)

# non-target stationary funneling
for sta_kind in ["bgr", "sil", "all"]:
    for rep in ["raw"]:
        compute_and_save_funneling_from_core_h5(
            snakemake.input[0], snakemake.output[0],
            early_window_s=(0.0, 1.0),
            late_window_s=(2.0, 3.0),
            episode_kind="sta",
            sta_kind=sta_kind,
            representation=rep,
            slice_mode="whole",
            slice_len_s=None,   # sta episodes already 3s
        )
