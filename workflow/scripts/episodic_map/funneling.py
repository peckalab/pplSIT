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

compute_and_save_funneling_from_core_h5(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],            
    group="funneling",
    early_window_s=(0.0, 1.5),
    late_window_s=(4.5, 6.0),
    n_time_shuffles=50,
    n_random_windows_pool=500,
    n_random_bootstrap=50,
    seed=0,
)
