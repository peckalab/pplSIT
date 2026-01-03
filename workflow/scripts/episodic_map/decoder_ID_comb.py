import os, sys
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.decoder_ID_comb import *

run_window_generalization_for_session(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],
    L_s_list=(0.25, 0.5, 1.0, 2.0, 3.0),
    representations=("raw", "resid", "resid_stim", "resid_ctx_stim"),
    mean_subtract_modes=("none", "episode_global", "train_window")
)
