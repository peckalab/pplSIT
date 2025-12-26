import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.regression import *


#cfg = snakemake.config['trajectories']

epsi_file = snakemake.input[0]

for i, s_type in enumerate(["S_latent", "S_pv"]):
    spec = RegressionSpec(
        y_key="S_latent",
        x_keys=("dt_s", "dspace_m", "dcenter_m", "dhd_rad", "dhd_rel_center_rad", "dturn", "dstill"),
        model="ridge",
        ridge_alpha=1.0,
        n_perm_episode_shuffle=200,
        n_perm_y_shuffle=200,
        n_perm_circshift=0,
        random_seed=0,
    )

    run_session_regression(
        h5_path=snakemake.input[0],
        out_h5=snakemake.output[i],
        spec=spec
    )
