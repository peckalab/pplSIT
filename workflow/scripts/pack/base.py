import os, json
import h5py
import numpy as np
from scipy import signal
from scipy.ndimage import median_filter

import sys
# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.pack import pack_base

# ---- Snakemake entry point ----
pack_base(
    snakemake.input.positions,
    snakemake.input.events,
    snakemake.input.sounds,
    snakemake.input.islands,
    snakemake.input.cfg,
    snakemake.input.manual,
    snakemake.output.base,
    drift_coeff=snakemake.config["pack"]["drift_coeff"],
)
