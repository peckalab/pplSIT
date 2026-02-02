import os, sys
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.decoder_ID_mean import *


run_episode_mean_only_decoding_default(
    core_h5_path=snakemake.input[0],
    out_h5_path=snakemake.output[0],
    D_use=10,
    strip_first_s=0.5,
    early_s=2.0,
    late_s=2.0,
)

specs = [
    MeanOnlyDecodeSpec(
        representation='raw', mean_window='whole',
        strip_first_s=0.5, early_s=2.0, late_s=2.0,
        zscore_across_episodes=True
    ),
    MeanOnlyDecodeSpec(
        representation='raw', mean_window='early',
        strip_first_s=0.5, early_s=2.0, late_s=2.0,
        zscore_across_episodes=True
    ),
    MeanOnlyDecodeSpec(
        representation='raw', mean_window='late',
        strip_first_s=0.5, early_s=2.0, late_s=2.0,
        zscore_across_episodes=True
    ),
]
run_episode_mean_only_decoding_single_session(
    core_h5_path=snakemake.input[1],
    out_h5_path=snakemake.output[1],
    subgroup='gate_lock_move_excluded',
    D_use=10, 
    specs=specs
)