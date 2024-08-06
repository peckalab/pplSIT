import h5py, os, sys
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.segmentation import manifold_segmentation

cfg     = snakemake.config['nMAP_segmentation']
s_path  = os.path.dirname(snakemake.input[0])
source  = os.path.dirname(os.path.dirname(s_path))
session = os.path.basename(s_path)
animal  = session.split('_')[0]

# reading datasets
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    events = np.array(f['processed']['sound_events'])
with h5py.File(snakemake.input[1], 'r') as f:
    speed = np.array(f['speed'])
    hd    = np.array(f['hd'])
with h5py.File(snakemake.input[2], 'r') as f:
    fit = np.array(f[cfg['ft']][str(cfg['fp'])])

d_map, fit_labels, fit_labels_sel, labels_ev, idxs_tgt_succ_state_ev, tgt_stats, tgt_stats_shuf = manifold_segmentation(events, tgt_mx, fit, cfg)

# saving everything
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('density_map', data=d_map)
    f.create_dataset('segmentation', data=fit_labels)
    f.create_dataset('segmentation_TGT_succ', data=fit_labels_sel)
    f.create_dataset('labels_ev', data=labels_ev)
    f.create_dataset('idxs_tgt_succ_state_ev', data=idxs_tgt_succ_state_ev)
    f.create_dataset('tgt_stats_shuf', data=tgt_stats_shuf)
    f.create_dataset('tgt_stats', data=tgt_stats)
