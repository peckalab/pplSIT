import h5py, os, sys, json
import numpy as np
from scipy import stats

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)


# read datasets
meta_file = snakemake.input[0]
unit_file = snakemake.input[1]
smk_cfg   = snakemake.config['ensembles']

s_path  = os.path.dirname(meta_file)
session = os.path.basename(s_path)
animal  = session.split('_')[0]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

inst_rates = {}
with h5py.File(unit_file, 'r') as f:
    for unit_name in f:
        inst_rates[unit_name] = np.array(f[unit_name]['inst_rate'])


# building ensebles - timeline resolution
ensembles = {}
for e_name, units in smk_cfg[session].items():
    ensemble_mx = np.zeros([len(units), len(tl)])
    for i, unit_id in enumerate(units):
        #ensemble_mx[i] = inst_rate(spike_times[unit_id], tl[:, 0], k_width=50)
        ensemble_mx[i] = inst_rates[unit_id]
        ensemble_mx[i] = stats.zscore(ensemble_mx[i])

    ensembles[e_name] = ensemble_mx.mean(axis=0)

# dump to H5
with h5py.File(snakemake.output[0], 'w') as out_file:
    for e_name, e_inst_rate in ensembles.items():
        out_file.create_dataset(e_name, data=e_inst_rate)