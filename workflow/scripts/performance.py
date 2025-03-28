import os, sys
import h5py
import json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.performance import calculate_performance


with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    tgt_matrix = np.array(f['processed']['target_matrix'])
    trials = np.array(f['processed']['trial_idxs'])
    cfg = json.loads(f['processed'].attrs['parameters'])

tr_succ = trials[trials[:, 5] == 1]
proportion_correct = len(tr_succ) / (len(trials))

performance = calculate_performance(tl, trials, cfg)

# write out
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('proportion_correct', data=proportion_correct)
    perf = f.create_dataset('performance', data=performance)
    perf.attrs['columns'] = 'performance_median, performance_upper_CI, performance_lower_CI, \
        chance_median, chance_upper_CI, chance_lower_CI, time'