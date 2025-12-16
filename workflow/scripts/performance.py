import os, sys
import h5py
import json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.performance import calculate_performance, plot_session_metrics, plot_performance


with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    tgt_matrix = np.array(f['processed']['target_matrix'])
    trials = np.array(f['processed']['trial_idxs'])
    cfg = json.loads(f['processed'].attrs['parameters'])
    islands = np.array(f['raw'].get('islands', None))

distractor_fail = cfg['experiment'].get('distractor_fail', False)

tr_succ = trials[trials[:, 5] == 1]
proportion_correct = len(tr_succ) / (len(trials))

if distractor_fail:
    # distractor fail is a trial the was incorrect (trial[:,5]==0) but trial time was shorter than trial_duration
    trial_time = tl[trials[:, 1].astype(np.int32)][:, 0] - tl[trials[:, 0].astype(np.int32)][:, 0]
    tr_distractor_fail = trials[(trials[:, 5] == 0) & 
                                (trial_time < cfg['experiment']['trial_duration'])]
    proportion_distractor_fail = len(tr_distractor_fail) / len(trials)

performance = calculate_performance(tl, trials, cfg, islands)


# write out
with h5py.File(snakemake.output.performance, 'w') as f:
    f.create_dataset('proportion_correct', data=proportion_correct)
    perf = f.create_dataset('performance', data=performance)
    if distractor_fail and performance.shape[1] > 7:
        f.create_dataset('proportion_distractor_fail', data=proportion_distractor_fail)
        perf.attrs['columns'] = 'time, chance_median, chance_upper_CI, chance_lower_CI, performance_median, performance_upper_CI, performance_lower_CI, \
            distractor_fail_median, distractor_fail_upper_CI, distractor_fail_lower_CI'
    else:
        perf.attrs['columns'] = 'time, chance_median, chance_upper_CI, chance_lower_CI, performance_median, performance_upper_CI, performance_lower_CI'

# plot performance figure
plot_performance(cfg, performance, snakemake.output.performance_figure)

# plot session metrics figure
plot_session_metrics(tl, trials, cfg, snakemake.output.session_metrics_figure)