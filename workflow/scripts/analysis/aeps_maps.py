import os, sys
import h5py
import numpy as np
import matplotlib.pyplot as plt


# TODO: move to params?
to_plot = ['P1', 'N1', 'P3']
sigma = 0.5
bin_size = 0.04
x_min, x_max = -0.5, 0.5
y_min, y_max = -0.5, 0.5

x_range = x_max - x_min
y_range = y_max - y_min
y_bin_count = int(np.ceil(y_range / bin_size))
x_bin_count = int(np.ceil(x_range / bin_size))

# load timeline / events
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])

# load AEP metrics
AEP_metrics = {}
with h5py.File(snakemake.input[1], 'r') as f:
    for area in f:
        metrics = {}
        for m in f[area]:
            metrics[m] = np.array(f[area][m])
        AEP_metrics[area] = metrics

# filter, e.g. BGR + TGT only
idxs_filt = np.where( (sound_events[:, 1] == 1) | (sound_events[:, 1] == 2) )[0]

# positions where AEPs happened
x_pos = tl[sound_events[idxs_filt][:, 2].astype(np.int32)][:, 1]
y_pos = tl[sound_events[idxs_filt][:, 2].astype(np.int32)][:, 2]

def get_aep_map(metric):
    aep_map = np.zeros([x_bin_count, y_bin_count])
    for x_bin_idx in range(x_bin_count - 1):
        for y_bin_idx in range(y_bin_count - 1):
            x_l = x_min + bin_size * x_bin_idx
            x_r = x_min + bin_size * (x_bin_idx + 1)
            y_l = y_min + bin_size * y_bin_idx
            y_r = y_min + bin_size * (y_bin_idx + 1)

            # select all AEP responses in that bin
            x_idxs = np.where( (x_pos > x_l) & (x_pos < x_r) )[0]
            y_idxs = np.where( (y_pos > y_l) & (y_pos < y_r) )[0]
            sel_idxs = np.intersect1d(x_idxs, y_idxs)
            if len(sel_idxs) > 0:
                aep_map[x_bin_idx, y_bin_idx] = metric[sel_idxs].mean()
    return aep_map

# plotting
cols = len(to_plot)
rows = 2 * len(AEP_metrics)
fig = plt.figure(figsize=(10, rows*3.5))

for i, (area, metrics) in enumerate(AEP_metrics.items()):
    for j, m_type in enumerate(['norm', 'raw']):
        for k, m_name in enumerate(to_plot):
            m_vals = metrics['%s_%s' % (m_name, m_type)]
            aep_map = get_aep_map(m_vals[idxs_filt])
            
            ax = fig.add_subplot(rows, cols, 6*i + 3*j + k + 1)
            ax.imshow(aep_map, cmap='jet')
            ax.set_title("%s: %s (%s)" % (area, m_name, m_type), fontsize=14)

fig.savefig(snakemake.output[0])