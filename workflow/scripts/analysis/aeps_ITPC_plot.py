import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_tf_matrix(fig, ax, matrix, title, times, freqs, vmin=None, vmax=None):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    
    im = ax.imshow(matrix, aspect='auto', origin='lower',
               extent=[times[0], times[-1], freqs[0], freqs[-1]],
               cmap='viridis', vmin=vmin, vmax=vmax)
    fig.colorbar(im, cax=cax, orientation='vertical')
    ax.set_xlabel('Time (ms)')
    #ax.set_ylabel('Frequency (Hz)')
    ax.set_title(title)


state_labels = ['Target', 'Background (sta)', 'Background (run)', 'No stimulus (sta)', 'No stimulus (run)']
state_colors = ['tab:orange', 'tab:blue', 'red', 'grey', 'black']
state_names  = ['idxs_tgt_sta_succ', 'idxs_bgr_sta', 'idxs_bgr_run', 'idxs_sil_sta', 'idxs_sil_run']


fig, axes = plt.subplots(4, len(state_names), figsize=(16, 16))

for idx, state_name in enumerate(state_names):
    with h5py.File(snakemake.input[0], 'r') as f:
        times  = np.array(f[state_name]['times'])
        freqs  = np.array(f[state_name]['frequencies'])
        power  = np.array(f[state_name]['power_real'])
        itpc   = np.array(f[state_name]['itpc_real'])
        itpc_z = np.array(f[state_name]['itpc_shuf_z'])
        mask   = np.array(f[state_name]['non_phaselocked_mask'])  # ['significant_mask']

    plot_tf_matrix(fig, axes[0][idx], power, f'Power - {state_labels[idx]}', times, freqs)
    plot_tf_matrix(fig, axes[1][idx], itpc, f'ITPC - {state_labels[idx]}', times, freqs, vmin=0, vmax=1)
    plot_tf_matrix(fig, axes[2][idx], itpc_z, f'ITPC_z - {state_labels[idx]}', times, freqs)
    plot_tf_matrix(fig, axes[3][idx], mask, f'ITPC (Mask) - {state_labels[idx]}', times, freqs, vmin=0, vmax=1)

fig.tight_layout()
fig.savefig(snakemake.output[0])