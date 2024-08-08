import os, sys
import h5py, json
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import get_unit_names_sorted


# reading some configs
with h5py.File(snakemake.input[0], 'r') as f:
    cfg = json.loads(f['processed'].attrs['parameters'])
with h5py.File(snakemake.input[1], 'r') as f:
    units_to_plot = get_unit_names_sorted([unit for unit in f['BGR']])  # should always exist in profiles

# configuration
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}
colors = {0: 'indigo', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
bgr_dur    = cfg['sound']['sounds']['background']['duration']  # in seconds
cols = len(event_types)
rows = len(units_to_plot)


# plot figures
fig, axes = plt.subplots(rows, cols, figsize=(3*cols, 3*rows))

for i, unit_name in enumerate(units_to_plot):
    for j, ev_type in enumerate(event_types):
    
        # read PSTH profile
        with h5py.File(snakemake.input[1], 'r') as f:
            # bins, shuffled mean, std, perc 5, perc 95
            profile_stats = np.array(f[ev_names[ev_type]][unit_name]['profile_stats'])
        
        # read shuffled stats
        with h5py.File(snakemake.input[2], 'r') as f:
            # bins, shuffled mean, std, perc 5, perc 95
            shuffled = np.array(f[ev_names[ev_type]][unit_name]['shuffled'])
        
        ax = axes[i][j]
        for k, c_stats in enumerate([shuffled, profile_stats]):
            clr = colors[ev_type] if k == 1 else 'black'
            label = ev_names[ev_type] if k == 1 else 'SHUF'
            ax.plot(c_stats[0], c_stats[1], color=clr, lw=2, label=label)
            ax.fill_between(c_stats[0], c_stats[3], c_stats[4], color=clr, alpha=0.4)
        
        ax.set_ylim(bottom=0)
        ax.axvline(0, color='black', ls='--')
        ax.set_title("%s" % unit_name, fontsize=14)
        ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
        ax.legend()
        ax.grid()
    
    axes[i][0].set_ylabel("Firing Rate, Hz", fontsize=14)
            
    fig.tight_layout()
    fig.savefig(snakemake.output[0])

