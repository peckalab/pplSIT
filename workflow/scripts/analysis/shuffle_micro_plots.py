import os
import sys
import h5py
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import get_unit_names_sorted


# -------------------------
# read configs / unit names
# -------------------------
with h5py.File(snakemake.input[0], 'r') as f:
    cfg = json.loads(f['processed'].attrs['parameters'])

with h5py.File(snakemake.input[1], 'r') as f:
    units_to_plot = get_unit_names_sorted([unit for unit in f['BGR']])  # should always exist in profiles

# configuration
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}
colors = {0: 'indigo', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}

bgr_dur = cfg['sound']['sounds']['background']['duration']  # in seconds
latency = cfg['sound']['latency']

cols = len(event_types)
units_per_page = 10  # adjust as you like

# -------------------------
# preload data once
# -------------------------
profile_cache = {}
with h5py.File(snakemake.input[1], 'r') as f:
    for ev_type in event_types:
        ev_name = ev_names[ev_type]
        profile_cache[ev_name] = {}
        for unit_name in units_to_plot:
            profile_cache[ev_name][unit_name] = np.array(f[ev_name][unit_name]['profile_stats'])

shuffled_cache = {}
with h5py.File(snakemake.input[2], 'r') as f:
    for ev_type in event_types:
        ev_name = ev_names[ev_type]
        shuffled_cache[ev_name] = {}
        for unit_name in units_to_plot:
            shuffled_cache[ev_name][unit_name] = np.array(f[ev_name][unit_name]['shuffled'])

# -------------------------
# multipage PDF
# -------------------------
with PdfPages(snakemake.output[0]) as pdf:
    for page_start in range(0, len(units_to_plot), units_per_page):
        units_chunk = units_to_plot[page_start:page_start + units_per_page]
        rows = len(units_chunk)

        fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 2.8 * rows), squeeze=False)

        for i, unit_name in enumerate(units_chunk):
            for j, ev_type in enumerate(event_types):
                ev_name = ev_names[ev_type]

                profile_stats = profile_cache[ev_name][unit_name]
                shuffled = shuffled_cache[ev_name][unit_name]

                ax = axes[i][j]

                for k, c_stats in enumerate([shuffled, profile_stats]):
                    bin_size = c_stats[0][1] - c_stats[0][0]
                    clr = colors[ev_type] if k == 1 else 'black'
                    label = ev_name if k == 1 else 'SHUF'

                    ax.plot(c_stats[0] + bin_size, c_stats[1], color=clr, lw=2, label=label)
                    ax.fill_between(c_stats[0] + bin_size, c_stats[3], c_stats[4], color=clr, alpha=0.4)

                ax.set_xlim(-latency, latency)
                ax.set_ylim(bottom=0)
                ax.axvline(0, color='black', ls='--')
                ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
                ax.axvspan(-latency, -latency + bgr_dur, alpha=0.3, color='gray')
                ax.set_title(f"{unit_name}", fontsize=12)
                ax.grid()

                if i == 0:
                    ax.legend(fontsize=8)

            axes[i][0].set_ylabel("Firing Rate, Hz", fontsize=11)

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)