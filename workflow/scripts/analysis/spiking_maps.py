import matplotlib.pyplot as plt
import os, sys
import h5py
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.hdf import H5NAMES


with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    
# load units
unit_names, single_units, spike_times, spike_idxs = [], {}, {}, {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = [x for x in f]
with h5py.File(snakemake.input[1], 'r') as f:
    for unit_name in unit_names:
        spike_idxs[unit_name] = np.array(f[unit_name][H5NAMES.spike_idxs['name']])

limits = (-0.5, 0.5, -0.5, 0.5)
cols = 5
rows = int(np.ceil(len(unit_names)/cols))
fig = plt.figure(figsize=(14, rows*3.2))

for i, unit_name in enumerate(unit_names):
    run_idxs = np.where(tl[:, 3] > 0.04)[0]
    
    # animal positions
    a_pos = tl[run_idxs][:, 1:3]

    # spike positions
    s_pos = tl[np.intersect1d(run_idxs, spike_idxs[unit_name])][:, 1:3]

    ax = fig.add_subplot(rows, cols, i+1)
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.add_patch(plt.Circle((0, 0), 0.46, color='r', fill=False))
    ax.scatter(a_pos[:, 0], a_pos[:, 1], alpha=0.1, s=2)
    ax.scatter(s_pos[:, 0], s_pos[:, 1], alpha=0.4, s=10, c='black')
    ax.set_aspect('equal')
    ax.set_title(unit_name, fontsize=14)

fig.savefig(snakemake.output[0])