import matplotlib.pyplot as plt
import os, sys
import h5py
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.hdf import H5NAMES
from utils.spatial import pol2cart


# load units
unit_names, spike_idxs, f_maps, f_PFRs, spat_info, peak_FRs = [], {}, {}, {}, {}, {}
with h5py.File(snakemake.input[0], 'r') as f:
    unit_names = [x for x in f]
with h5py.File(snakemake.input[0], 'r') as f:
    for unit_name in unit_names:
        spike_idxs[unit_name] = np.array(f[unit_name][H5NAMES.spike_idxs['name']])
        f_maps[unit_name]     = np.array(f[unit_name][H5NAMES.f_maps['name']])
        f_PFRs[unit_name]     = np.array(f[unit_name][H5NAMES.pfr_center['name']])
        spat_info[unit_name]  = np.array(f[unit_name][H5NAMES.spat_info['name']])
        peak_FRs[unit_name]   = np.array(f[unit_name][H5NAMES.peak_FR['name']])

limits = (-0.5, 0.5, -0.5, 0.5)
cols = 5
rows = int(np.ceil(len(unit_names)/cols))
fig = plt.figure(figsize=(14, rows*3.2))

for i, unit_name in enumerate(unit_names):
    x_PFR, y_PFR = pol2cart(f_PFRs[unit_name][0], f_PFRs[unit_name][1])

    ax = fig.add_subplot(rows, cols, i+1)
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.imshow(f_maps[unit_name].T, cmap='jet', origin='lower', extent=limits)
    #ax.plot([0, x_PFR], [0, y_PFR], color='white')
    ax.add_patch(plt.Circle((0, 0), 0.46, color='r', fill=False))
    ax.add_patch(plt.Circle((x_PFR, y_PFR), 0.02, color='white', fill=True))
    ax.text(-0.45, 0.4, "%.2f" % spat_info[unit_name], color='white', fontsize=12)  # spatial info in the corner
    #ax.text(-0.45, -0.45, round(np.rad2deg(f_PFR[i][1]), 1), color='white', fontsize=12)  # spatial info in the corner
    ax.text(x_PFR + 0.05, y_PFR + 0.05, "%.1f" % peak_FRs[unit_name], color='black', fontsize=12)  # peak firing rate
    #ax.set_aspect('equal')
    ax.set_title(unit_name, fontsize=14)

fig.savefig(snakemake.output[0])