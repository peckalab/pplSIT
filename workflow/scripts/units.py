import os, sys
import h5py
import json
import numpy as np
import scipy.ndimage as ndi

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import load_clu_res, XMLHero
from utils.spiketrain import instantaneous_rate, spike_idxs
from utils.hdf import create_dataset, H5NAMES
from utils.spatial import place_field_2D, map_stats, get_field_patches
from utils.spatial import bins2meters, cart2pol


# unit metrics to compute
metric_names = (H5NAMES.o_maps, H5NAMES.f_maps, H5NAMES.sparsity, H5NAMES.selectivity, \
                H5NAMES.spat_info, H5NAMES.peak_FR, H5NAMES.f_patches, H5NAMES.f_COM, \
                H5NAMES.pfr_center, H5NAMES.occ_info, H5NAMES.o_patches, H5NAMES.o_COM)

sorted_data_path = os.path.dirname(snakemake.input[1])

# loading unit data
units = load_clu_res(sorted_data_path)  # spikes are in samples, not seconds
if snakemake.config['units']['source'] == 'neurosuite':
    # neurosuite: read from XML
    xml_files = [f for f in os.listdir(sorted_data_path) if f.find('.xml') > 0]
    if len(xml_files) == 0:
        raise ValueError('Need XML settings file to store spiketrains sorted by Neurosuite')

    neurosuite_settings_file = os.path.join(sorted_data_path, xml_files[0])
    sampling_rate = XMLHero(neurosuite_settings_file).get_sampling_rate()

else:
    # kilosort: read from settings.json
    kilosort_settings_file = os.path.join(sorted_data_path, 'settings.json')
    with open(kilosort_settings_file, 'r') as json_file:
        sampling_rate = json.load(json_file)['fs']


# loading timeline
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])  # time, X, Y, speed, HD, trials, sounds
    run_idxs = np.where(tl[:, 3] > 0.04)[0]

for electrode_idx in units.keys():
    unit_idxs = units[electrode_idx]

    for unit_idx, spiketrain in unit_idxs.items():
        unit_name = '%s-%s' % (electrode_idx, unit_idx)

        # 3 main ways to store a spiketrain
        s_times = spiketrain/sampling_rate  # spike times in seconds
        i_rate  = instantaneous_rate(s_times, tl[:, 0])  # instantaneous rate, sampling to timeline
        s_idxs  = spike_idxs(s_times, tl[:, 0])  # timeline indices of individual spikes

        create_dataset(snakemake.output[0], unit_name, H5NAMES.spike_times, s_times)
        create_dataset(snakemake.output[0], unit_name, H5NAMES.inst_rate, i_rate)
        create_dataset(snakemake.output[0], unit_name, H5NAMES.spike_idxs, s_idxs)

        # spatial metrics
        xy_range = [-0.5, 0.5, -0.5, 0.5]  # make fixed for cross-comparisons
        bin_size = 0.02  # 2 cm
        s_rate_pos = round(1.0 / np.diff(tl[:, 0]).mean())

        # keep only spiking when running? > 4cm/s
        #s_idxs = np.intersect1d(s_idxs, run_idxs)

        # compute 2D maps: occupancy and firing rate (place fields)
        unit_pos = tl[s_idxs][:, 1:3]
        traj_pos = tl[:, 1:3]
        #xy_range = [tl[:, 1].min(), tl[:, 1].max(), tl[:, 2].min(), tl[:, 2].max()]
        o_map, s1_map, s2_map, f_map = place_field_2D(traj_pos, unit_pos, s_rate_pos, bin_size=bin_size, xy_range=xy_range)

        # firing map metrics
        sparsity, selectivity, spat_info, peak_FR = map_stats(f_map, o_map)

        # place field metrics
        patches = get_field_patches(f_map)  # 2D matrix, patches labeled according to the size
        #f_sizes = np.bincount(patches.flat)[1:]  # 1D array of field sizes, sorted
        if f_map.max() == 0:
            f_COM_rho, f_COM_phi, pfr_rho, pfr_phi = 0, 0, 0, 0
        else:
            x_in_b, y_in_b = ndi.center_of_mass(f_map, labels=patches, index=1)  # largest field COM, in bins
            f_COM_rho, f_COM_phi = cart2pol(*bins2meters(x_in_b, y_in_b, xy_range))  # largest field COM, in polar coords.
            x, y = np.where(f_map == np.max(f_map))  # location of the peak unit firing, in bins
            pfr_rho, pfr_phi = cart2pol(*bins2meters(x[0], y[0], xy_range))  # location of the peak unit firing, in polar
        
        # same for occupancy
        _, _, occ_info, _ = map_stats(o_map, o_map)
        o_patches = get_field_patches(o_map)  # 2D matrix, patches labeled according to the size
        x, y = ndi.center_of_mass(o_map, labels=o_patches, index=1)  # largest field COM, in bins
        o_COM_rho, o_COM_phi = cart2pol(*bins2meters(x, y, xy_range))     # largest field COM, in polar coords.

        # iterate over metrics, order should match metric_names defined above
        for i, ds in enumerate([o_map, f_map, sparsity, selectivity, spat_info, peak_FR, \
            patches, (f_COM_rho, f_COM_phi), (pfr_rho, pfr_phi), occ_info, \
            o_patches, (o_COM_rho, o_COM_phi)]):
            create_dataset(snakemake.output[0], unit_name, metric_names[i], np.array(ds))