import os, json
import pandas as pd
import numpy as np


def load_ks_units(ks_path):
    # path to the folder with kilosorted data

    with open(os.path.join(ks_path, 'probe.json'), 'r') as json_file:  
        probe = json.load(json_file)

    ch_shank_map = np.array(probe['kcoords'], dtype=np.int16) + 1
    shanks = np.unique(ch_shank_map)

    # load kilosorted spike times / clusters / templates / positions / labels
    s_times   = np.load(os.path.join(ks_path, 'spike_times.npy'))  # all spike times of all clusters (1D array)
    s_clust   = np.load(os.path.join(ks_path, 'spike_clusters.npy'))  # IDs of clusters for each spike
    templates = np.load(os.path.join(ks_path, 'templates.npy'))  # cluster (unit), timepoints, channel
    ch_pos    = np.load(os.path.join(ks_path, 'channel_positions.npy'))  # cluster (unit), timepoints, channel
    ks_labels = pd.read_csv(os.path.join(ks_path, 'cluster_KSLabel.tsv'), sep='\t', header=0)  # cluster, good / mua

    template_maxchans = np.abs(templates).max(axis=1).argmax(axis=1)  # channel with highest AP amplitude for each unit
    clu_ch_mapping    = ch_shank_map[template_maxchans]
    clu_pos_mapping   = ch_pos[template_maxchans]

    # 'good' units
    good_idxs = np.where(ks_labels['KSLabel'] == 'good')[0]

    all_units = {}
    all_pos = {}
    for shank in shanks:
        clu_idxs = np.where(clu_ch_mapping == shank)[0]
        sel_clusters = np.intersect1d(good_idxs, clu_idxs)
        
        spiketrains = {}
        sel_pos = {}
        for clu_id in sel_clusters:
            spiketrains[clu_id] = s_times[np.where(s_clust == clu_id)[0]]
            sel_pos[clu_id] = clu_pos_mapping[clu_id]
            
        all_units[shank] = spiketrains
        all_pos[shank] = sel_pos

    return all_units, all_pos
