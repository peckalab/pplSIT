import os, json
import pandas as pd
import numpy as np


def load_ks_units_before(ks_path):
    # load units from kilosort BEFORE manual curation.
    # uses KSLabel == good to load clusters.
    # path to the folder with kilosorted data

    with open(os.path.join(ks_path, 'probe.json'), 'r') as json_file:  
        probe = json.load(json_file)

    ch_shank_map = np.array(probe['kcoords'], dtype=np.int16) + 1
    shanks = np.unique(ch_shank_map)

    # load kilosorted spike times / clusters / templates / positions / labels
    s_times   = np.load(os.path.join(ks_path, 'spike_times.npy'))  # all spike times of all clusters (1D array)
    s_clust   = np.load(os.path.join(ks_path, 'spike_clusters.npy'))  # IDs of clusters for each spike
    templates = np.load(os.path.join(ks_path, 'templates.npy'))  # cluster (unit), timepoints, channel
    ch_pos    = np.load(os.path.join(ks_path, 'channel_positions.npy'))  # x, y
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


def load_ks_units_after(ks_path):
    # load kilosorted spike times / clusters / templates / positions / labels
    s_times  = np.load(os.path.join(ks_path, 'spike_times.npy'))  # all spike times of all clusters (1D array)
    s_clust  = np.load(os.path.join(ks_path, 'spike_clusters.npy'))  # IDs of clusters for each spike
    clu_info = pd.read_csv(os.path.join(ks_path, 'cluster_info.tsv'), sep='\t')

    all_units = {}
    unit_info = {}
    shanks = clu_info['sh'].unique()
    for shank in shanks:
        clu_info_sh = clu_info[clu_info['sh'] == shank]

        clu_info_sh_good = clu_info_sh[(clu_info_sh['KSLabel'] == 'good') | (clu_info_sh['group'] == 'good')]
        clu_info_sh_good = clu_info_sh_good[clu_info_sh_good['group'] != 'noise']

        spiketrains = {}
        u_info_sh = {}
        for i, record in clu_info_sh_good.iterrows():
            clu_id = record['cluster_id']
            spiketrains[clu_id] = s_times[np.where(s_clust == clu_id)[0]]

            u_info_sh[clu_id] = np.array([
                record['Amplitude'],
                record['ContamPct'],
                record['amp'],
                record['ch'],
                record['depth'],
                record['fr'],
                1 if record['KSLabel'] == 'good' else 2,
            ])
        
        all_units[int(shank)+1] = spiketrains
        unit_info[int(shank)+1] = u_info_sh

    return all_units, unit_info  # Amplitude, ContamPct, amp, channel, depth, FR, label