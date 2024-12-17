import os, json
import numpy as np


dst_path = os.path.dirname(snakemake.input[1])

# load probe configuration
with open(snakemake.input[0]) as json_file:  
    probe = json.load(json_file)

channels = probe['chanMap']
ch_shank_map = np.array(probe['kcoords'], dtype=np.int16) + 1
shanks = np.unique(ch_shank_map)

# load kilosorted spike times / clusters
s_times   = np.load(snakemake.input[1])  # all spike times of all clusters (1D array)
s_clust   = np.load(snakemake.input[2])  # IDs of clusters for each spike
templates = np.load(snakemake.input[3])  # cluster (unit), timepoints, channel

# assign clusters to shanks using templates
clu_mapping = []  # 1D, shank id for each cluster
for template in templates:  # template IDs should match cluster IDs
    amplitudes = [t.max() - t.min() for t in template.T]
    best_channel = np.argmax(amplitudes) #+ channel_map[0][0]  # 0-based indexing, adjust to starting channel
    clu_mapping.append(ch_shank_map[np.where(channels == best_channel)[0][0]])
clu_mapping = np.array(clu_mapping)

# split all spikes / cluster IDs into diff shanks
clu_per_shank, spk_per_shank = [], []
for shank in shanks:
    clu_ids = np.where(clu_mapping == shank)[0]
    bool_idxs = np.isin(s_clust, clu_ids)
    clu_per_shank.append(s_clust[bool_idxs])
    spk_per_shank.append(s_times[bool_idxs])

# create .res files with spike times in samples
for i, data in enumerate(spk_per_shank):
    np.savetxt(os.path.join(dst_path, 'NS.res.%d' % shanks[i]), np.array([data]).T, fmt='%d')

# create .clu files
for i, data in enumerate(clu_per_shank):
    f_path = os.path.join(dst_path, 'NS.clu.%d' % shanks[i])
    np.savetxt(f_path, np.array([data]).T, fmt='%d', header=str(len(np.unique(data))), comments='')

# fake .fet files with FAKE features
spread = 500
for i, spk_times in enumerate(spk_per_shank):
    f_path = os.path.join(dst_path, 'NS.fet.%d' % shanks[i])
    fet_count = len(np.where(ch_shank_map == i+1)[0]) * 3 + 1
    fet_mx = np.zeros([len(spk_times), fet_count])
    
    for j, spk_time in enumerate(spk_times):
        #spread = 5000 if clu_per_shank[i][j] == clu and j < 100 else 500
        fake_fet = np.random.randint(spread*2, size=fet_count - 1) - spread
        fake_fet = np.concatenate([fake_fet, [spk_time]])
        fet_mx[j] = fake_fet
    
    clu = clu_per_shank[i][0]  # some random cluster, induce larger spread for better visualization
    idxs = np.where(clu_per_shank[i] == clu)[0][:2]
    fet_mx[idxs[0]][:-1] = np.random.randint(spread*2*10, size=fet_count - 1) - spread*10
    fet_mx[idxs[1]][:-1] = np.random.randint(spread*2*10, size=fet_count - 1) - spread*10

    np.savetxt(f_path, fet_mx, fmt='%d', delimiter=' ', header=str(fet_count), comments='')