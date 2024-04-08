import h5py
import numpy as np

from umap import UMAP
from sklearn.manifold import TSNE


def get_ratio_matrix(moseq, tl, win_l=2, step=1, s_rate=100, syl_num=10):
    idxs_srm_tl = np.arange(0, len(tl), int(step*s_rate))
    syl_ratio_mx = np.zeros([len(idxs_srm_tl), syl_num])
    for k, idx in enumerate(idxs_srm_tl):
        curr_syls = moseq[:, 1][idx:idx + int(win_l*s_rate)]  # second column is syllables reindexed
        for j in np.arange(syl_num):
            syl_ratio_mx[k, j] = np.sum(curr_syls == j) / int(win_l*s_rate)

    # roll 1 step to match
    syl_ratio_mx = np.roll(syl_ratio_mx, 1, axis=0)
    syl_ratio_mx[0] = syl_ratio_mx[1]
    
    return syl_ratio_mx, idxs_srm_tl



# parameters for tSNE / UMAP. Should go to yaml?
umap_dists        = [0.1, 0.3, 0.5, 0.7]
perplexities      = [20, 50, 70, 100]
umap_fits         = {}
tsne_fits         = {}
win_l = 2    # in seconds. Window size to compute ratios
step  = 0.5  # in seconds. Step for syllable ratio matrix


# load data
with h5py.File(snakemake.input[1], 'r') as f:
    moseq = np.array(f['moseq'])
    headers = list(f['moseq'].attrs['headers'])

with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    
syl_ratio_mx, idxs_srm_tl = get_ratio_matrix(moseq, tl, win_l=win_l, step=step)

for perp in perplexities:
    tsne = TSNE(n_components=2, perplexity=perp, random_state=0)
    tsne_fit = tsne.fit_transform(syl_ratio_mx)
    tsne_fits[perp] = tsne_fit
    
for dist in umap_dists:
    umap_2d = UMAP(n_components=2, n_neighbors=30, min_dist=dist, random_state=0)
    umap_fit = umap_2d.fit_transform(syl_ratio_mx)
    umap_fits[dist] = umap_fit

with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('syl_ratio_mx', data=syl_ratio_mx)
    f.create_dataset('idxs_srm_tl', data=idxs_srm_tl)
    f.create_group('tSNE')
    for perp in perplexities:
        f['tSNE'].create_dataset(str(perp), data=tsne_fits[perp])
    f.create_group('UMAP')
    for dist in umap_dists:
        f['UMAP'].create_dataset(str(dist), data=umap_fits[dist])
