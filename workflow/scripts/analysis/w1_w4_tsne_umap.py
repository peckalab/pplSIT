import h5py, os, sys, json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from umap import UMAP
from sklearn.manifold import TSNE
from utils.population import activity_at_phase, pop_activity_phase_shifted


# parameters for tSNE / UMAP. Should go to yaml?
umap_dists        = [0.1, 0.3, 0.5, 0.7]
perplexities      = [20, 50, 70, 100]
umap_fits         = {}
tsne_fits         = {}

s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
animal  = session.split('_')[0]

with h5py.File(snakemake.input[0], 'r') as f:
    cfg = json.loads(f['processed'].attrs['parameters'])

# Population activity - old way
# w_mx = []
# for phase in [1, 2, 3, 4]:
#     w_pca = activity_at_phase(s_path, phase, do_pca=True)
#     w_mx.append(w_pca) # stay in events space
# w_mx = np.column_stack(w_mx)

# select binning for target / backgound
binning_mPFC_AEP = {
    50:  [8, 25, 45, 85, 105],
    75:  [8, 25, 45, 105, 127],
    100: [8, 25, 45, 128, 155],
}

binning_mPFC_PCA = {
    50:  [8, 16, 65, 85],
    75:  [8, 16, 90, 110],
    100: [8, 16, 115, 135],
}

binning_AC_AEP = {
    50:  [15, 28, 73, 100],
    90:  [15, 28, 107, 135],
    100: [15, 28, 118, 145],
    150: [15, 28, 168, 195],
}

# for the moment define binning scheme manually for each animal
all_binnings = {
    '009265': binning_AC_AEP,
    '009266': binning_AC_AEP,
    '57': binning_mPFC_PCA,
}
selected_binning = all_binnings[animal] if animal in all_binnings else binning_AC_AEP

# define actual binning for target / background
dur_bgr = int(cfg['sound']['sounds']['background']['duration'] * 1000)  # in ms
dur_tgt = int(cfg['sound']['sounds']['target']['duration'] * 1000)  # in ms

assert dur_bgr in selected_binning
assert dur_tgt in selected_binning

# electrodes in A1
electrodes = snakemake.config['nMAP_electrodes'][animal]

w_mx = pop_activity_phase_shifted(s_path, selected_binning[dur_bgr], selected_binning[dur_tgt], electrodes=electrodes, do_pca=True)

# reduce 4-to-2 dimensions
for perp in perplexities:
    tsne = TSNE(n_components=2, perplexity=perp, random_state=0)
    tsne_fit = tsne.fit_transform(w_mx)
    tsne_fits[perp] = tsne_fit
    
for dist in umap_dists:
    umap_2d = UMAP(n_components=2, n_neighbors=30, min_dist=dist, random_state=0)
    umap_fit = umap_2d.fit_transform(w_mx)
    umap_fits[dist] = umap_fit

with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('w_mx', data=w_mx)
    f.create_group('tSNE')
    for perp in perplexities:
        f['tSNE'].create_dataset(str(perp), data=tsne_fits[perp])
    f.create_group('UMAP')
    for dist in umap_dists:
        f['UMAP'].create_dataset(str(dist), data=umap_fits[dist])