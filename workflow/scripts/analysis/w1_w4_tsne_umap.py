import h5py, os, sys
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from umap import UMAP
from sklearn.manifold import TSNE
from utils.population import activity_at_phase


# parameters for tSNE / UMAP. Should go to yaml?
umap_dists        = [0.1, 0.3, 0.5, 0.7]
perplexities      = [20, 50, 70, 100]
umap_fits         = {}
tsne_fits         = {}


s_path = os.path.dirname(snakemake.input[0])

# Population activity
w_mx = []
for phase in [1, 2, 3, 4]:
    w_pca = activity_at_phase(s_path, phase, do_pca=True)
    w_mx.append(w_pca) # stay in events space
w_mx = np.column_stack(w_mx)

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