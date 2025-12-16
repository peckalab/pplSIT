import os, sys
import h5py
import json
import numpy as np

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import umap  # umap-learn


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.low_D_embeddings import *

#aep_type    = snakemake.config['state_decoder']['aep_type']
#kernel_size = snakemake.config['state_decoder']['kernel_size']

cfg = snakemake.config['low_D_embeddings']

h5_path = snakemake.output[0]
target_label = 0  # by default: STATE_KEYS order puts target first
X_is_neurons_by_time = True
keep_only_labeled = True

transform =        cfg['transform']
zscore =           cfg['zscore']
pca_n_components = cfg['pca_n_components']
tsne_seed =        cfg['tsne_seed']
tsne_perplexity =  cfg['tsne_perplexity']
k_list =           cfg['k_list']
n_perm =           cfg['n_perm']
umap_seeds =       cfg['umap_seeds']
umap_params = dict(
    n_neighbors=cfg['umap_n_neighbors'],
    min_dist=cfg['umap_min_dist'],
    metric=cfg['umap_metric']
)

meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
actm_file = snakemake.input[2]

STATE_KEYS = [
    "idxs_tgt_sta_succ",
    "idxs_bgr_sta",
    "idxs_bgr_run",
    "idxs_sil_sta",
    "idxs_sil_run",
]

STATE_NAMES  = ["tgt", "bgr_sta", "bgr_run", "sil_sta", "sil_run"]
STATE_LABELS = ['Target', 'BGR/sta', 'BGR/run', 'No stim./sta', 'No stim./run']
STATE_COLORS = ['tab:orange', 'mediumblue', 'deepskyblue', 'black', 'grey']

# -----------------------------
# Reading datasets
# -----------------------------

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

with h5py.File(actm_file, 'r') as f:
    X_counts = np.array(f['mx_250ms']['mx'])
    t_edges = np.array(f['mx_250ms']['bins'])
    
with h5py.File(segm_file, 'r') as f:
    state_idxs = {}
    for idxs_name in STATE_KEYS:
        state_idxs[idxs_name] = np.array(f[idxs_name]).astype(np.int32)

# -----------------------------
# Main pipeline
# -----------------------------

# Labels on original time grid
t_bins = (X_counts.shape[1] if X_is_neurons_by_time else X_counts.shape[0])
labels_full, overlaps = build_state_labels(t_bins, state_idxs)

# Preprocess X and remove any non-finite rows (rare)
Xz, stats, finite_mask = preprocess_counts(
    X_counts, X_is_neurons_by_time=X_is_neurons_by_time, transform=transform, zscore=zscore
)
# Map labels to finite-only rows
labels_fin = labels_full[finite_mask]
t_edges_fin = t_edges  # edges still refer to original grid; store original too.

# Keep only labeled bins if requested
if keep_only_labeled:
    keep_mask = labels_fin >= 0
else:
    keep_mask = np.ones_like(labels_fin, dtype=bool)

X_use = Xz[keep_mask]
labels_use = labels_fin[keep_mask].astype(np.int16)

# PCA (fit in the space used for embeddings)
pca_n = min(pca_n_components, X_use.shape[1])
pca = PCA(n_components=pca_n, random_state=0)
Zp = pca.fit_transform(X_use)  # (n_samples, pca_n)

Y_pca2 = Zp[:, :2].astype(np.float32)

# UMAP embeddings
umap_embeddings = {}
for seed in umap_seeds:
    reducer = umap.UMAP(
        n_components=2,
        random_state=int(seed),
        **umap_params
    )
    umap_embeddings[int(seed)] = reducer.fit_transform(X_use).astype(np.float32)

# t-SNE embedding (2D)
tsne = TSNE(
    n_components=2,
    perplexity=float(tsne_perplexity),
    init="pca",
    learning_rate="auto",
    random_state=int(tsne_seed),
)
Y_tsne2 = tsne.fit_transform(X_use).astype(np.float32)

# kNN purity stats vs circular shift null (do it in 2D embedding spaces)
knn = {
    "pca2": compute_knn_stats(
        Y_pca2, labels_use, target_label=target_label, k_list=k_list, n_perm=n_perm, seed=0
    ),
    "tsne2": compute_knn_stats(
        Y_tsne2, labels_use, target_label=target_label, k_list=k_list, n_perm=n_perm, seed=0
    ),
    "umap2": {}
}
for seed, Y_umap2 in umap_embeddings.items():
    knn["umap2"][str(seed)] = compute_knn_stats(
        Y_umap2, labels_use, target_label=target_label, k_list=k_list, n_perm=n_perm, seed=0
    )

# -----------------------------
# Write to HDF5
# -----------------------------

with h5py.File(h5_path, "w") as f:
    # metadata
    meta = {
        "STATE_KEYS": STATE_KEYS,
        "STATE_NAMES": STATE_NAMES,
        "keep_only_labeled": bool(keep_only_labeled),
        "transform": transform,
        "zscore": bool(zscore),
        "pca_n_components": int(pca_n),
        "umap_params": umap_params,
        "umap_seeds": list(map(int, umap_seeds)),
        "tsne_seed": int(tsne_seed),
        "tsne_perplexity": float(tsne_perplexity),
        "target_label": int(target_label),
        "k_list": list(map(int, k_list)),
        "n_perm": int(n_perm),
        "overlaps": overlaps,
    }
    f.attrs["meta_json"] = json.dumps(meta)

    # inputs
    g_in = f.create_group("inputs")
    g_in.create_dataset("t_edges", data=np.asarray(t_edges, dtype=np.float64), compression="gzip")
    g_in.create_dataset("labels_full", data=labels_full.astype(np.int16), compression="gzip")
    g_in.create_dataset("finite_mask", data=finite_mask.astype(np.uint8), compression="gzip")
    g_in.create_dataset("keep_mask", data=keep_mask.astype(np.uint8), compression="gzip")
    g_in.create_dataset("labels_use", data=labels_use.astype(np.int16), compression="gzip")

    # preprocessing stats
    g_pp = f.create_group("preproc")
    if "mu" in stats:
        g_pp.create_dataset("mu", data=stats["mu"], compression="gzip")
        g_pp.create_dataset("sd", data=stats["sd"], compression="gzip")

    # PCA
    g_pca = f.create_group("pca")
    g_pca.create_dataset("scores", data=Zp.astype(np.float32), compression="gzip")  # (n_samples, n_pc)
    g_pca.create_dataset("embedding2d", data=Y_pca2, compression="gzip")
    g_pca.create_dataset("components", data=pca.components_.astype(np.float32), compression="gzip")  # (n_pc, n_neurons)
    g_pca.create_dataset("explained_variance", data=pca.explained_variance_.astype(np.float32), compression="gzip")
    g_pca.create_dataset("explained_variance_ratio", data=pca.explained_variance_ratio_.astype(np.float32), compression="gzip")
    g_pca.create_dataset("mean", data=pca.mean_.astype(np.float32), compression="gzip")

    # UMAP
    g_umap = f.create_group("umap")
    for seed, emb in umap_embeddings.items():
        gg = g_umap.create_group(f"seed_{seed}")
        gg.create_dataset("embedding2d", data=emb, compression="gzip")

    # t-SNE
    g_tsne = f.create_group("tsne")
    g_tsne.create_dataset("embedding2d", data=Y_tsne2, compression="gzip")

    # kNN stats
    g_knn = f.create_group("knn")

    def write_knn_block(parent, block_dict):
        """
        block_dict is dict: k -> {obs, null_mean, null_std, z, p_one_sided, null_vals}
        """
        for k, d in block_dict.items():
            kg = parent.create_group(f"k_{k}")
            kg.create_dataset("obs", data=np.float32(d["obs"]))
            kg.create_dataset("null_mean", data=np.float32(d["null_mean"]))
            kg.create_dataset("null_std", data=np.float32(d["null_std"]))
            kg.create_dataset("z", data=np.float32(d["z"]))
            kg.create_dataset("p_one_sided", data=np.float32(d["p_one_sided"]))
            kg.create_dataset("null_vals", data=d["null_vals"], compression="gzip")

    # PCA2 and tSNE2
    write_knn_block(g_knn.create_group("pca2"), knn["pca2"])
    write_knn_block(g_knn.create_group("tsne2"), knn["tsne2"])

    # UMAP2 per seed
    g_umap_knn = g_knn.create_group("umap2")
    for seed_str, block in knn["umap2"].items():
        write_knn_block(g_umap_knn.create_group(f"seed_{seed_str}"), block)
