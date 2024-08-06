import numpy as np
from skimage.segmentation import watershed
from utils.behavior import get_idxs_in_patches, get_extent, density_map


def manifold_segmentation(events, tgt_mx, fit, cfg):

    # collecting successful target event indices
    idxs_succ_ev = []
    for tgt_rec in tgt_mx[tgt_mx[:, 4] == 1]:
        idxs_succ_ev += list(np.arange(tgt_rec[0], tgt_rec[1] + 1))
    idxs_succ_ev = np.array(idxs_succ_ev)

    extent = get_extent(fit, margin=cfg['margin'])

    # density map
    d_map  = density_map(fit, extent, sigma=cfg['sigma'], bin_count=cfg['bin_count'])

    # watershed segmentation
    mask = d_map > 0.1*d_map.max()
    fit_labels = watershed(-d_map, mask=mask)

    # using k-means - not that nice
    #     algo = KMeans(n_clusters=10, random_state=0)
    #     algo.fit(fit)
    #     centers = algo.cluster_centers_
    #     ax.scatter(fit[:, 0], fit[:, 1], s=10, c=algo.labels_)
    #     ax.scatter(centers[:, 0], centers[:, 1], c="r", s=20)

    # TGT cluster statistics
    labels = np.unique(fit_labels)[1:]
    fit_area = (fit_labels > 0).sum()
    label_areas = [(fit_labels == x).sum()/fit_area for x in labels]
    clu_tgt_counts = np.zeros(len(labels))
    x_bins = np.linspace(extent[0], extent[1], fit_labels.shape[0]+1)
    y_bins = np.linspace(extent[2], extent[3], fit_labels.shape[1]+1)

    # assigning labels (clusters) to events and collecting label counts for TGT
    labels_ev = np.zeros(len(fit))
    for i, record in enumerate(fit):  # [idxs_succ_ev]
        x_bin_idx = np.where(x_bins > record[0])[0][0] - 1
        y_bin_idx = np.where(y_bins > record[1])[0][0] - 1

        curr_label = fit_labels[x_bin_idx, y_bin_idx]
        labels_ev[i] = curr_label

        if i in idxs_succ_ev:  # collecting label counts for TGT success
            clu_tgt_counts[curr_label-1] += 1
    tgt_stats = clu_tgt_counts/label_areas

    # shuffle TGT and conf intervals
    iter_count = cfg['shuffle_iter_count']
    tgt_stats_shuf = np.zeros([iter_count, len(labels)])
    tgt_mx_succ = tgt_mx[tgt_mx[:, 4] == 1]
    for k in range(iter_count):
        clu_tgt_counts = np.zeros(len(labels))
        
        tgt_mx_succ_shuf = tgt_mx_succ.copy()
        rand_shift = np.random.randint(len(events))
        tgt_mx_succ_shuf[:, 0] = np.mod(tgt_mx_succ[:, 0] + rand_shift, len(events))
        tgt_mx_succ_shuf[:, 1] = np.mod(tgt_mx_succ[:, 1] + rand_shift, len(events))

        idxs_succ_shuf_ev = []
        for tgt_rec in tgt_mx_succ_shuf:
            idxs_succ_shuf_ev += list(np.arange(tgt_rec[0], tgt_rec[1] + 1))

        for record in fit[idxs_succ_shuf_ev]:
            x_bin_idx = np.where(x_bins > record[0])[0][0] - 1
            y_bin_idx = np.where(y_bins > record[1])[0][0] - 1

            curr_label = fit_labels[x_bin_idx, y_bin_idx]
            clu_tgt_counts[curr_label-1] += 1
        tgt_stats_shuf[k] = clu_tgt_counts/label_areas

    # mean and percentiles by cluster
    tgt_stats_shuf_mean = tgt_stats_shuf.mean(axis=0)
    confidence_low  = np.zeros(tgt_stats_shuf.shape[1])
    confidence_high = np.zeros(tgt_stats_shuf.shape[1])
    for k, col in enumerate(tgt_stats_shuf.T):
        confidence_low[k]  = np.percentile(col, 5)
        confidence_high[k] = np.percentile(col, 95)
        
    idxs_tgt_stats_low  = tgt_stats < confidence_high
    idxs_tgt_stats_high = tgt_stats > confidence_high

    # seleted areas by TGT (set non-TGT clusters to 0, leave TGT patches)
    labels_sel = labels[idxs_tgt_stats_high]
    fit_labels_sel = fit_labels.copy()
    for label in labels:
        if not label in labels_sel:
            idxs_x = np.where(fit_labels_sel == label)[0]
            idxs_y = np.where(fit_labels_sel == label)[1]
            for x, y in np.column_stack([idxs_x, idxs_y]):
                fit_labels_sel[x][y] = 0

    idxs_tgt_succ_state_ev = get_idxs_in_patches(fit, fit_labels_sel, extent, bin_count=cfg['bin_count'])

    return d_map, fit_labels, fit_labels_sel, labels_ev, idxs_tgt_succ_state_ev, tgt_stats, tgt_stats_shuf