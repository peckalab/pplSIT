from scipy import signal
import h5py
import numpy as np
import sys, os
from utils.spatial import gaussian_kernel_2D, get_field_patches



def get_extent(fit, margin=5):
    x_range = fit[:, 0].max() - fit[:, 0].min()
    y_range = fit[:, 1].max() - fit[:, 1].min()
    max_range = np.max([x_range, y_range])

    x_min = (fit[:, 0].min() + x_range/2) - max_range*(0.5 + margin/100)
    x_max = (fit[:, 0].min() + x_range/2) + max_range*(0.5 + margin/100)
    y_min = (fit[:, 1].min() + y_range/2) - max_range*(0.5 + margin/100)
    y_max = (fit[:, 1].min() + y_range/2) + max_range*(0.5 + margin/100)

    return x_min, x_max, y_min, y_max


def density_map(fit, extent, sigma=0.4, bin_count=100):
    pos_range = np.array([[extent[0], extent[1]], [extent[2], extent[3]]])

    d_map, xs_edges, ys_edges = np.histogram2d(fit[:, 0], fit[:, 1], bins=[bin_count, bin_count], range=pos_range)

    kernel = gaussian_kernel_2D(sigma)
    return signal.convolve2d(d_map, kernel, mode='same')


def get_idxs_in_patches(fit, patches, extent, bin_count=100):
    x_bins = np.linspace(extent[0], extent[1], bin_count)
    y_bins = np.linspace(extent[2], extent[3], bin_count)

    idxs_in = []
    for i, (x, y) in enumerate(fit):
        x_idx = np.argmin(np.abs(x - x_bins))
        y_idx = np.argmin(np.abs(y - y_bins))
        if patches[x_idx][y_idx] > 0:
            idxs_in.append(i)

    return np.array(idxs_in)


def get_idxs_behav_state(source, session, idxs_tl_sample, fit_type='tSNE', fit_parm=70, sigma=0.3, margin=10, bin_count=100):
    # returns idxs to timeline!
    animal = session.split('_')[0]
    meta_file        = os.path.join(source, animal, session, 'meta.h5')
    moseq_class_file = os.path.join(source, animal, session, 'analysis', 'MoSeq_tSNE_UMAP.h5')

    with h5py.File(meta_file, 'r') as f:
        tl = np.array(f['processed']['timeline'])
        tgt_mx = np.array(f['processed']['target_matrix'])
    with h5py.File(moseq_class_file, 'r') as f:
        idxs_srm_tl = np.array(f['idxs_srm_tl'])
        fit = np.array(f[fit_type][str(fit_parm)])

    idxs_state = np.array([i for i, x in enumerate(idxs_srm_tl) if x in idxs_tl_sample], dtype=np.int32)
    
    extent = get_extent(fit, margin=margin)
    behav_map      = density_map(fit[idxs_state], extent, sigma=sigma, bin_count=bin_count)
    state_patches  = get_field_patches(behav_map, 0.2)
    idxs_srm_state = get_idxs_in_patches(fit, state_patches, extent, bin_count=bin_count)
    
    # convert to timeline idxs
    bins_to_fill = int((idxs_srm_tl[1] - idxs_srm_tl[0])/2)
    idxs_res = []
    for idx in idxs_srm_tl[idxs_srm_state]:
        idxs_res += list(range(idx - bins_to_fill, idx + bins_to_fill))
    idxs_res = np.array(idxs_res, dtype=np.int32)
    idxs_res = idxs_res[idxs_res > 0]
    idxs_res = idxs_res[idxs_res < len(tl) - 1]
        
    return idxs_res
    
    
def get_idxs_neuro_state(source, session, idxs_ev_sample, fit_type='tSNE', fit_parm=70, sigma=0.3, margin=10, bin_count=100):
    # returns idxs in sound events space
    animal = session.split('_')[0]
    meta_file        = os.path.join(source, animal, session, 'meta.h5')
    umap_file = os.path.join(source, animal, session, 'analysis', 'W1-W4_tSNE_UMAP.h5')

    with h5py.File(meta_file, 'r') as f:
        tl = np.array(f['processed']['timeline'])
        tgt_mx = np.array(f['processed']['target_matrix'])
    with h5py.File(umap_file, 'r') as f:
        fit = np.array(f[fit_type][str(fit_parm)])  # already in event sampling

    extent = get_extent(fit, margin=margin)
    selected_map  = density_map(fit[idxs_ev_sample], extent, sigma=sigma, bin_count=bin_count)
    state_patches = get_field_patches(selected_map, 0.3)

    return get_idxs_in_patches(fit, state_patches, extent, bin_count=bin_count)