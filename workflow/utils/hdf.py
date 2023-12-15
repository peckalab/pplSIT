import h5py

class H5NAMES:
    inst_rate =   {'name': 'inst_rate', 'dims': ['instantaneous firing rate at ~100Hz']}
    spike_times = {'name': 'spike_times', 'dims': ['spike times in seconds']}
    spike_idxs =  {'name': 'spike_idxs', 'dims': ['indices to timeline when spikes occured']}
    mfr =         {'name': 'mean_firing_rate', 'dims': ['epochs: original, conflict, control and all']}
    isi_cv =      {'name': 'isi_coeff_var', 'dims': ['epochs: original, conflict, control and all']}
    isi_fano =    {'name': 'isi_fano_factor', 'dims': ['epochs: original, conflict, control and all']}
    o_maps =      {'name': 'occupancy_maps', 'dims': [
        'epochs: original, conflict, control and all', 'X, bins', 'Y, bins'
    ]}
    f_maps =      {'name': 'firing_rate_maps', 'dims': [
        'epochs: original, conflict, control and all', 'X, bins', 'Y, bins'
    ]}
    sparsity =    {'name': 'sparsity', 'dims': ['epochs: original, conflict, control and all']}
    selectivity = {'name': 'selectivity', 'dims': ['epochs: original, conflict, control and all']}
    spat_info =   {'name': 'spatial_information', 'dims': ['epochs: original, conflict, control and all']}
    peak_FR =     {'name': 'peak_firing_rate', 'dims': ['epochs: original, conflict, control and all']}
    f_patches =   {'name': 'field_patches', 'dims': [
        'epochs: original, conflict, control and all', 'X, bins', 'Y, bins'
    ]}
    f_sizes =     {'name': 'field_sizes', 'dims': ['epochs: original, conflict, control and all']}
    f_COM =       {'name': 'field_center_of_mass', 'dims': ['epochs: original, conflict, control and all', 'rho, phi in polar coords.']}
    pfr_center =  {'name': 'field_center_of_firing', 'dims': ['epochs: original, conflict, control and all', 'rho, phi in polar coords.']}
    occ_info =    {'name': 'occupancy_information', 'dims': ['epochs: original, conflict, control and all']}
    o_patches =   {'name': 'occupancy_patches', 'dims': [
        'epochs: original, conflict, control and all', 'X, bins', 'Y, bins'
    ]}
    o_COM =       {'name': 'occupancy_center_of_mass', 'dims': ['epochs: original, conflict, control and all', 'rho, phi in polar coords.']}
    best_m_rot  = {'name': 'best_match_rotation', 'dims': ['match between: A-B, B-C, A-C', 'correlation profile']}


def create_dataset(h5name, where, descriptor, dataset):
    """
    h5name       path to an HDF5 file
    where        path inside the file
    descriptor   H5NAMES style descriptor of the dataset
    dataset      numpy array to store
    """
    with h5py.File(h5name, 'a') as f:
        if not where in f:
            f.create_group(where)
            
        target_group = f[where]

        if descriptor['name'] in target_group:  # overwrite mode
            del target_group[descriptor['name']]
            
        ds = target_group.create_dataset(descriptor['name'], data=dataset)
        for i, dim in enumerate(descriptor['dims']):
            ds.attrs['dim%s' % i] = dim