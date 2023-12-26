import numpy as np
import scipy.ndimage as ndi
from scipy import signal
#from skimage.registration import phase_cross_correlation
#from skimage.transform import warp_polar, rotate


def bins2meters(x, y, xy_range, bin_size=0.02):
    x_in_m = x*bin_size - (xy_range[1] - xy_range[0])/2
    y_in_m = y*bin_size - (xy_range[3] - xy_range[2])/2
    return x_in_m, y_in_m


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def gaussian_kernel_2D(sigma=0.1):
    lin_profile = np.linspace(-10, 10, 50)
    bump = np.exp(-sigma * lin_profile**2)
    bump /= np.trapz(bump)  # normalize to 1
    return bump[:, np.newaxis] * bump[np.newaxis, :]


def place_field_2D(pos, pos_firing, sampling_rate, bin_size=0.02, sigma=0.15, xy_range=None):
    """
    :param pos:             x, y positions sampled at sampling_rate
    :param pos_firing:      positions when spikes occured
    :param sampling_rate:   sampling rate of above in Hz
    :param bin_size:        size of the squared bin to calculate firing / occupancy,
                                same units as in pos, e.g. meters
    :param sigma:           standard deviation, for smoothing
    :param range:           [xmin, xmax, ymin, ymax] - array of X and Y boundaries to limit the map
    :return:
                            occupancy_map, spiking_map, firing_map, s_firing_map
    """
    x_min = xy_range[0] if xy_range is not None else pos[:, 0].min()
    x_max = xy_range[1] if xy_range is not None else pos[:, 0].max()
    y_min = xy_range[2] if xy_range is not None else pos[:, 1].min()
    y_max = xy_range[3] if xy_range is not None else pos[:, 1].max()

    x_range = x_max - x_min
    y_range = y_max - y_min

    y_bin_count = int(np.ceil(y_range / bin_size))
    x_bin_count = int(np.ceil(x_range / bin_size))

    pos_range = np.array([[x_min, x_max], [y_min, y_max]])

    occup, x_edges, y_edges = np.histogram2d(pos[:, 0], pos[:, 1], bins=[x_bin_count, y_bin_count], range=pos_range)
    occupancy_map = occup / sampling_rate

    # spiking map
    spiking_map, xs_edges, ys_edges = np.histogram2d(pos_firing[:, 0], pos_firing[:, 1], bins=[x_bin_count, y_bin_count], range=pos_range)

    # firing = spiking / occupancy
    firing_map = np.divide(spiking_map, occupancy_map, out=np.zeros_like(spiking_map, dtype=float), where=occupancy_map!=0)

    # apply gaussial smoothing
    kernel = gaussian_kernel_2D(sigma)
    s_firing_map  = signal.convolve2d(firing_map, kernel, mode='same')
    occupancy_map = signal.convolve2d(occupancy_map, kernel, mode='same')

    return occupancy_map, spiking_map, firing_map, s_firing_map


def map_stats(f_map, o_map):
    """
    f_map:  2D matrix of unit firing rate map
    o_map:  2D matrix of animal occupancy map
    """
    o_map_norm = o_map / o_map.sum()

    meanrate = np.nansum(np.nansum(np.multiply(f_map, o_map_norm)))
    meansquarerate = np.nansum(np.nansum(np.multiply(f_map ** 2, o_map_norm)))
    maxrate = np.max(f_map)

    sparsity = 0 if meansquarerate == 0 else meanrate**2 / meansquarerate
    selectivity = 0 if meanrate == 0 else maxrate / meanrate

    peak_FR = f_map.max()
    spatial_info = 0  # default value

    if peak_FR > 0:
        f_map_norm = f_map / meanrate
        tmp = np.multiply(o_map_norm, f_map_norm)
        tmp = np.multiply(tmp, np.log2(f_map_norm, where=f_map_norm > 0))
        spatial_info = np.nansum(tmp)

    return sparsity, selectivity, spatial_info, peak_FR


def get_field_patches(f_map, threshold=0.5):
    """
    f_map:      2D matrix of unit firing rate map

    return:     sorted field patches, where the largest field has index 1

    se also: https://exeter-data-analytics.github.io/python-data/skimage.html
    """
    patch_idxs = f_map > threshold*f_map.max()

    # label individual patches
    p_labels, _ = ndi.label(patch_idxs)

    # sort patches according to the patch size
    sort_idxs = (-1*np.bincount(p_labels.flat)[1:]).argsort()
    p_idxs = np.unique(p_labels)[1:]

    p_labels_sorted = np.zeros(p_labels.shape, dtype=np.int32)
    for i, idx in enumerate(sort_idxs):
        p_labels_sorted[p_labels == p_idxs[idx]] = i + 1

    return p_labels_sorted


def best_match_rotation_polar(map_A, map_B):
    """ From 
    https://scikit-image.org/docs/stable/auto_examples/registration/plot_register_rotation.html

    TODO explore log-polar transform
    """
    radius = map_A.shape[0]/2
    
    A_polar = warp_polar(map_A, radius=radius)
    B_polar = warp_polar(map_B, radius=radius)
    
    # this way doesn't work good
    #shifts, error, phasediff = phase_cross_correlation(A_polar, B_polar)
    
    ccf = signal.correlate(A_polar.T, B_polar.T, mode='same').sum(axis=0)
    
    return ccf, -(np.argmax(ccf) - 180)


def best_match_rotation_pearson(map_A, map_B, delta=6):
    if not int(delta) == delta:
        raise ValueError('delta should be of type int')
    
    n = map_A.shape[0]
    a, b, r = n/2, n/2, int(0.93*n/2)
    y,x = np.ogrid[-a:n-a, -b:n-b]
    mask = x*x + y*y <= r*r
    
    angles_num = int(360/delta)
    angles = np.linspace(0, (angles_num - 1) * delta, angles_num)  # in degrees
    
    corrs = np.zeros(len(angles))
    for i, alpha in enumerate(angles):
        corr_2D = np.corrcoef(map_A, ndi.rotate(map_B, alpha, reshape=False))
        #print(corr_2D.shape, mask.shape)
        #corrs[i] = corr_2D[mask].mean()
        corrs[i] = corr_2D.mean()

    # normalize before correlation?
        
    return angles, corrs, np.argmax(corrs)*6  # actual value


def get_positions_relative_to(pos_alo, HD, poi):
    """
    pos_allo    2D array of alloentric animal positions X, Y    (shape Nx2)
    HD          vector of HD angles (in rad.) for each pos_allo (shape N)
    poi         point of interest (x, y) relative to which to compute egocentric coords.
    """
    pos_poi = poi - pos_alo
    R = np.linalg.norm(pos_poi, axis=1)  # distance from animal to the POI
    phi_alo = (np.degrees(np.arctan2(pos_poi[:, 1], pos_poi[:, 0])) - HD) % 360  # angle to POI in allocentric frame, in deg.
    phi_ego = phi_alo - np.rad2deg(HD)  # angle to POI in egocentric frame, in deg.
    phi_ego = np.deg2rad(phi_ego)

    return np.array([np.multiply(R, np.cos(phi_ego)), np.multiply(R, np.sin(phi_ego))]).T