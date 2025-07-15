import scipy.cluster.hierarchy as sch
import numpy as np


def cluster_corr(corr_array, threshold=None, inplace=False, n_clusters=None):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    corr_array : a NxN correlation matrix with the columns and rows rearranged
    linkage    : linkage of distances
    labels     : cluster labels
    idx        : sorted incides for original labels
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2 if threshold is None else threshold
    if n_clusters:
        labels = sch.fcluster(linkage, t=n_clusters, criterion='maxclust')
    else:
        labels = sch.fcluster(linkage, cluster_distance_threshold, criterion='distance')
    idx = np.argsort(labels)
    
    if not inplace:
        corr_array = corr_array.copy()
    
    return corr_array[idx, :][:, idx], linkage, labels, idx