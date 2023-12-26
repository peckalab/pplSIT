import os, sys
import h5py
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.aeps import compute_metric, AEP_metrics_lims, outlier_lims, AEP_metrics_methods



with h5py.File(snakemake.input[1], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])

with h5py.File(snakemake.input[0], 'r') as f:
    areas = [x for x in f]

collected = dict([(area, {}) for area in areas])
for area in areas:
    with h5py.File(snakemake.input[0], 'r') as f:
        ds_name = [x for x in f[area]][0]
        aeps = np.array(f[area][ds_name])

    # remove outliers
    aeps[aeps >  outlier_lims[area]]  = outlier_lims[area]
    aeps[aeps < -outlier_lims[area]] = -outlier_lims[area]

    metrics_raw, metrics_norm  = {}, {}
    for m_name, bounds in AEP_metrics_lims[area].items():
        method = AEP_metrics_methods[area][m_name]
        m_raw, m_norm = compute_metric(aeps, method, bounds[0], bounds[1])
        metrics_raw[m_name]  = m_raw
        metrics_norm[m_name] = m_norm
        
    collected[area]['raw']  = metrics_raw
    collected[area]['norm'] = metrics_norm

# save metrics to file
with h5py.File(snakemake.output[0], 'w') as f:
    for area in areas:
        grp = f.create_group(area)

        for m_name in collected[area]['raw'].keys():
            lims_as_str = ','.join([str(x) for x in AEP_metrics_lims[area][m_name]]) 
            d = grp.create_dataset(m_name + '_raw', data=collected[area]['raw'][m_name])
            d.attrs['limits'] = lims_as_str
            d = grp.create_dataset(m_name + '_norm', data=collected[area]['norm'][m_name])
            d.attrs['limits'] = lims_as_str
            grp.attrs['limits'] = lims_as_str