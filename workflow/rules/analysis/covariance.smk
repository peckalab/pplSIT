rule covariance_matrices_and_plots:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        act_mxs=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        segments=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'covariance', 'covariance.h5'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'covariance', 'covariance.pdf'),
    script:
        "../../scripts/analysis/covariance.py"