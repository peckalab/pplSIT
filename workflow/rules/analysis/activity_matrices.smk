rule activity_matrices:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
    threads:
        96
    script:
        "../../scripts/analysis/activity_matrices.py"