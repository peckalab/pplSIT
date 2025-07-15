rule ensembles:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'ensembles.h5')
    script:
        "../../scripts/analysis/ensembles.py"

