rule segments:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    script:
        "../../scripts/analysis/segments.py"