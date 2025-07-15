rule vectors_progression:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        act_mxs=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        segm=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'vectors', 'progression.h5')
    script:
        "../../scripts/analysis/vectors_progression.py"

rule vectors_averages:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        act_mxs=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        segm=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'vectors', 'averages.h5')
    script:
        "../../scripts/analysis/vectors_averages.py"