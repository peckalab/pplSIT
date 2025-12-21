import os


rule EM_trajectories:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5')
    output:
        lowD = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5')
    script:
        "../../scripts/episodic_map/trajectories.py"