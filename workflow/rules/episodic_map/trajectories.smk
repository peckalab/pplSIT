import os


rule EM_trajectories:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        behf = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'behavior_features.h5')
    output:
        lowD = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5')
    script:
        "../../scripts/episodic_map/trajectories.py"


rule EM_decoder_ID_comb:
    input:
        traj = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5')
    output:
        decI = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'decoder_ID_comb.h5')
    script:
        "../../scripts/episodic_map/decoder_ID_comb.py"


rule EM_decoder_ID_mean:
    input:
        traj = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5')
    output:
        decM = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'decoder_ID_mean.h5')
    script:
        "../../scripts/episodic_map/decoder_ID_mean.py"