import os


rule EM_episode_similarity:
    input:
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        traj = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5'),
        behf = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'behavior_features.h5')
    output:
        epsim = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'episode_similarity.h5')
    script:
        "../../scripts/episodic_map/episode_similarity.py"


rule EM_regression:
    input:
        epsim = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'episode_similarity.h5')
    output:
        reg_latent = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'regression_latent.h5'),
        reg_pv = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'regression_pv.h5'),
    script:
        "../../scripts/episodic_map/regression.py"