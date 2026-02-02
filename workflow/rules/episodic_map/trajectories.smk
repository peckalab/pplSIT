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


rule EM_drift_PCA:
    input:
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        traj = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5')
    output:
        dPCA = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'drift_PCA.h5')
    script:
        "../../scripts/episodic_map/drift_PCA.py"


rule EM_per_cell_GLM:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        traj = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5'),
        behf = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'behavior_features.h5'),
        dPCA = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'drift_PCA.h5'),
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
    output:
        cglm1 = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'per_cell_GLM.h5'),
        cglm2 = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'per_cell_GLM_sta.h5')
    script:
        "../../scripts/episodic_map/per_cell_GLM.py"


rule EM_ensembles:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        unit = os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        sphl = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'sound_phase_lock.h5'),
        cglm = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'per_cell_GLM_sta.h5'),
    output:
        ens  = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'ensembles.yaml')
    script:
        "../../scripts/episodic_map/ensembles.py"


rule EM_trajectories_filtered:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        behf = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'behavior_features.h5'),
        ens  = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'ensembles.yaml')
    output:
        out  = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories_filtered.h5')
    script:
        "../../scripts/episodic_map/trajectories_filtered.py"


rule EM_masks:
    input:
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        ens  = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'ensembles.yaml')
    output:
        msks = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'masks.h5')
    script:
        "../../scripts/episodic_map/masks.py"


rule EM_decoder_ID_comb:
    input:
        traj = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5')
    output:
        decI = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'decoder_ID_comb.h5')
    script:
        "../../scripts/episodic_map/decoder_ID_comb.py"


rule EM_decoder_ID_mean:
    input:
        traj1 = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5'),
        traj2 = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories_filtered.h5')
    output:
        decM1 = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'decoder_ID_mean.h5'),
        decM2 = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'decoder_ID_mean_filtered.h5')
    script:
        "../../scripts/episodic_map/decoder_ID_mean.py"
