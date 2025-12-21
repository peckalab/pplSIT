import os


rule EM_behavior_features:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm = os.path.join(config['dst_path'], '{animal}', '{session}', 'DLC_100Hz.h5')
    output:
        lowD = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'behavior_features.h5')
    script:
        "../../scripts/episodic_map/behavior_features.py"