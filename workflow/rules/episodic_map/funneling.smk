import os


rule EM_funneling:
    input:
        traj = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'trajectories.h5')
    output:
        funn = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'funneling.h5')
    script:
        "../../scripts/episodic_map/funneling.py"