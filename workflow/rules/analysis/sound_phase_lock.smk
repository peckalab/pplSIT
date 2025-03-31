import os


rule sound_phase_lock:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
    output:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'sound_phase_lock.h5')
    script:
        "../../scripts/analysis/sound_phase_lock.py"