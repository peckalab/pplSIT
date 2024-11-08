import os


rule coherence:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
    output:
        coh=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'coherence.h5')
    script:
        "../../scripts/analysis/coherence.py"