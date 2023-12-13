import os

rule pack:
    input:
        os.path.join(config['src_path'], '{animal}', '{session}', 'positions.csv'),
        os.path.join(config['src_path'], '{animal}', '{session}', 'events.csv'),
        os.path.join(config['src_path'], '{animal}', '{session}', 'sounds.csv'),
        os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.json'),
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    script:
        "../scripts/pack.py"