import os

def conditional_input_isl(wildcards):
    f = os.path.join(config['src_path'], wildcards.animal, wildcards.session, 'islands.csv')
    return f if os.path.exists(f) else None

rule pack:
    input:
        positions = os.path.join(config['src_path'], '{animal}', '{session}', 'positions.csv'),
        events = os.path.join(config['src_path'], '{animal}', '{session}', 'events.csv'),
        sounds = os.path.join(config['src_path'], '{animal}', '{session}', 'sounds.csv'),
        islands = conditional_input_isl,
        cfg = os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.json'),
        manual = os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json')
    output:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    script:
        "../scripts/pack.py"