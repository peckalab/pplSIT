import os

def conditional_input_isl(wildcards):
    f = os.path.join(config['src_path'], wildcards.animal, wildcards.session, 'islands.csv')
    return f if os.path.exists(f) else None

# this function checks if continuous.dat exists at all, 
# if so the final path of dat in the raw folder is used as input,
# which will trigger its creation in copy.smk (if needed) as well as the adc files if it makes sense
def sync_inputs(wildcards):
    base = os.path.join(config['src_path'], wildcards.animal, wildcards.session)

    for dirpath, _, filenames in os.walk(base):
        if "continuous.dat" in filenames:
            return [
                os.path.join(config['src_path'], '{animal}', '{session}', '{session}.dat'),
                os.path.join(config['src_path'], '{animal}', '{session}', 'timestamps.npy'),
            ]

    return []  # important: return [] if not found

rule pack:
    input:
        positions = os.path.join(config['src_path'], '{animal}', '{session}', 'positions.csv'),
        events = os.path.join(config['src_path'], '{animal}', '{session}', 'events.csv'),
        sounds = os.path.join(config['src_path'], '{animal}', '{session}', 'sounds.csv'),
        islands = conditional_input_isl,
        cfg = os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.json'),
        manual = os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json'),
        
        optional_sync_inputs = sync_inputs,
    output:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    script:
        "../scripts/pack.py"