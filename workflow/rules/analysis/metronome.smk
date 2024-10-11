
rule metronome:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        periods=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'metronome', 'periods.h5'),
        psth=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'metronome', 'psth.h5'),
        shuffle=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'metronome', 'shuffle.h5')
    script:
        "../../scripts/analysis/metronome.py"