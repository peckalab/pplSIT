

rule psth_micro:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_tgt.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sil.pdf')
    script:
        "../../scripts/analysis/psth_micro.py"

rule psth_macro:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_tgt_onset.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_tgt_offset.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_trial_onset.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_noise_offset.pdf')
    script:
        "../../scripts/analysis/psth_macro.py"