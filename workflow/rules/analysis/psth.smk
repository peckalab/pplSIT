

rule psth_micro_profiles:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        psths=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_micro.h5')
    script:
        "../../scripts/analysis/psth_micro_profiles.py"

rule psth_micro_plots:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        psths=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_micro.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_tgt_line.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sil_line.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_tgt_bar.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sil_bar.pdf')
    script:
        "../../scripts/analysis/psth_micro_plots.py"

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