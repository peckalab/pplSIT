

rule psth_bootstrap_profiles:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        segms=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        psths=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_micro.h5')
    threads: 64
    script:
        "../../scripts/analysis/psth_bootstrap_profiles.py"

rule psth_bootstrap_plots:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        psths=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_micro.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_tgt_line.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sil_line.pdf'),
    threads: 64
    script:
        "../../scripts/analysis/psth_bootstrap_plots.py"