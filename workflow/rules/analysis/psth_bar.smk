

rule psth_micro:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_tgt_bar.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sil_bar.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_tgt_sta_bar.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sta_run_bar.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_sil_sta_run_bar.pdf')
        #os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_distractors.pdf')
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
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_distractors_onset.pdf')
        #os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_noise_offset.pdf')
    script:
        "../../scripts/analysis/psth_macro.py"


rule psth_micro_state:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        ensembles=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'ensembles.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sta_AL_PH.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bgr_sta_AL_tgt.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_tgt_AL_PH.pdf'),
    script:
        "../../scripts/analysis/psth_micro_state.py"


rule psth_micro_bnMAPs:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        desc=os.path.join(config['dst_path'], '{animal}', '{session}', 'descendants.h5'),
        #nMAP=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_with_PSTH.h5'),
        bMAP=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bMAP_segmentation.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bnMAPs_sil_sta_bU_sta_bE.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bnMAPs_sil_sta_bU_run_bU.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bnMAPs_bgr_sta_bU_sta_bE.pdf'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_bnMAPs_bgr_sta_bU_run_bU.pdf')
    script:
        "../../scripts/analysis/psth_micro_bnMAPs.py"