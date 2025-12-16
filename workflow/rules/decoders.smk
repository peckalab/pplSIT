import os


rule state_decoder_EV_SU:
    input:
        meta =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        nmap =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_with_PSTH.h5'),
        aeps =os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_EV_SU_metrics.h5')
    output:
        state_d=os.path.join(config['dst_path'], '{animal}', '{session}', 'decoders', 'behav_states.h5')
    script:
        "../scripts/decoders/state_decoder_ev_su.py"


rule tgt_non_tgt_decoder_ev_su:
    input:
        meta =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        nmap =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_with_PSTH.h5'),
        aeps =os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_EV_SU_metrics.h5')
    output:
        state_d=os.path.join(config['dst_path'], '{animal}', '{session}', 'decoders', 'tgt_non-tgt_ev_su.h5')
    script:
        "../scripts/decoders/tgt_non-tgt_decoder_ev_su.py"


rule state_decoder_unit_mx:
    input:
        meta =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        nmap =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5')
    output:
        state_d=os.path.join(config['dst_path'], '{animal}', '{session}', 'decoders', 'behav_states_unit_mx.h5')
    script:
        "../scripts/decoders/state_decoder_unit_mx.py"
