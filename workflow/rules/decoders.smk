import os


rule state_decoder:
    input:
        meta =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        nmap =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_with_PSTH.h5'),
        aeps =os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_EV_SU_metrics.h5')
    output:
        state_d=os.path.join(config['dst_path'], '{animal}', '{session}', 'decoders', 'behav_states.h5')
    script:
        "../scripts/decoders/state_decoder.py"
