rule CCR_LFP_pop_conds:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segments=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        resp_aep=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_EV_SU_metrics.h5'),
        resp_pop=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_noconv.h5')
        #resp_pop=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_with_PSTH.h5'),
    output:
        resp_mx=os.path.join(config['dst_path'], '{animal}', '{session}', 'CCR', 'CCR_LFP_pop_conds.h5')
    script:
        "../../scripts/analysis/CCR_LFP_pop_conds.py"