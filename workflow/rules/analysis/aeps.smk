import os


rule aeps_lfp_metrics:
    input:
        meta      =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        lfp       =os.path.join(config['dst_path'], '{animal}', '{session}', 'LFP', 'lfp.h5'),
        baseline  =os.path.join(config['dst_path'], '{animal}', '{session}', 'LFP', 'baseline.h5'),
        artifacts =os.path.join(config['dst_path'], '{animal}', '{session}', 'LFP', 'artifacts.h5')
    output:
        aeps_lfp_metrics=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_lfp_metrics.h5')
    script:
        "../../scripts/analysis/aeps_lfp_metrics.py"


rule extract_aeps:
    input:
        manual      =os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json'),
        meta        =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        lfp         =os.path.join(config['dst_path'], '{animal}', '{session}', 'LFP', 'lfp.h5'),
        lfp_base    =os.path.join(config['dst_path'], '{animal}', '{session}', 'LFP', 'baseline.h5'),
        artifacts   =os.path.join(config['dst_path'], '{animal}', '{session}', 'LFP', 'artifacts.h5'),
        aeps_metrics=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_lfp_metrics.h5')
    output:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'AEPs.h5')
    script:
        "../../scripts/analysis/aeps_extract.py"
        

rule aeps_ITPC:
    input:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'AEPs.h5'),
        segm=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        ITPC=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_ITPC.h5')
    script:
        "../../scripts/analysis/aeps_ITPC.py"
        

rule aeps_ITPC_plot:
    input:
        ITPC=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_ITPC.h5')
    output:
        ITPC_plots=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_ITPC.pdf')
    script:
        "../../scripts/analysis/aeps_ITPC_plot.py"
        

rule aeps_EV_SU_metrics:
    input:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'AEPs.h5')
    output:
        metrics=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_EV_SU_metrics.h5')
    script:
        "../../scripts/analysis/aeps_EV_SU_metrics.py"
        

rule plot_AEP_profiles:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'AEPs.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'AEP', 'aeps_profiles.pdf')
    script:
        "../../scripts/analysis/aeps_profiles.py"


# rule compute_aep_components:
#     input:
#         aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEPs.h5'),
#         meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
#     output:
#         os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'AEP_components.h5')
#     script:
#         "../../scripts/analysis/aeps_comps.py"
        

# rule plot_AEP_component_maps:
#     input:
#         os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
#         os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'AEP_components.h5')
#     output:
#         os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'aeps_maps.pdf')
#     script:
#         "../../scripts/analysis/aeps_maps.py"
