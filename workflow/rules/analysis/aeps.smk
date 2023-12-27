import os


rule extract_aeps:
    input:
        manual=os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json'),
        lfp=os.path.join(config['dst_path'], '{animal}', '{session}', 'lfp.h5'),
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
    output:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEPs.h5')
    script:
        "../../scripts/analysis/aeps_extract.py"
        

rule compute_aep_components:
    input:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEPs.h5'),
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'AEP_components.h5')
    script:
        "../../scripts/analysis/aeps_comps.py"
        

rule plot_AEP_component_maps:
    input:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'AEP_components.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'aeps_maps.pdf')
    script:
        "../../scripts/analysis/aeps_maps.py"