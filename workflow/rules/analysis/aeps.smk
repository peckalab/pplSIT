

rule compute_aep_components:
    input:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEPs.h5'),
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'AEP_components.h5')
    script:
        "../../scripts/analysis/aeps.py"
        