import os


rule w1_w4_tsne_umap:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
    output:
        w1_w4_classified=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'W1-W4_tSNE_UMAP.h5')
    script:
        "../../scripts/analysis/w1_w4_tsne_umap.py"


rule population_matrices:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
    output:
        resp_mx=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'population_matrices.h5')
    script:
        "../../scripts/analysis/population_matrices.py"


rule nMAP_EV_SU_no_convolution:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        psths=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_micro.h5'),
        segm=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        nMAP_EV_SU=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_noconv.h5')
    script:
        "../../scripts/analysis/nMAP_EV_SU_noconv.py"


rule nMAP_EV_SU_with_PSTH:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        resp_mx=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'population_matrices.h5'),
        psths=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_micro.h5')
    output:
        nMAP_EV_SU=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_EV_SU_with_PSTH.h5')
    script:
        "../../scripts/analysis/nMAP_EV_SU_with_PSTH.py"