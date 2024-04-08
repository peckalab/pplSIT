import os


rule w1_w4_tsne_umap:
    input:
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
    output:
        w1_w4_classified=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'W1-W4_tSNE_UMAP.h5')
    script:
        "../../scripts/analysis/w1_w4_tsne_umap.py"