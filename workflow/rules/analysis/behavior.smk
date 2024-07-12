import os


rule moseq_tsne_umap:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        moseq=os.path.join(config['dst_path'], '{animal}', '{session}', 'MoSeq.h5'),
    output:
        moseq_classified=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'MoSeq_tSNE_UMAP.h5')
    script:
        "../../scripts/analysis/moseq_tsne_umap.py"


rule moseq_w1_w4_GLM:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        moseq=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'MoSeq_tSNE_UMAP.h5'),
    output:
        moseq_w1_w4_GLM=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'MoSeq_W1-W4_GLM.h5')
    script:
        "../../scripts/analysis/moseq_w1_w4_GLM.py"