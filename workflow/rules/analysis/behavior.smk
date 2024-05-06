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


rule AL_behav_and_neuro_state_idxs:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        descendants=os.path.join(config['dst_path'], '{animal}', '{session}', 'descendants.h5'),
        moseq_tsne=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'MoSeq_tSNE_UMAP.h5'),
        neuro_tsne=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'W1-W4_tSNE_UMAP.h5')
    output:
        state_idxs=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'state_idxs.h5')
    script:
        "../../scripts/analysis/state_idxs.py"