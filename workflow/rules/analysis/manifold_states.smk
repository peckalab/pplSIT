import os


# legacy AL / PH state inference
# rule AL_behav_and_neuro_state_idxs:
#     input:
#         meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
#         descendants=os.path.join(config['dst_path'], '{animal}', '{session}', 'descendants.h5'),
#         moseq_tsne=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'MoSeq_tSNE_UMAP.h5'),
#         neuro_tsne=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'W1-W4_tSNE_UMAP.h5')
#     output:
#         state_idxs=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'state_idxs.h5')
#     script:
#         "../../scripts/analysis/state_idxs.py"

rule nMAP_segmentation:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        descendants=os.path.join(config['dst_path'], '{animal}', '{session}', 'descendants.h5'),
        neuro_tsne=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'W1-W4_tSNE_UMAP.h5')
    output:
        nMAP_segmentation=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'nMAP_segmentation.h5')
    script:
        "../../scripts/analysis/nMAP_segmentation.py"

rule bMAP_segmentation:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        descendants=os.path.join(config['dst_path'], '{animal}', '{session}', 'descendants.h5'),
        behav_tsne=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'MoSeq_tSNE_UMAP.h5')
    output:
        nMAP_segmentation=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bMAP_segmentation.h5')
    script:
        "../../scripts/analysis/bMAP_segmentation.py"