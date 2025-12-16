import os


rule low_D_embeddings:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        segm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5'),
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5')
    output:
        lowD = os.path.join(config['dst_path'], '{animal}', '{session}', 'episodic_map', 'low_D_embeddings.h5')
    script:
        "../scripts/episodic_map/low_D_embeddings.py"