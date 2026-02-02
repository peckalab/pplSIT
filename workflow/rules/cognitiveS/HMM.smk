import os


rule cognitiveS_HMM:
    input:
        meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        actm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'activity_matrices.h5'),
        covm = os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'covariance', 'covariance.h5')
    output:
        HMM  = os.path.join(config['dst_path'], '{animal}', '{session}', 'cognitiveS', 'HMM.h5')
    script:
        "../../scripts/cognitiveS/HMM.py"
