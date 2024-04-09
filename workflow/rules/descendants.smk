import os

rule descendants:
    input:
        moseq=os.path.join(config['dst_path'], '{animal}', '{session}', 'MoSeq.h5')
    output:
        descendants=os.path.join(config['dst_path'], '{animal}', '{session}', 'descendants.h5')
    script:
        "../scripts/descendants.py"