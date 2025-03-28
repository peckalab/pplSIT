import os

rule performance:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    output:
        performance=os.path.join(config['dst_path'], '{animal}', '{session}', 'performance.h5')
    script:
        "../scripts/performance.py"