import os

rule performance:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    output:
        performance=os.path.join(config['dst_path'], '{animal}', '{session}', 'performance', 'performance.h5'),
        performance_figure=os.path.join(config['dst_path'], '{animal}', '{session}', 'performance', 'performance.pdf'),
        session_metrics_figure=os.path.join(config['dst_path'], '{animal}', '{session}', 'performance', 'session_metrics.png')
    script:
        "../scripts/performance.py"