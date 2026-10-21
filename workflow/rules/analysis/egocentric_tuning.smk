from utils.egocentric_tuning.core import settings as tuning_settings


rule egocentric_tuning:
    input:
        geometry=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'geometry', 'geometry.h5'),
        responses=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'responses.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        implementation=[os.path.join(workflow.basedir, 'utils', 'egocentric_tuning', name)
                        for name in ('core.py', 'session.py', 'low_spike_thresholds.json')]
    output:
        h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'egocentric', 'tuning.h5'),
        csv=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'egocentric', 'summary.csv'),
        comparisons=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'egocentric', 'comparisons.csv')
    params:
        settings=lambda wc: tuning_settings(config['egocentric_tuning'])
    threads: 1
    script:
        '../../scripts/egocentric_tuning/tuning.py'


rule plot_egocentric_tuning:
    input:
        h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'egocentric', 'tuning.h5'),
        implementation=[os.path.join(workflow.basedir, 'utils', 'egocentric_tuning', name)
                        for name in ('plots.py', 'core.py', 'low_spike_thresholds.json')]
    output:
        pdf=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'egocentric', 'figures', 'tuning.pdf'),
        preview=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'egocentric', 'figures', 'tuning.png')
    script:
        '../../scripts/egocentric_tuning/plot.py'
