from utils.egocentric_geometry.core import settings as geometry_settings


def geometry_config():
    return geometry_settings(config['egocentric_geometry'])


def geometry_moseq_input(wildcards):
    settings = geometry_config()
    if 'moseq_heading' in settings['references']:
        return [os.path.join(config['dst_path'], wildcards.animal, wildcards.session, 'MoSeq.h5')]
    return []


rule egocentric_geometry:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        dlc=os.path.join(config['dst_path'], '{animal}', '{session}', 'DLC_100Hz.h5'),
        moseq=geometry_moseq_input,
        implementation=[os.path.join(workflow.basedir, 'utils', name) for name in (
            'egocentric_geometry/core.py', 'egocentric_geometry/session.py', 'setup.py', 'session.py')]
    output:
        h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'geometry', 'geometry.h5'),
        csv=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'geometry', 'summary.csv')
    params:
        settings=lambda wc: geometry_config()
    threads: 1
    script:
        '../../scripts/egocentric_geometry/geometry.py'


rule plot_egocentric_geometry:
    input:
        h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'geometry', 'geometry.h5'),
        implementation=[os.path.join(workflow.basedir, 'utils', 'egocentric_geometry', name)
                        for name in ('plots.py', 'core.py')]
    output:
        pdf=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'geometry', 'figures', 'coverage.pdf'),
        preview=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spatial_tuning', 'geometry', 'figures', 'coverage.png')
    script:
        '../../scripts/egocentric_geometry/plot.py'
