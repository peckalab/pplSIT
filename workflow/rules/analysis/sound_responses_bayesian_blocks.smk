import os
from utils.sound_responses.settings import fitting_config, timing_config

# One GPU job at a time, including when --resources is omitted. The backend also
# holds a per-device lock across this user's separate local workflow invocations.
workflow.global_resources.setdefault('sound_responses_gpu', 1)
if workflow.global_resources['sound_responses_gpu'] > 1:
    raise ValueError('sound_responses_gpu must be <= 1 (single-GPU scheduling)')


def sound_responses_environment(wildcards):
    settings = config['sound_responses_bayesian_blocks']
    backend = settings.get('backend', 'cpu')
    if backend not in ('cpu', 'gpu'):
        raise ValueError('sound_responses_bayesian_blocks.backend must be cpu or gpu')
    return settings['gpu_conda_env' if backend == 'gpu' else 'conda_env']


rule fit_sound_responses_bayesian_blocks:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        implementation=[os.path.join(workflow.basedir, 'utils', 'sound_responses', name)
                        for name in ('settings.py', 'segmentation.py', 'fits_session.py', 'metadata.py',
                                     'gpu.py', 'gpu_kernel.cu')]
    output:
        h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'fits.h5')
    params:
        settings=lambda wc: fitting_config(config['sound_responses_bayesian_blocks'])
    threads: 1
    resources:
        sound_responses_gpu=lambda wc: int(config['sound_responses_bayesian_blocks'].get('backend', 'cpu') == 'gpu')
    conda:
        sound_responses_environment
    script:
        '../../scripts/sound_responses/fit_bayesian_blocks.py'


rule sound_responses_bayesian_blocks:
    input:
        fits=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'fits.h5'),
        implementation=[os.path.join(workflow.basedir, 'utils', 'sound_responses', name)
                        for name in ('settings.py', 'bayesian_blocks.py', 'session.py', 'metadata.py')]
    output:
        h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'responses.h5'),
        csv=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'summary.csv')
    params:
        settings=lambda wc: timing_config(config['sound_responses_bayesian_blocks'])
    threads: 1
    conda:
        config['sound_responses_bayesian_blocks']['conda_env']
    script:
        '../../scripts/sound_responses/bayesian_blocks.py'


rule plot_sound_responses_bayesian_blocks:
    input:
        h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'responses.h5'),
        implementation=os.path.join(workflow.basedir, 'utils', 'sound_responses', 'plots.py')
    output:
        pdf=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'figures', 'responses.pdf')
    conda:
        config['sound_responses_bayesian_blocks']['conda_env']
    script:
        '../../scripts/sound_responses/plot_bayesian_blocks.py'


rule plot_sound_responses_psth_micro_overlay:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        responses=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'responses.h5'),
        phase_lock=conditional_input,
        implementation=[os.path.join(workflow.basedir, 'utils', 'sound_responses', 'psth_overlay.py'),
                        os.path.join(workflow.basedir, 'utils', 'psth.py')]
    output:
        pdf=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'figures', 'psth_bgr_tgt_bar.pdf'),
        preview=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'bayesian_blocks', 'figures', 'psth_bgr_tgt_bar_preview.png')
    params:
        micro=lambda wc: config['psth']['micro']
    # Uses the same base environment (including Numba) as the existing psth_micro rule.
    script:
        '../../scripts/sound_responses/plot_psth_micro_overlay.py'
