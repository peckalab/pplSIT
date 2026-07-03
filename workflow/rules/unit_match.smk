import os
from pathlib import Path

# find structure.oebin file in sub(sub)folders) of the given paths
def find_structure_oebin(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path {path} does not exist.")
    structure_files = list(path.rglob('structure.oebin'))
    if not structure_files:
        raise FileNotFoundError(f"No structure.oebin found in {path}.")
    return structure_files

def unit_match_methods():
    unit_match_config = config.get('unit_match', {})
    methods = unit_match_config.get('methods', unit_match_config.get('method', 'unitmatchpy'))
    if isinstance(methods, str):
        methods = [methods]
    return [str(method).strip().lower() for method in methods]

def unit_match_conda_env(wildcards=None):
    return config.get('unit_match', {}).get(
        'conda_env',
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/unit_match"
    )

def unit_match_waveform_source():
    return config.get('unit_match', {}).get('waveform_source', 'unitmatch').strip().lower()

def unit_match_waveform_dir(stream_path):
    source = unit_match_waveform_source()
    if source == 'unitmatch':
        return os.path.join(stream_path, 'RawWaveforms')
    if source == 'bombcell':
        return os.path.join(stream_path, 'bombcell', 'RawWaveforms')
    raise ValueError("Unsupported unit_match.waveform_source: " + repr(source))

def get_rawwaveforms_npy_files(wildcards):
    """Get waveform files and rerun triggers from configured stream-level RawWaveforms folders."""
    animal = wildcards.animal
    dst_path = config['dst_path']
    require_phy_folder = config.get('unit_match', {}).get('require_phy_folder', False)

    npy_files = []

    for session in sorted(config.get('session_IDs', [])):
        if session.split('_')[0] != animal:
            continue

        for stream in streams_for_session(config, animal, session):
            stream_path = os.path.join(dst_path, animal, session, 'kilosort', stream)

            rawwaveforms_path = unit_match_waveform_dir(stream_path)
            phy_path = os.path.join(stream_path, '.phy')

            if require_phy_folder and not os.path.isdir(phy_path):
                continue
            if unit_match_waveform_source() == 'bombcell':
                npy_files.append(os.path.join(stream_path, 'bombcell', 'RawWaveforms', 'RawWaveforms.ready'))
            if not os.path.isdir(rawwaveforms_path):
                continue

            for npy_file in Path(rawwaveforms_path).glob('*.npy'):
                npy_files.append(str(npy_file))

            # Include metadata as rerun triggers; the matching script filters
            # this list back down to RawWaveforms/*.npy before running the method.
            for trigger_name in ('spike_clusters.npy', 'cluster_group.tsv', 'probe.json'):
                trigger_path = os.path.join(stream_path, trigger_name)
                if os.path.exists(trigger_path):
                    npy_files.append(trigger_path)
            if unit_match_waveform_source() == 'bombcell':
                for trigger_name in ('bombcell/RawWaveforms/RawWaveforms.ready', 'cluster_bc_unitType.tsv'):
                    trigger_path = os.path.join(stream_path, trigger_name)
                    if os.path.exists(trigger_path):
                        npy_files.append(trigger_path)

    return npy_files

rule extract_raw_waveforms:
    input:
        dat_file=lambda w: stream_dat_path(w.animal, w.session, w.stream),
        spk_times=k_path('{animal}', '{session}', '{stream}', 'spike_times.npy'),
        spk_clusters=k_path('{animal}', '{session}', '{stream}', 'spike_clusters.npy'),
        unit_labels=k_path('{animal}', '{session}', '{stream}', 'cluster_group.tsv'),
        channel_positions=k_path('{animal}', '{session}', '{stream}', 'channel_positions.npy')
    output:
        ready=k_path('{animal}', '{session}', '{stream}', 'RawWaveforms', 'RawWaveforms.ready')
    params:
        spike_width=lambda w: config['unit_match']['spike_width'],
        samples_before=lambda w: config['unit_match']['samples_before'],
        sample_amount=lambda w: config['unit_match']['sample_amount'],
        extract_good_units_only=lambda w: config['unit_match']['extract_good_units_only'],
        KS4_data=lambda w: config['unit_match']['KS4_data']
    conda:
        unit_match_conda_env
    threads:
        int(os.environ.get(
            'RAW_WAVEFORM_THREADS',
            config.get('unit_match', {}).get('raw_waveform_threads', 128)
        ))
    resources:
        io_heavy=int(config.get('unit_match', {}).get('raw_waveform_io_heavy', 1))
    script:
        "../scripts/unit_match/extract_raw_waveforms.py"

rule unit_match_across_sessions:
    input:
        npy_files=get_rawwaveforms_npy_files # includes raw waveforms plus rerun triggers; script filters raw waveforms
    output:
        group_summary=os.path.join(config['prj_path'], 'unit_match', '{animal}', '{method}', 'probe_group_summary.csv')
    conda:
        unit_match_conda_env
    script:
        "../scripts/unit_match/unit_match_across_sessions.py"
