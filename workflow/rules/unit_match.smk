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

def get_rawwaveforms_npy_files(wildcards):
    """Get waveform files and rerun triggers from stream-level RawWaveforms folders."""
    animal = wildcards.animal
    dst_path = config['dst_path']
    require_phy_folder = config.get('unit_match', {}).get('require_phy_folder', False)
    
    # Get all session directories for this animal
    animal_path = os.path.join(dst_path, animal)
    if not os.path.exists(animal_path):
        return []
    
    npy_files = []
    
    for session in os.listdir(animal_path):
        session_path = os.path.join(animal_path, session)
        kilosort_path = os.path.join(session_path, 'kilosort')
        if not os.path.isdir(kilosort_path):
            continue

        for stream in sorted(os.listdir(kilosort_path)):
            stream_path = os.path.join(kilosort_path, stream)
            if not os.path.isdir(stream_path):
                continue

            rawwaveforms_path = os.path.join(stream_path, 'RawWaveforms')
            phy_path = os.path.join(stream_path, '.phy')

            if require_phy_folder and not os.path.isdir(phy_path):
                continue
            if not os.path.isdir(rawwaveforms_path):
                continue

            for npy_file in Path(rawwaveforms_path).glob('*.npy'):
                npy_files.append(str(npy_file))

            # Include label/spike files as rerun triggers; the matching script filters
            # this list back down to RawWaveforms/*.npy before running UnitMatchPy.
            for trigger_name in ('spike_clusters.npy', 'cluster_group.tsv'):
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
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/unit_match"
    threads: 256
    script:
        "../scripts/unit_match/extract_raw_waveforms.py"

rule unit_match_across_sessions:
    input:
        npy_files=get_rawwaveforms_npy_files # npy files includes raw waveforms, spike_clusters.npy and cluster_group.tsv files (see comment in the function)
    output:
        group_summary=os.path.join(config['prj_path'], 'unit_match', '{animal}', 'probe_group_summary.csv')
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/unit_match"
    script:
        "../scripts/unit_match/unit_match_across_sessions.py"
