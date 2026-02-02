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
    """Get all npy files from RawWaveforms folders for all sessions of an animal."""
    animal = wildcards.animal
    dst_path = config['dst_path']
    
    # Get all session directories for this animal
    animal_path = os.path.join(dst_path, animal)
    if not os.path.exists(animal_path):
        return []
    
    npy_files = []
    
    # Iterate through all sessions
    for session in os.listdir(animal_path):
        session_path = os.path.join(animal_path, session)
        rawwaveforms_path = os.path.join(session_path, 'kilosort', 'RawWaveforms')
        
        # Check if RawWaveforms folder exists
        if os.path.isdir(rawwaveforms_path):
            # Get all npy files in that folder
            for npy_file in Path(rawwaveforms_path).glob('*.npy'):
                npy_files.append(str(npy_file))
    
    return npy_files

rule extract_raw_waveforms:
    input:
        oebin_file=lambda w: find_structure_oebin(os.path.join(config['src_path'], w.animal, w.session))[0],
        dat_file=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', '{session}.dat'),
        spk_times=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'spike_times.npy'),
        spk_clusters=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'spike_clusters.npy'),
        unit_labels=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'cluster_group.tsv')
    output:
        ready=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'RawWaveforms', 'RawWaveforms.ready')
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/unit_match"
    script:
        "../scripts/unit_match/extract_raw_waveforms.py"

rule unit_match_across_sessions:
    input:
        npy_files=get_rawwaveforms_npy_files
    output:
        match_table=os.path.join(config['prj_path'], 'unit_match', '{animal}', 'MatchTable.csv')
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/unit_match"
    script:
        "../scripts/unit_match/unit_match_across_sessions.py"