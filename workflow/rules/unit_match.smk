import os
from pathlib import Path
from utils.unitmatch import update_cluster_group_with_KSlabel

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
    require_phy_folder = config.get('unit_match', {}).get('require_phy_folder', False)
    
    # Get all session directories for this animal
    animal_path = os.path.join(dst_path, animal)
    if not os.path.exists(animal_path):
        return []
    
    npy_files = []
    
    # Iterate through all sessions
    for session in os.listdir(animal_path):
        session_path = os.path.join(animal_path, session)
        kilosort_path = os.path.join(session_path, 'kilosort')
        rawwaveforms_path = os.path.join(kilosort_path, 'RawWaveforms')
        phy_path = os.path.join(kilosort_path, '.phy')

        if require_phy_folder and not os.path.isdir(phy_path):
            continue
        
        # Check if RawWaveforms folder exists
        if os.path.isdir(rawwaveforms_path):
            # Get all npy files in that folder
            for npy_file in Path(rawwaveforms_path).glob('*.npy'):
                npy_files.append(str(npy_file))
            
            # add spike_clusters.npy file as an input as well
            # this is a bit of a hack to trigger re-running if good/mua/noise labels change in cluster_group.tsv
            spike_clusters_file = os.path.join(session_path, 'kilosort',  'unit_match', 'spike_clusters.npy')
            npy_files.append(spike_clusters_file)
            # add cluster_group.tsv file as an input as well
            cluster_group_file = os.path.join(session_path, 'kilosort',  'unit_match', 'cluster_group.tsv')
            npy_files.append(cluster_group_file)
            
    
    return npy_files

rule extract_raw_waveforms:
    input:
        # these inputs are better defined, but for now we will have to use the others which do not
        # trigger rerun
        oebin_file=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'structure.oebin'),
        dat_file=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', '{session}.dat'),
        spk_times=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'spike_times.npy'),
        spk_clusters=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'spike_clusters.npy'),
        unit_labels=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'cluster_group.tsv')
        # oebin_file=lambda w: find_structure_oebin(os.path.join(config['src_path'], w.animal, w.session))[0],
        # dat_file=ancient(os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', '{session}.dat')),
        # spk_times=ancient(os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'spike_times.npy')),
        # spk_clusters=ancient(os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'spike_clusters.npy')),
        # unit_labels=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'cluster_group.tsv')
    output:
        ready=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'RawWaveforms', 'RawWaveforms.ready')
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


rule create_symlinks_for_unit_match:
    input:
        oebin_file=lambda w: find_structure_oebin(os.path.join(config['src_path'], w.animal, w.session))[0],
        dat_file=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', '{session}.dat'),
        spk_times=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'spike_times.npy'),
        spk_clusters=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'spike_clusters.npy'),
        cluster_group=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'cluster_group.tsv'),
        cluster_KSlabel=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'cluster_KSLabel.tsv'),
        channel_positions=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'channel_positions.npy')
    output:
        oebin_symlink=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'structure.oebin'),
        dat_symlink=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', '{session}.dat'),
        spk_times_symlink=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'spike_times.npy'),
        spk_clusters_symlink=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'spike_clusters.npy'),
        cluster_group_symlink=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'cluster_group.tsv'),
        channel_positions_symlink=os.path.join(config['dst_path'], '{animal}', '{session}', 'kilosort', 'unit_match', 'channel_positions.npy')
    run:
        # before loading data, update cluster_group.tsv
        # we have to "merge" the labels in cluster_KSlabel.tsv with the ones in cluster_group.tsv
        # cluster_group.tsv has priority, so only units that are not labeled in cluster_group.tsv will get labels from cluster_KSlabel.tsv
        ks_dir = os.path.dirname(input.cluster_group)
        cluster_group_df_updated = update_cluster_group_with_KSlabel(ks_dir)
        if cluster_group_df_updated is None:
            # this means no update was needed, so we can just create a symlink to the original cluster_group.tsv
            os.makedirs(os.path.dirname(output.cluster_group_symlink), exist_ok=True)
            if not os.path.exists(output.cluster_group_symlink):
                os.symlink(input.cluster_group, output.cluster_group_symlink)
                print(f"Created symlink for {input.cluster_group} at {output.cluster_group_symlink}.")
        else:
            # save the updated cluster_group_df to the output path
            cluster_group_symlink_dir = os.path.dirname(output.cluster_group_symlink)
            os.makedirs(cluster_group_symlink_dir, exist_ok=True)
            cluster_group_df_updated.to_csv(output.cluster_group_symlink, sep='\t', index=False)
            print(f"Updated {output.cluster_group_symlink} with labels from {input.cluster_KSlabel}.")
        # now create symlinks for the other files
        os.makedirs(os.path.dirname(output.oebin_symlink), exist_ok=True)
        if not os.path.exists(output.oebin_symlink):
            os.symlink(input.oebin_file, output.oebin_symlink)
            print(f"Created symlink for {input.oebin_file} at {output.oebin_symlink}.")
        os.makedirs(os.path.dirname(output.dat_symlink), exist_ok=True)
        if not os.path.exists(output.dat_symlink):
            os.symlink(input.dat_file, output.dat_symlink)
            print(f"Created symlink for {input.dat_file} at {output.dat_symlink}.")
        os.makedirs(os.path.dirname(output.spk_times_symlink), exist_ok=True)
        if not os.path.exists(output.spk_times_symlink):
            os.symlink(input.spk_times, output.spk_times_symlink)
            print(f"Created symlink for {input.spk_times} at {output.spk_times_symlink}.")
        os.makedirs(os.path.dirname(output.spk_clusters_symlink), exist_ok=True)
        if not os.path.exists(output.spk_clusters_symlink):
            os.symlink(input.spk_clusters, output.spk_clusters_symlink)
            print(f"Created symlink for {input.spk_clusters} at {output.spk_clusters_symlink}.")
        os.makedirs(os.path.dirname(output.channel_positions_symlink), exist_ok=True)
        if not os.path.exists(output.channel_positions_symlink):
            os.symlink(input.channel_positions, output.channel_positions_symlink)
            print(f"Created symlink for {input.channel_positions} at {output.channel_positions_symlink}.")
        # finally create a symlink for RawWaveforms folder if it exists
        rawwaveforms_folder = os.path.join(os.path.dirname(input.spk_times), 'RawWaveforms')
        rawwaveforms_symlink = os.path.join(os.path.dirname(output.spk_times_symlink), 'RawWaveforms')
        if os.path.isdir(rawwaveforms_folder):
            if not os.path.exists(rawwaveforms_symlink):
                os.symlink(rawwaveforms_folder, rawwaveforms_symlink)
                print(f"Created symlink for {rawwaveforms_folder} at {rawwaveforms_symlink}.")