import os


rule move_dat_from_subfolder:
    output:
        dat=os.path.join(config['src_path'], '{animal}', '{session}', '{session}.dat')
    run:
        import os
        import shutil

        # Define the source path
        session_path = os.path.join(config['src_path'], wildcards.animal, wildcards.session)

        # Get a list of all items in the source directory
        items = os.listdir(session_path)

        # Filter this list to only include directories
        dirs = [item for item in items if os.path.isdir(os.path.join(session_path, item))]

        # Check if there is only one directory
        if len(dirs) != 1:
            raise ValueError("There should be one and only one subdirectory in the session path")
        
        recording_root_directory = os.path.join(session_path, dirs[0])

        # Check whether there are multiple recordings
        if len(os.listdir(os.path.join(recording_root_directory,'Record Node 117','experiment1'))) > 1:
            raise ValueError("There should be only one recording in the source path")

        # Define the source and destination paths for the .dat file
        src_dat = os.path.join(recording_root_directory,'Record Node 117','experiment1', 'recording1', 'continuous', 'Rhythm_FPGA-114.0', 'continuous.dat')
        dest_dat = output.dat

        # Copy the .dat file
        shutil.copy(src_dat, dest_dat)

        # TODO: rewrite by searching .dat in the hierarchy + create hard link!


# HARD-linking raw data to destination folder
rule copy_ephys:
    input:
        xml=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')),
        dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
    output:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat_ns=protected(n_path('{animal}', '{session}', '{session}.dat')),
        dat_ks=protected(k_path('{animal}', '{session}', '{session}.dat'))
    shell:
        "ln {input.xml} {output.xml}; ln {input.dat} {output.dat_ns}; ln {input.dat} {output.dat_ks}"


rule create_xml_from_template:
    input:
        template=ancient(config['template_xml'])
    output:
        xml=os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')
    shell:
        "cp {input.template} {output.xml}"


rule create_manual_json_from_template:
    input:
        template=ancient(config['template_manual_json'])
    output:
        man_json=os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json')
    shell:
        "cp {input.template} {output.man_json}"


rule create_kilosort_settings_from_template:
    input:
        template=ancient(config['kilosort']['settings_path'])
    output:
        kilo=os.path.join(config['src_path'], '{animal}', '{session}', 'kilosort.json')
    shell:
        "cp {input.template} {output.kilo}"


rule create_probe_from_template:
    input:
        template=ancient(config['kilosort']['probe_path'])
    output:
        kilo=os.path.join(config['src_path'], '{animal}', '{session}', 'probe.json')
    shell:
        "cp {input.template} {output.kilo}"