import os
import subprocess


# Creates a hard link in the session folder to the actual raw files with recorded ephys data.
# For Neuropixels recordings, also creates a hard link to the ADC raw data file.
rule move_dat_from_subfolder:
    output:
        dat=os.path.join(config['src_path'], '{animal}', '{session}', '{session}.dat')
    run:
        import os
        import shutil

        # Define the source path
        session_path = os.path.join(config['src_path'], wildcards.animal, wildcards.session)

        dat_path = None
        for dirpath, dirnames, filenames in os.walk(session_path):
            for filename in [f for f in filenames if f.endswith('.dat')]:
                parent_dirname = os.path.basename(dirpath)
                dat_path = os.path.join(dirpath, filename)

                if parent_dirname.find('OneBox-ADC') > 0:  # this is ADC dat file, special case for NP
                    subprocess.run(['ln', dat_path, os.path.join(session_path, 'ADC.dat')])

                    # assume here should be timestamps file too - need to move it up as well
                    adc_ts_path = os.path.join(dirpath, 'timestamps.npy')
                    subprocess.run(['ln', adc_ts_path, os.path.join(session_path, 'ADC_timestamps.npy')])
                else:
                    subprocess.run(['ln', dat_path, output.dat])

        if dat_path is None:
            raise ValueError("There should be at least one .dat file in the session path")


# HARD-linking raw data to destination folder
rule copy_ephys_ns:
    input:
        xml=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')),
        dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
    output:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat=n_path('{animal}', '{session}', '{session}.dat')
    shell:
        "ln {input.xml} {output.xml}; ln {input.dat} {output.dat}"

rule copy_ephys_ks:
    input:
        dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
    output:
        dat=k_path('{animal}', '{session}', '{session}.dat')
    shell:
        "ln {input.dat} {output.dat}"


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
