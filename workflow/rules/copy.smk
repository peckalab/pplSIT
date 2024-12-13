import os
import subprocess


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
                dat_path = os.path.join(dirpath, filename)
                break

        if dat_path is None:
            raise ValueError("There should be one and only one subdirectory in the session path")

        subprocess.run(['ln', dat_path, output.dat])


# HARD-linking raw data to destination folder
rule copy_ephys:
    input:
        xml=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')),
        dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
    output:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat_ns=n_path('{animal}', '{session}', '{session}.dat'),
        dat_ks=k_path('{animal}', '{session}', '{session}.dat')
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