import os


# TODO: mark all raw data files as read-only

# copying raw data to destination folder
rule copy_ephys:
    input:
        xml=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')),
        dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
    output:
        xml=protected(n_path('{animal}', '{session}', '{session}.xml')),
        dat=protected(n_path('{animal}', '{session}', '{session}.dat'))
    shell:
        "cp {input.xml} {output.xml}; cp {input.dat} {output.dat}"