import os


# TODO: mark all raw data files as read-only

# copying raw data to destination folder
rule copy:
    input:
        xml=os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml'),
        dat=os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat')
    output:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat=n_path('{animal}', '{session}', '{session}.dat')
    shell:
        "cp {input.xml} {output.xml}; cp {input.dat} {output.dat}"