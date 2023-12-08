import os


src_path = config['src_path']
dst_path = config['dst_path']
session  = config['sessions'][0]

# wrappers to build absolute paths to source / destination folders
abs_src = lambda file_name: os.path.join(src_path, session.split('_')[0], session, file_name)
abs_dst = lambda file_name: os.path.join(dst_path, session.split('_')[0], session, file_name)


# copying raw data to destination folder
rule copy_dat_xml:
    input:
        xml=abs_src('%s.xml' % session),
        dat=abs_src('%s.dat' % session)
    output:
        xml=abs_dst('%s.xml' % session),
        dat=abs_dst('%s.dat' % session)
    shell:
        "cp {input.xml} {output.xml}; cp {input.dat} {output.dat}"


# high-pass filtering
rule hipass:
    input:
        xml=abs_dst('%s.xml' % session),
        dat=abs_dst('%s.dat' % session)
    output:
        abs_dst('%s.fil' % session)
    shell:
        "cd %s; %s %s" % (
            abs_dst(''),
            os.path.join(config['ndm_path'], "ndm_hipass"),
            '%.xml' % session
        )

# extract spikes
rule extractspikes:
    input:
        abs_dst('%s.fil' % session)
    output:
        # ???
    shell:
        "cd %s; %s %s" % (
            abs_dst(''),
            os.path.join(config['ndm_path'], "ndm_extractspikes"),
            '%.xml' % session
        )

# compute PCA

# execute KlustaKwik

# start from here actually..