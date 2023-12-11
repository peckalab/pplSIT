import os
import xml.etree.ElementTree as ET


# ------ configuration ----------------------

src_path = config['src_path']
dst_path = config['dst_path']
session  = config['sessions'][0]

# wrappers to build absolute paths to source / destination folders
abs_src = lambda file_name: os.path.join(src_path, session.split('_')[0], session, file_name)
abs_dst = lambda file_name: os.path.join(dst_path, session.split('_')[0], session, file_name)

# parse original .xml file to get channel groups for spikesorting
xml_path = abs_src('%s.xml' % session)
root = ET.parse(xml_path).getroot()
ch_groups_xml = root.findall('spikeDetection')[0].findall('channelGroups')[0]

channel_groups = {}
for i, grp in enumerate(ch_groups_xml):
    for channels in grp:
        if not channels.tag == 'channels':
            continue
        channel_groups[i+1] = [int(ch.text) for ch in channels]

# inputs to klustakwik
electrodes = list(channel_groups.keys())
num_fet = ["".join(['1' for x in range(len(channel_groups[el]) * 3 + 1)]) for el in electrodes]

# ------ rules -------------------------------

# copying raw data to destination folder
rule copy_dat_xml:
    input:
        xml=abs_src('%s.xml' % session),
        dat=abs_src('%s.dat' % session)
    output:
        xml=abs_dst('%s.xml' % session),
        dat=temp(abs_dst('%s.dat' % session))
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
            session + '.xml'
        )

# extract spikes
rule extractspikes:
    input:
        xml=abs_dst('%s.xml' % session),
        fil=abs_dst('%s.fil' % session)
    output:
        [abs_dst('%s.res.%d' % (session, el)) for el in channel_groups.keys()],
        [abs_dst('%s.spk.%d' % (session, el)) for el in channel_groups.keys()],
    shell:
        "cd %s; %s %s" % (
            abs_dst(''),
            os.path.join(config['ndm_path'], "ndm_extractspikes"),
            session + '.xml'
        )

# compute PCA
rule PCA:
    input:
        [abs_dst('%s.spk.%d' % (session, el)) for el in channel_groups.keys()]
    output:
        [abs_dst('%s.fet.%d' % (session, el)) for el in channel_groups.keys()]
    shell:
        "cd %s; %s %s" % (
            abs_dst(''),
            os.path.join(config['ndm_path'], "ndm_pca"),
            session + '.xml'
        )

# sort clusters with KlustaKwik
rule kkwik:
    input:
        [abs_dst('%s.fet.%d' % (session, el)) for el in channel_groups.keys()],
        [abs_dst('%s.res.%d' % (session, el)) for el in channel_groups.keys()]
    output:
        [abs_dst('%s.clu.%d' % (session, el)) for el in channel_groups.keys()]
    shell:
        "cd %s; " % abs_dst('') + "".join(["%s %s %s %s %s; " % (
            os.path.join(config['kkwik_path'], "KlustaKwik"),
            session,
            el,
            config['kwik_args'],
            fet
        ) for el, fet in zip(electrodes, num_fet)])

# do some clean up!


