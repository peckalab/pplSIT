import os
import xml.etree.ElementTree as ET


# parse original .xml file to get channel groups for spikesorting
def get_channel_groups(animal, session):
    xml_path = os.path.join(config['src_path'], animal, session, '%s.xml' % session)
    root = ET.parse(xml_path).getroot()
    ch_groups_xml = root.findall('spikeDetection')[0].findall('channelGroups')[0]

    channel_groups = {}
    for i, grp in enumerate(ch_groups_xml):
        for channels in grp:
            if not channels.tag == 'channels':
                continue
            channel_groups[i+1] = [int(ch.text) for ch in channels]
    return channel_groups

def get_electrodes(animal, session):
    channel_groups = get_channel_groups(animal, session)
    return list(channel_groups.keys())

def get_features(animal, session):
    channel_groups = get_channel_groups(animal, session)
    electrodes = list(channel_groups.keys())
    return ["".join(['1' for x in range(len(channel_groups[el]) * 3 + 1)]) for el in electrodes]

def get_fet_string(animal, session, electrode):
    electrodes = get_electrodes(animal, session)
    el_idx = electrodes.index(int(electrode))
    return get_features(animal, session)[el_idx]


# main function to create all inputs and wildcards
def get_clu_inputs(w):
    electrodes = get_electrodes(w.animal, w.session)
    return expand(n_path(w.animal, w.session, w.session + '.clu_copy.' + '{electrode}'), electrode=electrodes)


# -------- rules ----------------------


# high-pass filtering
rule hipass:
    input:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat=n_path('{animal}', '{session}', '{session}.dat')
    output:
        temp(n_path('{animal}', '{session}', '{session}.fil'))
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "cd %s; %s %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['ndm_path'], "ndm_hipass"),
            '{params.session}.xml'
        )

# extract spikes
rule extractspikes:
    input:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        fil=n_path('{animal}', '{session}', '{session}.fil')
    output:
        n_path('{animal}', '{session}', '{session}.spk.{electrode}'),
        n_path('{animal}', '{session}', '{session}.res.{electrode}')
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "cd %s; %s %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['ndm_path'], "ndm_extractspikes"),
            '{params.session}.xml'
        )

# compute PCA
rule PCA:
    input:
        n_path('{animal}', '{session}', '{session}.spk.{electrode}')
    output:
        n_path('{animal}', '{session}', '{session}.fet.{electrode}')
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "cd %s; %s %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['ndm_path'], "ndm_pca"),
            '{params.session}.xml'
        )


# sort clusters with KlustaKwik
rule kkwik:
    input:
        n_path('{animal}', '{session}', '{session}.fet.{electrode}'),
        n_path('{animal}', '{session}', '{session}.res.{electrode}')
    output:
        n_path('{animal}', '{session}', '{session}.clu.{electrode}')
    params:
        session="{session}",
        animal="{animal}",
        electrode="{electrode}",
        fet_string=lambda w, output: get_fet_string(w.animal, w.session, w.electrode)
    shell:
        "cd %s; %s %s %s %s %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['kkwik_path'], "KlustaKwik"),
            '{params.session}',
            '{params.electrode}',
            config['kwik_args'],
            '{params.fet_string}'
        )


# backup .clu files
rule clu_backup:
    input:
        n_path('{animal}', '{session}', '{session}.clu.{electrode}')
    output:
        n_path('{animal}', '{session}', '{session}.clu_copy.{electrode}')
    shell:
        "cp {input} {output}"


 # dump units to HDF5
rule dump_units:
    input:
        get_clu_inputs
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'units.txt')
    shell:
        "touch {output}"


# TODO export spikes to units.h5