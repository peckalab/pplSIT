import os
import xml.etree.ElementTree as ET


# parse original .xml file to get channel groups for spikesorting
def get_channel_groups(where, animal, session):
    xml_path = os.path.join(where, animal, session, '%s.xml' % session)
    root = ET.parse(xml_path).getroot()
    ch_groups_xml = root.findall('spikeDetection')[0].findall('channelGroups')[0]

    channel_groups = {}
    for i, grp in enumerate(ch_groups_xml):
        for channels in grp:
            if not channels.tag == 'channels':
                continue
            channel_groups[i+1] = [int(ch.text) for ch in channels]
    return channel_groups

def get_electrodes(where, animal, session):
    channel_groups = get_channel_groups(where, animal, session)
    return list(channel_groups.keys())

def get_features(where, animal, session):
    channel_groups = get_channel_groups(where, animal, session)
    electrodes = list(channel_groups.keys())
    return ["".join(['1' for x in range(len(channel_groups[el]) * 3 + 1)]) for el in electrodes]

def get_fet_string(where, animal, session, electrode):
    electrodes = get_electrodes(where, animal, session)
    el_idx = electrodes.index(int(electrode))
    return get_features(where, animal, session)[el_idx]