import os
import numpy as np
import xml.etree.ElementTree as ET


class XMLHero:
    # class to work with the Neurosuite XML file

    def __init__(self, xml_path):
        if not os.path.exists(xml_path):
            raise ValueError('Given XML file does not exist')
        self.xml_path = xml_path

    def get_channel_groups(self):
        # parse original .xml file to get channel groups for spikesorting
        root = ET.parse(self.xml_path).getroot()
        ch_groups_xml = root.findall('spikeDetection')[0].findall('channelGroups')[0]

        channel_groups = {}
        for i, grp in enumerate(ch_groups_xml):
            for channels in grp:
                if not channels.tag == 'channels':
                    continue
                channel_groups[i+1] = [int(ch.text) for ch in channels]
        return channel_groups

    def get_electrodes(self):
        channel_groups = self.get_channel_groups()
        return list(channel_groups.keys())

    def get_fet_string(self, electrode):
        channel_groups = self.get_channel_groups()
        features = ["".join(['1' for x in range(len(chs) * 3 + 1)]) for el, chs in channel_groups.items()]
        el_idx = list(channel_groups.keys()).index(int(electrode))
        return features[el_idx]

    def get_sampling_rate(self):
        root = ET.parse(self.xml_path).getroot()
        sampling_rate = root.findall('acquisitionSystem')[0].findall('samplingRate')[0]
        return int(sampling_rate.text)


def load_clu_res(where):
    """
    Neurosuite files:
    
    dat     - raw signal in binary (usually int16) format as a matrix channels x signal
    lfp     - raw signal, downsampled (historically to 1250Hz)
    fet     - list of feature vectors for every spike for a particular electrode
    spk     - list of spike waveforms for every spike for a particular electrode, binary
    res     - spike times in samples for all clusters (units) from a particular electrode
    clu     - list of cluster (unit) numbers for each spike from 'res'

    Load spike times from 'clu' (clusters) and 'res' (spike times) files generated by KlustaKwik.

    :param where:       path to the session folder
    :return:            a dict in a form like {<clustered_unit_no>: <spike_times>, ...}
    """
    filebase = os.path.basename(where)
    clu_files = [f for f in os.listdir(where) if f.find('.clu.') > 0]
    if not len(clu_files) > 0:
        return {}
    
    idxs = [int(x.split('.')[2]) for x in clu_files]  # electrode indexes
    
    all_units = {}
    for idx in idxs:
        clu_file = os.path.join(where, '.'.join([filebase, 'clu', str(idx)]))
        res_file = os.path.join(where, '.'.join([filebase, 'res', str(idx)]))

        if not os.path.isfile(clu_file) or not os.path.isfile(res_file):
            continue

        cluster_map = np.loadtxt(clu_file, dtype=np.uint16)  # uint16 for clusters
        all_spikes = np.loadtxt(res_file, dtype=np.uint64)   # uint64 for spike times

        cluster_map = cluster_map[1:]  # remove the first element - number of clusters

        result = {}
        for cluster_no in np.unique(cluster_map)[1:]:  # already sorted / remove 1st cluster - noise
            result[cluster_no] = all_spikes[cluster_map == cluster_no]
            
        all_units[idx] = result

    return all_units