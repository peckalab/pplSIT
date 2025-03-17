import numpy as np
import os, json
import xml.etree.ElementTree as ET


def oe2kilosort(oe_xml_path):
    # extract probe information from the OpenEphys XML file and 
    # store in kilosort probe format dictionnary.
    # Works for Neuropixels OneBox ONLY for now.
    root = ET.parse(oe_xml_path).getroot()

    try:
        onebox = [x for x in root.findall('SIGNALCHAIN')[0].findall('PROCESSOR') if dict(x.items())['name'] == 'OneBox'][0]
    except IndexError:
        raise ValueError('OneBox processor not found in the OpenEphys signal chain')

    # get sampling rate from probe settings - just in case
    probe_stream = [x for x in onebox.findall('STREAM') if dict(x.items())['name'] == 'ProbeA'][0]
    sample_rate = float(probe_stream.attrib['sample_rate'])

    probe = onebox.findall('EDITOR')[0].findall('NP_PROBE')[0]
    channels = probe.findall('CHANNELS')[0]
    x_pos = probe.findall('ELECTRODE_XPOS')[0]
    y_pos = probe.findall('ELECTRODE_YPOS')[0]

    xc = [float(x) for x in x_pos.attrib.values()]
    yc = [float(y) for y in y_pos.attrib.values()]
    kcoords = [int(desc.split(':')[1]) for ch, desc in channels.attrib.items()]
    ch_count = len(channels.attrib)

    # get channel map
    ch_map = [x for x in root.findall('SIGNALCHAIN')[0].findall('PROCESSOR') if dict(x.items())['name'] == 'Channel Map'][0]
    channels_xml = list(list(ch_map.findall('CUSTOM_PARAMETERS')[0])[0])  # taking the first "stream"

    # compute mapping indices
    channel_map_list = [int(ch.attrib['index']) for ch in channels_xml]
    if 384 in channel_map_list: channel_map_list.remove(384)

    channels_probe = [int(ch[2:]) for ch in channels.attrib.keys()]
    idxs_channels = [channels_probe.index(ch) for ch in channel_map_list]

    return {
        'chanMap': [x for x in range(ch_count)],
        'xc': [xc[i] for i in idxs_channels],
        'yc': [yc[i] for i in idxs_channels],
        'kcoords': [kcoords[i] for i in idxs_channels],
        'n_chan': ch_count
    }


def kilosort_probe_manually():
    # requires conda kilosort env
    from kilosort.io import save_probe

    # for the dual A1 / PPC implant with a single linear probe in A1 (32 channels)
    # and a 4-shank probe in PPC (another 32 channels)

    cn_count = 64  # channels to sort (no ADC/DAC channels)
    channels = np.arange(cn_count)
    kcoords  = np.zeros(cn_count)  # shank numbers
    for i in range(4):
        kcoords[32 + i*8:32 + (i+1)*8] = i+1

    # probes configurations (site coordinates)
    xc, yc = np.zeros(cn_count), np.zeros(cn_count)

    # A1 probe is 1 shank, -5mm ML, Y-spacing is 35um
    xc[0:32] = -5000
    yc[0:32] = np.array([35*i for i in range(32)])

    # PPC probe is 4 shank, 3mm ML, V-shape
    for i in range(4):
        xc[32 + i*8:32 + (i+1)*8] = -4000 + np.array([i*200 + 3*j for j in range(8)])
        yc[32 + i*8:32 + (i+1)*8] = np.abs(np.array([20*j for j in range(8)]) - 20*4)

    probe = {
        'chanMap': channels,
        'xc': xc,
        'yc': yc,
        'kcoords': kcoords,
        'n_chan': cn_count
    }

    save_probe(probe, 'probe.json')
