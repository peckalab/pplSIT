import numpy as np


def imro2numpy(imro_path):
    """
    read the IMRO text file and convert it
    into a numpy matrix (int32) of the IMRO format:
    virtual channel ID, shank ID, bank ID, reference ID, 
    electrode ID (0-1279, sequential numbering from the tip of the probe)
    """
    with open(imro_path, 'r') as f:
        imro_text = f.read()

    idx_last_br = imro_text.rfind(')')
    imro_text = imro_text[1:idx_last_br]  # remove leading / trailing brackets
    imro_recs = imro_text.split(')(')[1:]  # ignore first element - not a channel

    imro_mx = np.zeros([len(imro_recs), 5], dtype=np.int32)
    for i, channel_text in enumerate(imro_recs):
        imro_mx[i] = np.array(channel_text.split(' '))

    return imro_mx