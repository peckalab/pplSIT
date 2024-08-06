import os
import json
import h5py
from datetime import datetime


# North, south-east, south-west - CW
channel_positions_1 = {  # before 03.2023 inclusive
    1: (0., 0.5),
    2: (0.37, -0.37),
    3: (-0.37, -0.37)
}
channel_positions_2 = {  # after 04.2023 till 05.2024 inclusive
    3: (0., 0.5),
    4: (0.37, -0.37),
    5: (-0.37, -0.37)
}
channel_positions_3 = {  # after 04.2023 till 05.2024 inclusive (TODO insert correct positions)
    3: (0., 0.5),
    4: (0.37, -0.37),
    5: (-0.37, -0.37)
}

# setup change dates
date1 = datetime.strptime('2023-03-31_23-59-00', "%Y-%m-%d_%H-%M-%S")
date2 = datetime.strptime('2024-05-31_23-59-00', "%Y-%m-%d_%H-%M-%S")


def get_speaker_positions(source, session):
    animal = session.split('_')[0]

    s_path     = os.path.join(source, animal, session)
    meta_file  = os.path.join(source, animal, session, 'meta.h5')
    json_file  = os.path.join(source, animal, session, session + '.json')

    if os.path.exists(json_file):
        with open(json_file, 'r') as f:  # read from raw JSON
            cfg = json.loads(f.read())

    elif os.path.exists(meta_file):
        with h5py.File(meta_file, 'r') as f:  # read from meta
            cfg = json.loads(f['processed'].attrs['parameters'])

    else:
        raise IOError('Config file not found')

    session_datetime = datetime.strptime(session[-19:], "%Y-%m-%d_%H-%M-%S")

    if session_datetime < date1:
        ch_pos = channel_positions_1
    elif session_datetime < date2:
        ch_pos = channel_positions_2
    else:
        print('TODO: update channel / speaker positions!')
        ch_pos = channel_positions_3
        
    chs_tgt = cfg['sound']['sounds']['target']['channels']
    chs_bgr = cfg['sound']['sounds']['background']['channels']

    try:
        ch_tgt = [ch for ch in chs_tgt if ch in ch_pos][0]
        ch_bgr = [ch for ch in chs_bgr if ch in ch_pos][0]
    except IndexError:
        raise ValueError('Channel not found in the specification')
    
    return {'BGR': ch_pos[ch_bgr], 'TGT': ch_pos[ch_tgt]}