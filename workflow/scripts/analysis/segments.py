import h5py, os, sys, json
import numpy as np
import itertools

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.states import get_state_as_periods
from utils.behavior import get_idxs_as_periods


# read datasets
meta_file = snakemake.input[0]
smk_cfg   = snakemake.config['segments']

s_path  = os.path.dirname(meta_file)
session = os.path.basename(s_path)
animal  = session.split('_')[0]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    if 'distractor_matrix' in f['processed']:
        dis_mx = np.array(f['processed']['distractor_matrix'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])

# indices in event space
x_pos_ev = tl[sound_events[:, 2].astype(np.int32)][:, 1]
y_pos_ev = tl[sound_events[:, 2].astype(np.int32)][:, 2]
# speed of an event is the maximum speed during that sound event
# from start idx sound_events[:, 2].astype(np.int32) to the start of the next event (or end of timeline)
speed_ev = np.zeros(len(sound_events))
for i in range(len(sound_events)):
    idx_start = int(sound_events[i][2])

    if i < len(sound_events) - 1:
        idx_end = int(sound_events[i + 1][2])
    else:
        idx_end = len(tl)

    # guard against zero-length or reversed slices
    if idx_end <= idx_start:
        # fallback: use single-sample speed at idx_start if possible
        if 0 <= idx_start < len(tl):
            speed_ev[i] = tl[idx_start, 3]
        else:
            speed_ev[i] = 0.0
        continue

    speed_ev[i] = tl[idx_start:idx_end][:, 3].max()

speed_max = 0.04
idxs_sta_ev = np.where(speed_ev < speed_max)[0]
idxs_run_ev = np.where(speed_ev > speed_max)[0]
idxs_bgr_ev = np.where(sound_events[:, 1] == 1)[0]
idxs_sil_ev = np.where(sound_events[:, 1] == 0)[0]
idxs_tgt_ev = np.where(sound_events[:, 1] == 2)[0]
idxs_di1_ev = np.where(sound_events[:, 1] == 3)[0]
idxs_di2_ev = np.where(sound_events[:, 1] == 4)[0]


# times of first tgt success pulses
tgt_mx_succ = tgt_mx[tgt_mx[:, 4] == 1]
idxs_tgt_first_ev = tgt_mx_succ[:, 0]
tgt_first_t = sound_events[idxs_tgt_first_ev][:, 0]

# success stays
idxs_tgt_sta_succ = []
for tgt_rec in tgt_mx_succ:
    idxs_tgt_sta_succ += list(np.arange(tgt_rec[0], tgt_rec[1] + 1))
idxs_tgt_sta_succ = np.array(idxs_tgt_sta_succ, dtype=np.int64)

tgt_sta_succ_mx = np.zeros([len(tgt_mx_succ), 4])
for i, tgt_rec in enumerate(tgt_mx_succ):
    x_pos = tl[np.arange(tgt_rec[2], tgt_rec[3])][:, 1]
    y_pos = tl[np.arange(tgt_rec[2], tgt_rec[3])][:, 2]
    
    tgt_sta_succ_mx[i] = [tgt_rec[0], tgt_rec[1], x_pos.mean(), y_pos.mean()]

# distractor fails
if cfg['experiment']['distractor_fail']:
    dis_mx_fail = dis_mx[dis_mx[:, 4] == 1]
    idxs_dis_first_ev = dis_mx_fail[:, 0]
    dis_first_t = sound_events[idxs_dis_first_ev][:, 0]
    idxs_dis_fail = []
    for dis_rec in dis_mx_fail:
        idxs_dis_fail += list(np.arange(dis_rec[0], dis_rec[1] + 1))
    idxs_dis_fail = np.array(idxs_dis_fail, dtype=np.int64)
    dis_fail_mx = np.zeros([len(dis_mx_fail), 4])
    for i, dis_rec in enumerate(dis_mx_fail):
        x_pos = tl[np.arange(dis_rec[2], dis_rec[3])][:, 1]
        y_pos = tl[np.arange(dis_rec[2], dis_rec[3])][:, 2]
        
        dis_fail_mx[i] = [dis_rec[0], dis_rec[1], x_pos.mean(), y_pos.mean()]
    


# stationary states, including TGT (alternative to success stays above)
tgt_sta_mx, idxs_tgt_sta = get_state_as_periods(s_path, 'TGT', None,  None, smk_cfg['tgt_sta_min_pulses'], strip_l=smk_cfg['strip_l']['tgt'], strip_r=smk_cfg['strip_r']['tgt'])
bgr_sta_mx, idxs_bgr_sta = get_state_as_periods(s_path, 'BGR', 'STA', None, smk_cfg['bgr_sta_min_pulses'], strip_l=smk_cfg['strip_l']['bgr'], strip_r=smk_cfg['strip_r']['bgr'])
sil_sta_mx, idxs_sil_sta = get_state_as_periods(s_path, 'SIL', 'STA', None, smk_cfg['sil_sta_min_pulses'], strip_l=smk_cfg['strip_l']['sil'], strip_r=smk_cfg['strip_r']['sil'])

# running states
bgr_run_mx, idxs_bgr_run = get_state_as_periods(s_path, 'BGR', 'RUN', None, smk_cfg['bgr_run_min_pulses'], strip_l=smk_cfg['strip_l']['bgr'], strip_r=smk_cfg['strip_r']['bgr'])
sil_run_mx, idxs_sil_run = get_state_as_periods(s_path, 'SIL', 'RUN', None, smk_cfg['sil_run_min_pulses'], strip_l=smk_cfg['strip_l']['sil'], strip_r=smk_cfg['strip_r']['sil'])

# distractors
di1_sta_mx, di2_sta_mx = None, None
distr_count = int(cfg['experiment']['distractor_islands'])
if cfg['sound']['sounds']['distractor1']['enabled'] and distr_count > 0:
    di1_sta_mx, idxs_di1_sta = get_state_as_periods(s_path, 'DI1', 'STA', None, smk_cfg['dis_sta_min_pulses'], strip_l=smk_cfg['strip_l']['dis'], strip_r=smk_cfg['strip_r']['dis'])
if cfg['sound']['sounds']['distractor2']['enabled'] and distr_count > 1:
    di2_sta_mx, idxs_di2_sta = get_state_as_periods(s_path, 'DI2', 'STA', None, smk_cfg['dis_sta_min_pulses'], strip_l=smk_cfg['strip_l']['dis'], strip_r=smk_cfg['strip_r']['dis'])

# target visits (no matter run or stationary, important is where)
r_max = smk_cfg['visits']['radius']  # in meters
visits_mxs  = {'SIL': [], 'BGR': []}  # visits as ranges of indices to sound events and positions
visits_idxs = {'SIL': [], 'BGR': []}  # visits as ranges of indices to sound events
titles = ['SIL', 'BGR']

for i, tgt_xy in enumerate(tgt_sta_succ_mx):
    idxs_around_ev = np.where( (x_pos_ev - tgt_xy[2])**2 + (y_pos_ev - tgt_xy[3])**2 < r_max**2 )[0]
    idxs_around_sil_ev = np.intersect1d(idxs_around_ev, idxs_sil_ev)
    idxs_around_bgr_ev = np.intersect1d(idxs_around_ev, idxs_bgr_ev)

    for j, idxs_events in enumerate([idxs_around_sil_ev, idxs_around_bgr_ev]):
        idxs_coll = list(visits_mxs[titles[j]])
        idxs_sil_as_per = get_idxs_as_periods(idxs_events)
        if len(idxs_sil_as_per.shape) == 1:
            idxs_sil_as_per = np.array([idxs_sil_as_per], dtype=np.int64)
        idxs_sil_as_per = idxs_sil_as_per[np.where(np.diff(idxs_sil_as_per, axis=1) > smk_cfg['visits']['min_pulses'] - 2)[0]]
        for per in idxs_sil_as_per:
            idx_tl_s = int(sound_events[per[0]][2])
            idx_tl_e = int(sound_events[per[1]][2])
            x_pos = tl[np.arange(idx_tl_s, idx_tl_e)][:, 1]
            y_pos = tl[np.arange(idx_tl_s, idx_tl_e)][:, 2]
            idxs_coll.append([int(per[0]), int(per[1]), x_pos.mean(), y_pos.mean(), i])

        visits_mxs[titles[j]] = np.array(idxs_coll, dtype=np.float64)

        idxs_flat = []
        for rec in idxs_coll:
            idxs_flat += list(np.arange(int(rec[0]), int(rec[1])))
        visits_idxs[titles[j]] = np.array(idxs_flat, dtype=np.int64)

results = {
    'tgt_sta_succ_mx': tgt_sta_succ_mx,
    'idxs_tgt_sta_succ': idxs_tgt_sta_succ,

    'tgt_sta_mx': tgt_sta_mx,
    'idxs_tgt_sta': idxs_tgt_sta,
    'bgr_sta_mx': bgr_sta_mx,
    'idxs_bgr_sta': idxs_bgr_sta,
    'sil_sta_mx': sil_sta_mx,
    'idxs_sil_sta': idxs_sil_sta,

    'bgr_run_mx': bgr_run_mx,
    'idxs_bgr_run': idxs_bgr_run,
    'sil_run_mx': sil_run_mx,
    'idxs_sil_run': idxs_sil_run,

    'bgr_vis_mx': visits_mxs['BGR'],
    'idxs_bgr_vis': visits_idxs['BGR'],
    'sil_vis_mx': visits_mxs['SIL'],
    'idxs_sil_vis': visits_idxs['SIL'],
}

if di1_sta_mx is not None:
    results['di1_sta_mx'] = di1_sta_mx
    results['idxs_di1_sta'] = idxs_di1_sta
if di2_sta_mx is not None:
    results['di2_sta_mx'] = di2_sta_mx
    results['idxs_di2_sta'] = idxs_di2_sta
if cfg['experiment']['distractor_fail']:
    results['dis_fail_mx'] = dis_fail_mx
    results['idxs_dis_fail'] = idxs_dis_fail

# AL / PH states
if smk_cfg['ensembles']:
    bgr_sta_AL_mx, idxs_bgr_sta_AL = get_state_as_periods(s_path, 'BGR', 'STA', 'AL', 12)
    bgr_sta_PH_mx, idxs_bgr_sta_PH = get_state_as_periods(s_path, 'BGR', 'STA', 'PH', 2, strip_l=1, strip_r=0)
    sil_sta_AL_mx, idxs_sil_sta_AL = get_state_as_periods(s_path, 'SIL', 'STA', 'AL', 4)
    sil_sta_PH_mx, idxs_sil_sta_PH = get_state_as_periods(s_path, 'SIL', 'STA', 'PH', 2, strip_l=1, strip_r=0)

    results['bgr_sta_AL_mx'] = bgr_sta_AL_mx
    results['idxs_bgr_sta_AL'] = idxs_bgr_sta_AL
    results['bgr_sta_PH_mx'] = bgr_sta_PH_mx
    results['idxs_bgr_sta_PH'] = idxs_bgr_sta_PH
    results['sil_sta_AL_mx'] = sil_sta_AL_mx
    results['idxs_sil_sta_AL'] = idxs_sil_sta_AL
    results['sil_sta_PH_mx'] = sil_sta_PH_mx
    results['idxs_sil_sta_PH'] = idxs_sil_sta_PH

# dump to H5
with h5py.File(snakemake.output[0], 'w') as out_file:
    for name, segment_element in results.items():
        # Explicitly save index arrays as int64 to prevent HDF5 dtype inference issues
        if isinstance(segment_element, np.ndarray):
            if segment_element.dtype == np.object_:
                # Object arrays must be converted - use float64 for _mx arrays, int64 for idxs_
                if name.startswith('idxs_'):
                    segment_element = segment_element.astype(np.int64)
                else:
                    segment_element = segment_element.astype(np.float64)
            elif name.startswith('idxs_'):
                segment_element = segment_element.astype(np.int64)
        out_file.create_dataset(name, data=segment_element)
