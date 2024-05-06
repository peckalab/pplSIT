import h5py, os, sys
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.behavior import get_idxs_behav_state, get_idxs_neuro_state


# move to YAML?
ft = 'tSNE'
fp = 70
speed_thresh = 0.04  # m/s

s_path  = os.path.dirname(snakemake.input[0])
source  = os.path.dirname(os.path.dirname(s_path))
session = os.path.basename(s_path)

with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    events = np.array(f['processed']['sound_events'])
with h5py.File(snakemake.input[1], 'r') as f:
    speed = np.array(f['speed'])
    hd    = np.array(f['hd'])

# auditory state (BGR, SIL etc.) and speed filter
idxs_sta_ev = np.where(speed[events[:, 2].astype(np.int32)] < speed_thresh)[0]  # define speed filter here
idxs_bgr_ev  = np.where(events[:, 1] == 1)[0]
idxs_sil_ev  = np.where(events[:, 1] == 0)[0]
idxs_tgt_ev  = np.where(events[:, 1] == 2)[0]
idxs_noi_ev  = np.where(events[:, 1] == -1)[0]
idxs_tri_ev  = np.where( (events[:, 1] == 1) | (events[:, 1] == 2) )[0]  # in the trial

# behavioral state
idxs_tl_tgt_succ = []
idxs_tl_pellet   = []
for tgt_rec in tgt_mx[tgt_mx[:, 4] == 1]:
    idxs_tl_tgt_succ += list(range(tgt_rec[2], tgt_rec[3]))
    idxs_tl_pellet   += list(range(tgt_rec[3] + 1*100, tgt_rec[3] + 6*100))
idxs_tl_tgt_succ = np.array(idxs_tl_tgt_succ)
idxs_tl_pellet   = np.array(idxs_tl_pellet)

idxs_AL_tl = get_idxs_behav_state(source, session, idxs_tl_tgt_succ, fit_type=ft, fit_parm=fp, sigma=0.3, margin=10, bin_count=100)
idxs_PC_tl = get_idxs_behav_state(source, session, idxs_tl_pellet, fit_type=ft, fit_parm=fp, sigma=0.3, margin=10, bin_count=100)

idxs_AL_ev = np.array([k for k, x in enumerate(events) if x[2] in idxs_AL_tl], dtype=np.int32)
idxs_PC_ev = np.array([k for k, x in enumerate(events) if x[2] in idxs_PC_tl], dtype=np.int32)


# final separation
idxs_AL_bgr_ev  = np.intersect1d(idxs_AL_ev, idxs_bgr_ev)
idxs_AL_tgt_ev  = np.intersect1d(idxs_AL_ev, idxs_tgt_ev)
idxs_AL_sil_ev  = np.intersect1d(idxs_AL_ev, idxs_sil_ev)

idxs_PH_ev = np.array([x for x in range(len(events)) if not x in idxs_AL_ev])
idxs_PH_bgr_ev = np.intersect1d(idxs_PH_ev, idxs_bgr_ev)
idxs_PH_sil_ev = np.intersect1d(idxs_PH_ev, idxs_sil_ev)

idxs_AL_bgr_sta_ev = np.intersect1d(idxs_AL_bgr_ev, idxs_sta_ev)
idxs_PH_bgr_sta_ev = np.intersect1d(idxs_PH_bgr_ev, idxs_sta_ev)
idxs_AL_sil_sta_ev = np.intersect1d(idxs_AL_sil_ev, idxs_sta_ev)
idxs_PH_sil_sta_ev = np.intersect1d(idxs_PH_sil_ev, idxs_sta_ev)

# ------- computing AL / PH from neural state
idxs_neuro_AL_b_ev = get_idxs_neuro_state(source, session, idxs_AL_bgr_sta_ev)  # basically for BGR / TGT, when sound is ON
idxs_neuro_AL_bgr_ev = np.intersect1d(idxs_neuro_AL_b_ev, idxs_bgr_ev)  # ??
#idxs_neuro_AL_bgr_sta_ev = np.intersect1d(idxs_neuro_AL_bgr_ev, idxs_sta_ev)
idxs_neuro_AL_s_ev = get_idxs_neuro_state(source, session, idxs_AL_sil_sta_ev)  # basically for BGR / TGT, when sound is ON
idxs_neuro_AL_sil_ev = np.intersect1d(idxs_neuro_AL_s_ev, idxs_sil_ev)  # ??
#idxs_neuro_AL_sil_sta_ev = np.intersect1d(idxs_neuro_AL_sil_ev, idxs_sta_ev)

# PH would be all not AL in BGR / TGT
idxs_neuro_PH_bgr_ev = np.array([x for x in np.arange(len(events)) if not x in idxs_neuro_AL_b_ev])
idxs_neuro_PH_bgr_ev = np.intersect1d(idxs_neuro_PH_bgr_ev, idxs_bgr_ev)
#idxs_neuro_PH_bgr_sta_ev = np.intersect1d(idxs_neuro_PH_bgr_ev, idxs_sta_ev)
idxs_neuro_PH_sil_ev = np.array([x for x in np.arange(len(events)) if not x in idxs_neuro_AL_s_ev])
idxs_neuro_PH_sil_ev = np.intersect1d(idxs_neuro_PH_sil_ev, idxs_sil_ev)
#idxs_neuro_PH_sil_sta_ev = np.intersect1d(idxs_neuro_PH_sil_ev, idxs_sta_ev)


# saving states
with h5py.File(snakemake.output[0], 'w') as f:
    # behavioral states (AL - active listening, PC - pellet chasing)
    f.create_dataset('idxs_AL_tl', data=idxs_AL_tl)
    f.create_dataset('idxs_PC_tl', data=idxs_PC_tl)
    f.create_dataset('idxs_AL_ev', data=idxs_AL_ev)
    f.create_dataset('idxs_PC_ev', data=idxs_PC_ev)

    # neural states, based on behavioral states ()
    f.create_dataset('idxs_neuro_AL_bgr_ev', data=idxs_neuro_AL_bgr_ev)
    f.create_dataset('idxs_neuro_AL_sil_ev', data=idxs_neuro_AL_sil_ev)
    f.create_dataset('idxs_neuro_PH_bgr_ev', data=idxs_neuro_PH_bgr_ev)
    f.create_dataset('idxs_neuro_PH_sil_ev', data=idxs_neuro_PH_sil_ev)