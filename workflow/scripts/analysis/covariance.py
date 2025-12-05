import h5py, os, sys, json
import numpy as np
import itertools
from scipy import stats, signal
from sklearn import decomposition

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

# from utils.population import unit_activity_matrix
# from utils.neurosuite import get_unit_names_sorted
# from utils.psth import staple_spike_times
# from utils.spiketrain import get_shuffled
from utils.states import get_state_as_periods
from utils.correlation import cluster_corr
from utils.maths import pval2text
from scipy.spatial.distance import cosine
from sklearn.manifold import TSNE

# read datasets
s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
animal  = session.split('_')[0]

meta_file = snakemake.input[0]
unit_file = snakemake.input[1]
actm_file = snakemake.input[2]
segm_file = snakemake.input[3]

spike_times, unit_sic, unit_sic_nostim = {}, {}, {}
with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])
with h5py.File(unit_file, 'r') as f:
    unit_names = ([name for name in f])
    for unit_name in unit_names:
        spike_times[unit_name]  = np.array(f[unit_name]['spike_times'])
        unit_sic[unit_name] = float(np.array(f[unit_name]['spatial_information']))
        unit_sic_nostim[unit_name] = float(np.array(f[unit_name]['spatial_information_nostim']))

# mean firing rates
FR_mx = np.zeros([len(unit_names), 4])  # mean rate, median ISI, SIC
for i, unit_id in enumerate(unit_names):
    spiketrain = spike_times[unit_id]
    
    mean_rate = len(spiketrain) / (tl[-1][0] - tl[0][0])
    
    isis = np.diff(spiketrain)
    robust_rate = 1 / np.median(isis)

    FR_mx[i] = np.array([mean_rate, robust_rate, unit_sic[unit_id], unit_sic_nostim[unit_id]])


# getting experimental states
idxs_states = ['idxs_tgt_sta', 'idxs_bgr_sta', 'idxs_sil_sta', 'idxs_bgr_run', 'idxs_sil_run', \
               'idxs_bgr_sta_AL', 'idxs_bgr_sta_PH', 'idxs_sil_sta_AL', 'idxs_sil_sta_PH']

with h5py.File(segm_file, 'r') as f:
    event_idxs = {}
    for idxs_name in idxs_states:
        if idxs_name in f:
            event_idxs[idxs_name] = np.array(f[idxs_name])

# # positions and speed
# x_pos_ev = tl[sound_events[:, 2].astype(np.int32)][:, 1]
# y_pos_ev = tl[sound_events[:, 2].astype(np.int32)][:, 2]
# speed_ev = tl[sound_events[:, 2].astype(np.int32)][:, 3]

# # experimental states
# speed_max = 0.04
# idxs_sta_ev = np.where(speed_ev < speed_max)[0]
# idxs_run_ev = np.where(speed_ev > speed_max)[0]
# idxs_bgr_ev = np.where(sound_events[:, 1] == 1)[0]
# idxs_sil_ev = np.where(sound_events[:, 1] == 0)[0]
# idxs_tgt_ev = np.where(sound_events[:, 1] == 2)[0]
# idxs_di1_ev = np.where(sound_events[:, 1] == 3)[0]
# idxs_di2_ev = np.where(sound_events[:, 1] == 4)[0]

# # times of first tgt success pulses
# tgt_mx_succ = tgt_mx[tgt_mx[:, 4] == 1]
# idxs_tgt_first_ev = tgt_mx_succ[:, 0]
# tgt_first_t = sound_events[idxs_tgt_first_ev][:, 0]

# # success stays
# idxs_tgt_succ = []
# for tgt_rec in tgt_mx_succ:
#     idxs_tgt_succ += list(np.arange(tgt_rec[0], tgt_rec[1] + 1))
# idxs_tgt_succ = np.array(idxs_tgt_succ)

# tgt_succ_pos_mx = np.zeros([len(tgt_mx_succ), 4])
# for i, tgt_rec in enumerate(tgt_mx_succ):
#     x_pos = tl[np.arange(tgt_rec[2], tgt_rec[3])][:, 1]
#     y_pos = tl[np.arange(tgt_rec[2], tgt_rec[3])][:, 2]
    
#     tgt_succ_pos_mx[i] = [tgt_rec[0], tgt_rec[1], x_pos.mean(), y_pos.mean()]
    
# bgr_sta_long_mx, idxs_bgr_sta_long = get_state_as_periods(s_path, 'BGR', 'STA', None, 12, strip_l=1)
# sil_sta_long_mx, idxs_sil_sta_long = get_state_as_periods(s_path, 'SIL', 'STA', None, 12, strip_l=1)
# bgr_sta_AL_mx, idxs_bgr_sta_AL = get_state_as_periods(s_path, 'BGR', 'STA', 'AL', 12)
# bgr_sta_PH_mx, idxs_bgr_sta_PH = get_state_as_periods(s_path, 'BGR', 'STA', 'PH', 2, strip_l=1, strip_r=0)
# sil_sta_AL_mx, idxs_sil_sta_AL = get_state_as_periods(s_path, 'SIL', 'STA', 'AL', 4)
# sil_sta_PH_mx, idxs_sil_sta_PH = get_state_as_periods(s_path, 'SIL', 'STA', 'PH', 2, strip_l=1, strip_r=0)

# idxs_bgr_run_ev = np.intersect1d(idxs_run_ev, idxs_bgr_ev)
# idxs_sil_run_ev = np.intersect1d(idxs_run_ev, idxs_sil_ev)
# bgr_run_long_mx, idxs_bgr_run_long = get_state_as_periods(s_path, 'BGR', 'RUN', None, 2, strip_l=1)
# sil_run_long_mx, idxs_sil_run_long = get_state_as_periods(s_path, 'SIL', 'RUN', None, 4, strip_l=1)

# # distractors
# di1_sta_mx, idxs_di1_sta = get_state_as_periods(s_path, 'DI1', 'STA', None, 4, strip_l=1)
# di2_sta_mx, idxs_di2_sta = get_state_as_periods(s_path, 'DI2', 'STA', None, 4, strip_l=1)


n_clusters = 8
# combs = [
#     ['mx_50ms', 'mx_z'],
#     ['mx_50ms', 'mx_sm10_z'],
#     ['mx_50ms', 'mx_sm20_z'],
#     ['mx_50ms', 'mx_sm50_z'],
#     ['mx_50ms', 'mx_sm100_z'],
# ]
combs = [  # do z-scoring within conditions
    ['mx_50ms', 'mx'],
    ['mx_50ms', 'mx_sm10'],
    ['mx_50ms', 'mx_sm20'],
    ['mx_50ms', 'mx_sm50'],
    ['mx_50ms', 'mx_sm100'],
]
suptitles = [
    '50ms binning, no smoothing',
    '50ms binning, 0.5s gauss. sm.',
    '50ms binning, 1.0s gauss. sm.',
    '50ms binning, 2.5s gauss. sm.',
    '50ms binning, 5.0s gauss. sm.',
]
conds = {
    'all': np.arange(len(sound_events) - 1),
    'tgt_succ':     event_idxs['idxs_tgt_sta'],
    'bgr_sta': event_idxs['idxs_bgr_sta'],
    'sil_sta': event_idxs['idxs_sil_sta'],
    'bgr_run':      event_idxs['idxs_bgr_run'],
    'sil_run':      event_idxs['idxs_sil_run'],
    #'bgr_sta_AL':   event_idxs['idxs_bgr_sta_AL'],
    #'bgr_sta_PH':   event_idxs['idxs_bgr_sta_PH'],
    #'sil_sta_AL':   event_idxs['idxs_sil_sta_AL'],
    #'sil_sta_PH':   event_idxs['idxs_sil_sta_PH'],
}
vmin, vmax = -1, 1
colors = []
for i, (k, v) in enumerate(mcolors.TABLEAU_COLORS.items()):
    if i == n_clusters:
        break
    colors.append(k)




# computing correlations
with h5py.File(snakemake.output[0], 'w') as cov_file:

    for group, name in combs:
        with h5py.File(actm_file, 'r') as f:
            #mx_z = np.array(f[group][name])[:, ::5]  # !! take every 5th element to match the event resolution
            mx = np.array(f[group][name])[:, ::5]

        grp = cov_file.create_group(group + '_' + name)

        for cond_id, idxs_cond in conds.items():

            # z-scoring across conditions
            #cov_mx = np.cov(mx_z.T[idxs_cond].T)

            # z-scoring within conditions
            mx_z = stats.zscore(mx.T[idxs_cond].T, axis=1)
            mx_z_clean = np.nan_to_num(mx_z, nan=0.0)
            cov_mx = np.cov(mx_z_clean)

            cov_mx_srt, Z, labels, idxs_sort = cluster_corr(cov_mx, n_clusters=8)
            labels_dict = dict([(u, l) for u, l in zip(unit_names, labels)])

            grp_cond = grp.create_group(cond_id)
            grp_cond.create_dataset('cov_mx', data=cov_mx)
            grp_cond.create_dataset('cov_mx_srt', data=cov_mx_srt)
            grp_cond.create_dataset('linkage', data=Z)
            grp_cond.create_dataset('labels', data=labels)
            grp_cond.create_dataset('idxs_sort', data=idxs_sort)


# plotting correlation matrices
with PdfPages(snakemake.output[1]) as pdf:
    for k, (group, name) in enumerate(combs):

        # --------------- loading ------------------

        # load unit activity matrix
        with h5py.File(actm_file, 'r') as f:
            mx_z = np.array(f[group][name + '_z'])[:, ::5]  # !! take every 5th element to match the event resolution
            mx = np.array(f[group][name])[:, ::5]

        # take labels / sorting from the 'all' condition
        with h5py.File(snakemake.output[0], 'r') as f:  
            labels    = np.array(f[group + '_' + name]['all']['labels'])
            idxs_sort = np.array(f[group + '_' + name]['all']['idxs_sort'])
            
        # --------------- computing ------------------

        # compute things first
        eignvectors, exp_ratios, cov_mxs = {}, {}, {}
        for cond_id, cond_idxs in conds.items():
            if cond_id == 'all':  # skip 'all' state
                continue

            # PCA on that state - reading and plotting in the same place, not super nice tbh
            #unit_mx_state = mx_z.T[cond_idxs]

            mx_z_cond  = stats.zscore(mx.T[cond_idxs].T, axis=1)
            mx_z_clean = np.nan_to_num(mx_z_cond, nan=0.0)
            unit_mx_state = mx_z_clean.T

            ev_pca = decomposition.PCA(n_components=5)
            ev_X   = ev_pca.fit_transform(unit_mx_state)
            eignvectors[cond_id] = ev_pca.components_[0]
            exp_ratios[cond_id]  = ev_pca.explained_variance_ratio_

            # sorted covariance matrix
            with h5py.File(snakemake.output[0], 'r') as f:  
                cov_mx = np.array(f[group + '_' + name][cond_id]['cov_mx'])
                #labels = np.array(f[group + '_' + name][cond_id]['labels'])
            cov_mx_srt = cov_mx[np.ix_(idxs_sort, idxs_sort)]
            cov_mxs[cond_id] = cov_mx_srt

        # eigenvector correlations / cosince similarity
        corr_mx = np.zeros([len(eignvectors), len(eignvectors)])
        coss_mx = np.zeros([len(eignvectors), len(eignvectors)])

        for i, (cond_id_1, v1) in enumerate(eignvectors.items()):
            for j, (cond_id_2, v2) in enumerate(eignvectors.items()):
                # Cosine similarity (1 - cosine distance)
                cos_sim = 1 - cosine(np.abs(v1), np.abs(v2))
        
                # standard pearsons corr coeff
                #corr, pv = scipystats.pearsonr(v1, v2)
                corr, pv = stats.pearsonr(np.abs(v1), np.abs(v2))

                corr_mx[i][j] = corr
                coss_mx[i][j] = cos_sim

        # --------------- plotting ------------------

        # --------------- sorted heatmaps
        cols = 2
        rows = len(conds) // cols + 1
        fig = plt.figure(figsize=(rows*6, cols*6))
        fig.suptitle(suptitles[k], fontsize=18)
        for i, (cond_id, idxs_cond) in enumerate(conds.items()):
            
            with h5py.File(snakemake.output[0], 'r') as f:  
                cov_mx = np.array(f[group + '_' + name][cond_id]['cov_mx'])
                #labels = np.array(f[group + '_' + name][cond_id]['labels'])

            cov_mx_srt = cov_mx[np.ix_(idxs_sort, idxs_sort)]
            cov_mxs[cond_id] = cov_mx_srt

            ax = fig.add_subplot(rows, cols, i+1)
            ax.imshow(cov_mx_srt, vmin=-1, vmax=1, cmap='seismic')
            title = cond_id if not cond_id == 'all' else f'{cond_id} - {group} ({name})'
            ax.set_title(title, fontsize=14)
            for clu_no in set(labels):
                l_pos = np.where(labels[idxs_sort] == clu_no)[0][-1] + 0.5
                ax.axhline(l_pos, color='black', lw=1.5)
                ax.axvline(l_pos, color='black', lw=1.5)

            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        pdf.savefig(fig)

        # ----------------- correlations, MFR / SIC, Scree, PC1 corrs
        diffs = []
        bins = np.linspace(-1, 1, 51)
        fig, axes = plt.subplots(1, 5, figsize=(16, 5))

        # corr histogram plot
        ax = axes[0]
        for cond_id, cov_mx_srt in cov_mxs.items():
            if cond_id == 'all':
                continue
            all_corrs = cov_mx_srt.flatten()
            dist, _ = np.histogram(all_corrs, density=True, bins=bins)
            #pos = dist[25:].sum()
            #neg = dist[:25].sum()
            #diff = np.abs(pos) - np.abs(neg)
            #diffs.append(diff)
            st, pv = stats.wilcoxon(all_corrs)
            med = np.mean(all_corrs)

            ax.plot(bins[1:], dist, label=f'{cond_id} {med:.3f}')
        ax.axvline(0, ls='--', color='black')
        ax.set_xlim(-0.5, 0.5)
        ax.legend()

        # MFR / SIC
        ax = axes[1]
        for i, clu_id in enumerate(list(set(labels))):
            #sel_units = [unit for unit, label in labels_dict.items() if label in (clu_id,)]
            #unit_idxs = np.array([unit_names.index(unit_id) for unit_id in sel_units])
            idxs_clu = np.where(labels == clu_id)[0]
            ax.scatter(FR_mx[idxs_clu][:, 0], FR_mx[idxs_clu][:, 2], color=colors[i], alpha=0.7, label=str(clu_id))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('MFR, Hz', fontsize=14)
        ax.set_ylabel('SIC, bits', fontsize=14)
        ax.axhline(0.1, ls='--', color='black')
        ax.legend()
        
        # MFR / SIC - no stimulus
        ax = axes[2]
        for i, clu_id in enumerate(list(set(labels))):
            idxs_clu = np.where(labels == clu_id)[0]
            ax.scatter(FR_mx[idxs_clu][:, 0], FR_mx[idxs_clu][:, 3], color=colors[i], alpha=0.7, label=str(clu_id))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('MFR, Hz', fontsize=14)
        ax.set_ylabel('SIC (no stim.), bits', fontsize=14)
        ax.axhline(0.1, ls='--', color='black')
        ax.legend()

        # PC1 corr matrix
        ax = axes[3]
        ax.imshow(corr_mx, cmap='coolwarm', vmin=-1, vmax=1)
        ax.set_xticks(range(corr_mx.shape[0]))
        ax.set_xticklabels([key for key in eignvectors.keys()], rotation=60)
        ax.set_yticks(range(corr_mx.shape[0]))
        ax.set_yticklabels([key for key in eignvectors.keys()])

        # Scree on PCA
        ax = axes[4]
        for cond_name, exp_r in exp_ratios.items():
            ax.plot(np.arange(len(exp_r)), exp_r, linewidth=2, label=cond_name)
        ax.legend()

        fig.tight_layout()
        pdf.savefig(fig)

        # -------------- PCA on states
        idxs_states_A = [  # stationary and running
            conds['bgr_sta'],
            conds['sil_sta'],
            conds['bgr_run'],
            conds['sil_run'],
            conds['tgt_succ'],
        ]
        idxs_states_B = [  # just stationary
            conds['bgr_sta'],
            conds['sil_sta'],
            conds['tgt_succ'],
        ]

        clrs = [
            ['tab:blue', 'grey', 'red', 'black', 'tab:orange'],
            ['tab:blue', 'grey', 'tab:orange']
        ]
        ax_labels = [
            ['bgr_sta', 'sil_sta', 'bgr_run', 'sil_run', 'tgt_succ'],
            ['bgr_sta', 'sil_sta', 'tgt_succ'],
        ]
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        for i, idxs_state in enumerate([idxs_states_A, idxs_states_B]):
        
            idxs_states_flat = np.concatenate(idxs_state)
            lbls = []
            for j, idxs in enumerate(idxs_state):
                lbls += list(j*np.ones(len(idxs)))
            lbls = np.array(lbls)

            w_mx = mx_z.T[idxs_states_flat]

            tsne = TSNE(n_components=2, perplexity=30, random_state=0)
            tsne_fit = tsne.fit_transform(w_mx)

            ax = axes[i]
            for j, label in enumerate(set(lbls)):
                idxs_label = np.where(lbls == j)[0]
                ax.scatter(tsne_fit[idxs_label][:, 0], tsne_fit[idxs_label][:, 1], color=clrs[i][j], alpha=0.2, label=ax_labels[i][j])
            ax.legend()
            ax.set_xlabel('t-SNE 1', fontsize=14)
            if i == 0:
                ax.set_ylabel('t-SNE 2', fontsize=14)

        fig.tight_layout()
        pdf.savefig(fig)

        # -------------- diffs
        diff_combs = [
            ['bgr_run', 'sil_run'],
            ['bgr_sta', 'sil_sta'],
            ['bgr_sta', 'bgr_run'],
            ['sil_sta', 'sil_run'],
            ['bgr_sta', 'tgt_succ'],
            ['bgr_run', 'tgt_succ'],
            # ['bgr_sta', 'bgr_sta_AL'],
            # ['bgr_sta', 'bgr_sta_PH'],
            # ['sil_sta', 'sil_sta_AL'],
            # ['sil_sta', 'sil_sta_PH'],
            # ['bgr_sta_PH', 'tgt_succ'],
            # ['bgr_sta_PH', 'sil_sta_PH'],
            # ['bgr_sta_PH', 'bgr_run'],
            # ['sil_sta_PH', 'sil_run'],
        ]
        rows = len(diff_combs)
        cols = 5
        fig, axes = plt.subplots(rows, cols, figsize=(cols*3, rows*3))

        for i, diff_comb in enumerate(diff_combs):
            with h5py.File(snakemake.output[0], 'r') as f:  
                cov_mx_A = np.array(f[group + '_' + name][diff_comb[0]]['cov_mx'])
                cov_mx_B = np.array(f[group + '_' + name][diff_comb[1]]['cov_mx'])

            cov_mx_A_srt = cov_mx_A[np.ix_(idxs_sort, idxs_sort)]
            cov_mx_B_srt = cov_mx_B[np.ix_(idxs_sort, idxs_sort)]
            cov_mx_C_srt = cov_mx_B_srt - cov_mx_A_srt

            titles = [diff_comb[0], diff_comb[1], 'Difference']
            for j, cov_mx_srt in enumerate([cov_mx_A_srt, cov_mx_B_srt, cov_mx_C_srt]):
                ax = axes[i][j]
                ax.imshow(cov_mx_srt, vmin=vmin, vmax=vmax, cmap='seismic')
                ax.set_title(titles[j], fontsize=14)
                for clu_no in set(labels):
                    l_pos = np.where(labels[idxs_sort] == clu_no)[0][-1] + 0.5
                    ax.axhline(l_pos, color='black', lw=1.5)
                    ax.axvline(l_pos, color='black', lw=1.5)
                ax.set_xticks([])
                ax.set_yticks([])

            # eigenvectors
            v1 = eignvectors[diff_comb[0]]
            v2 = eignvectors[diff_comb[1]]
            corr, pv = stats.pearsonr(np.abs(v1), np.abs(v2))
            ax = axes[i][3]
            #ax.plot(v1 + 0.0, np.flip(np.arange(len(v1))), color='red', label=diff_comb[0])
            #ax.plot(v2 + 0.5, np.flip(np.arange(len(v2))), color='navy', label=diff_comb[1])
            #ax.plot(v2-v1 +1, np.flip(np.arange(len(v2))), color='black', label='PC1 diff')
            ax.barh(np.flip(np.arange(len(v1))), v1, left=0.0*np.ones(len(v1)), color='red', label=diff_comb[0])
            ax.barh(np.flip(np.arange(len(v1))), v2, left=0.5*np.ones(len(v1)), color='navy', label=diff_comb[1])
            ax.barh(np.flip(np.arange(len(v1))), v2-v1, left=1.0*np.ones(len(v1)), color='black', label='PC1 diff')

            ax.set_title(f'Corr.: {corr:.2f}')
            ax.set_xlim(-0.25, 1.25)
            poss = [0, 0.5, 1]
            clrs = ['red', 'navy', 'black']
            for l in range(len(poss)):
                ax.axvline(poss[l], lw=1, color=clrs[l])

            # histograms
            bins = np.linspace(-1, 1, 51)
            clrs = ['red', 'navy']
            ax = axes[i][4]
            for j, cov_mx_srt in enumerate([cov_mx_A_srt, cov_mx_B_srt]):
                all_corrs = cov_mx_srt.flatten()
                dist, _ = np.histogram(all_corrs, density=True, bins=bins)
                st, pv = stats.wilcoxon(all_corrs)
                med = np.mean(all_corrs)

                ax.plot(bins[1:], dist, label=f'{diff_comb[j]} {med:.3f}', color=clrs[j])
            ax.axvline(0, ls='--', color='black')
            ax.set_xlim(-0.5, 0.5)
            ax.legend()

        fig.tight_layout()
        pdf.savefig(fig)