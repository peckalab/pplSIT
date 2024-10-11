import os, sys
import numpy as np
import h5py, json
import cv2

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.color import rgba2rgb

# import util functions from utils module - if needed
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)
#from utils.aeps import compute_metric, AEP_metrics_lims, outlier_lims, AEP_metrics_methods


# loading data
with h5py.File(snakemake.input[0], 'r') as f:
    tl     = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    events = np.array(f['processed']['sound_events'])
    trials = np.array(f['processed']['trial_idxs'])

# TODO: read unit_names from clustered corr matrix?
with open(snakemake.input[2], 'r') as f:
    units_data  = json.load(f)
    unit_names  = units_data[0]  # order matters!
    unit_labels = units_data[1]  # order matters!

single_units, spike_times = {}, {}
#with h5py.File(snakemake.input[1], 'r') as f:
#    unit_names = [x for x in f]
with h5py.File(snakemake.input[1], 'r') as f:
    for unit_name in unit_names:
        spike_times[unit_name]  = np.array(f[unit_name]['spike_times'])
        single_units[unit_name] = np.array(f[unit_name]['inst_rate'])

# timeline indices
idxs_target = np.where(tl[:, 6] == 2)[0]
idxs_backgr = np.where(tl[:, 6] == 1)[0]
idxs_noise  = np.where(tl[:, 6] ==-1)[0]
#idxs_reward = trials[trials[:, 5] == 1][:, 1].astype(np.int32)  # another way
idxs_reward = tgt_mx[tgt_mx[:, 4] == 1][:, 3]

# events
idxs_to_idx = np.where(np.diff(idxs_backgr) > 5)[0] + 1
idxs_bgr_start = idxs_backgr[:-1][idxs_to_idx]
idxs_to_idx = np.where(np.diff(idxs_noise) > 5)[0] + 1
idxs_nos_start = np.concatenate([[idxs_noise[0]], idxs_noise[:-1][idxs_to_idx]])
idxs_tgt_start = tgt_mx[:, 2]


# fragments of experimental activity
out_path = snakemake.output[0]
time_range = os.path.basename(out_path).split('.')[0]  # in seconds
t_l = int(time_range.split('_')[0])
t_r = int(time_range.split('_')[1])
u_id_min, u_id_max = 0, len(unit_names)
u_id_diff = u_id_max - u_id_min

colors = list((plt.rcParams['axes.prop_cycle'].by_key()['color']))
colors = colors + colors

fig = plt.figure(figsize=(15, 15))#, frameon=False, dpi=100)
gs = fig.add_gridspec(2, 1, height_ratios=(20, 1),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

# units
ax1 = fig.add_subplot(gs[0])
for i, unit_name in enumerate(unit_names[u_id_min:u_id_max]):
    i_rate = single_units[unit_name]
    s_times = spike_times[unit_name]
    s_vals = np.random.rand(len(s_times))
    clr = colors[unit_labels[i] - 1]
    
    idxs_ts = np.where((s_times > t_l) & (s_times < t_r))[0]
    idxs_tl = np.where((tl[:, 0] > t_l) & (tl[:, 0] < t_r))[0]
    
    ax1.scatter(s_times[idxs_ts], s_vals[idxs_ts] + u_id_diff - (i+1), s=1, color=clr)
    ax1.axhline(u_id_diff - i, color=clr, lw=1)
    
    #i_max = i_rate[idxs_tl].max()
    #ax.plot(tl[:, 0][idxs_tl], i_rate[idxs_tl]/i_max + (i+1), lw=1)
    
ax1.set_xlim(t_l, t_r)
ax1.set_ylim(1, u_id_diff)
ax1.set_yticks(np.arange(len(unit_names[u_id_min:u_id_max])) + 0.5)
_ = ax1.set_yticklabels(list(reversed(unit_names[u_id_min:u_id_max])), fontsize=10)
#ax1.margins(x=0)

# experimental timeline 
ax2 = fig.add_subplot(gs[1], sharex=ax1)
ax2.scatter(tl[idxs_target][:, 0], 2*np.ones(len(idxs_target)), s=1, color=colors[1])
ax2.scatter(tl[idxs_backgr][:, 0], 1*np.ones(len(idxs_backgr)), s=1, color=colors[0])
ax2.scatter(tl[idxs_noise][:, 0],  0*np.ones(len(idxs_noise)),  s=1, color=colors[3])
ax2.scatter(tl[idxs_reward][:, 0], 2 * np.ones(len(idxs_reward)), s=50, color=colors[2])
ax2.set_ylim(-0.5, 2.5)
#ax2.margins(x=0)
#fig.tight_layout()

# vertical lines
v_min, v_max = -0.5, u_id_diff
lines_to_plot = [idxs_bgr_start, idxs_tgt_start, idxs_reward, idxs_nos_start]
for j, idxs_var in enumerate(lines_to_plot):
    xy1 = np.vstack([ tl[idxs_var][:, 0], v_max * np.ones(len(idxs_var)) ]).T
    xy3 = np.vstack([ tl[idxs_var][:, 0], v_min * np.ones(len(idxs_var)) ]).T
    for i in range(len(idxs_var)):
        if xy1[i][0] < t_l or xy1[i][0] > t_r:
            continue
        con = ConnectionPatch(xyA=xy1[i], coordsA=ax1.transData, xyB=xy3[i], coordsB=ax2.transData, color=colors[j])
        fig.add_artist(con)

# ----- time animation line
# https://matplotlib.org/stable/gallery/animation/multiple_axes.html#sphx-glr-gallery-animation-multiple-axes-py
# https://stackoverflow.com/questions/31252107/how-to-draw-vertical-lines-interactively-in-matplotlib
anim_line = ConnectionPatch(xyA=[t_l, v_max], coordsA=ax1.transData, xyB=[t_l, v_min], \
                            coordsB=ax2.transData, color='black')
anim_line.set(lw=2, ls='--')
fig.add_artist(anim_line)


# save frames
video_path = os.path.dirname(snakemake.input[3])
frame_path = os.path.join(video_path, time_range)
if not os.path.exists(frame_path):
    os.makedirs(frame_path)

# https://www.futurelearn.com/info/courses/introduction-to-image-analysis-for-plant-phenotyping/0/steps/305359
# open video stream 
cap = cv2.VideoCapture(snakemake.input[3])
fps = cap.get(cv2.CAP_PROP_FPS)

for i in range(int((t_r - t_l)*fps + 1)):  # iteration over frames
    filename = os.path.join(frame_path, '%05d.png' % i)
    if os.path.exists(filename) and not snakemake.config['video']['overwrite']:
        continue

    t = i/fps
    anim_line.xy1 = (t + t_l, v_max)
    anim_line.xy2 = (t + t_l, v_min)
    fig.savefig(filename, dpi=100, bbox_inches='tight')
    
    # cut to even resolution, make RGB
    img = plt.imread(filename)
    if not img.shape[0] % 2 == 0:
        img = img[:img.shape[0]-1, :, :]
    if not img.shape[1] % 2 == 0:
        img = img[:, :img.shape[1]-1, :]
    if img.shape[2] > 3:
        img = rgba2rgb(img)
    
    #int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    frame_id = i + int(fps * t_l)  # global frame id
    cap.set(cv2.CAP_PROP_POS_FRAMES, frame_id)
    ret, frame = cap.read()
    
    result = np.hstack([frame/255., img])
    plt.imsave(filename, result)
    
    print('\rFrame at time %.2f rendered' % t, end='')
    
# compile video from all frames
ff_path = os.path.join(frame_path, "%05d.png")
os.system("ffmpeg -r " + str(fps) + " -i " + ff_path + " -vcodec mpeg4 -y " + snakemake.output[0])

# delete frames