import h5py, os
import numpy as np
from scipy import signal

"""
MoSeq bodyparts:

- left_eye
- right_eye
- left_ear
- right_ear
- nose
- red_dot
- tail_base
- lower_spine
- neck
- left_shoulder
- right_shoulder
- left_hip
- right_hip
"""


# load data
with h5py.File(snakemake.input[0], 'r') as f:
    moseq = np.array(f['moseq'])
    headers = list(f['moseq'].attrs['headers'])

# compute speed from MoSeq centroid, smooth with 1 sec gaussian
width = 75  # 100 points ~= 1 sec with at 100Hz
kernel = signal.gaussian(width, std=(width) / 7.2)
dx = np.sqrt(np.square(np.diff(moseq[:, 3])) + np.square(np.diff(moseq[:, 4])))
dt = np.diff(moseq[:, 0])
speed = np.concatenate([dx/dt, [dx[-1]/dt[-1]]])
speed_smooth = np.convolve(speed, kernel, 'same') / kernel.sum()

# compute head direction as Neck-to-Nose vector (10, 11 Nose, 18, 19 Neck)
diff_x = moseq[:, 10] - moseq[:, 18]
diff_y = moseq[:, 11] - moseq[:, 19]

hd = np.arctan2(diff_y, diff_x)
hd = np.concatenate([np.array([hd[0]]), hd])  # make the same length as timeline

# write out
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('speed', data=speed_smooth)
    f.create_dataset('hd', data=hd)