import numpy as np
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
