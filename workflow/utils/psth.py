import numpy as np
import h5py
from scipy import stats


def get_shuffled(spiketrain):
    # shuffle spike times preserving inter-spike intervals
    ISIs = np.diff(spiketrain)
    np.random.shuffle(ISIs)
    return np.concatenate([[spiketrain[0]], spiketrain[0] + np.cumsum(ISIs)])


def get_spike_counts(spk_times, pulse_times, hw=0.25, bin_count=51):
    collected = []
    for t_pulse in pulse_times:
        selected = spk_times[(spk_times > t_pulse - hw) & (spk_times < t_pulse + hw)]
        collected += [x for x in selected - t_pulse]
    collected = np.array(collected)

    bins = np.linspace(-hw, hw, bin_count)
    counts, _ = np.histogram(collected, bins=bins)
    counts = counts / len(pulse_times) # * 1/((2. * hw)/float(bin_count - 1))
    counts = counts / (bins[1] - bins[0])  # divide by bin size to get firing rate
    
    return bins, counts


def compute_shuffled_metrics(strain, event_times, offset, bin_count, iter_count=1000):
    # psth with shuffled spike data across a given sound condition
    psth_shuffled = np.zeros((iter_count, bin_count - 1))
    for i in range(iter_count):  # shuffle 1000 times
        shuffled = get_shuffled(strain)
        bins, psth = get_spike_counts(shuffled, event_times, hw=offset, bin_count=bin_count)
        psth_shuffled[i] = psth
        
    # percentiles
    confidence_5_0_low  = np.zeros(psth_shuffled.shape[1])
    #confidence_2_5_low  = np.zeros(psth_shuffled.shape[1])
    confidence_95_0_high = np.zeros(psth_shuffled.shape[1])
    #confidence_97_5_high = np.zeros(psth_shuffled.shape[1])
    for i, col in enumerate(psth_shuffled.T):
        confidence_5_0_low[i]  = np.percentile(col, 5)
        confidence_95_0_high[i] = np.percentile(col, 95)
        #confidence_2_5_low[i]  = np.percentile(col, 2.5)
        #confidence_97_5_high[i] = np.percentile(col, 97.5)
        
    # bins, shuffled mean, shuffled std, 0.05, 0.95 percentiles (p=0.05, p=0.0)
    return np.vstack([
        bins[:-1],  
        psth_shuffled.mean(axis=0),
        psth_shuffled.std(axis=0),
        #confidence_2_5_low, 
        #confidence_97_5_high,
        confidence_5_0_low, 
        confidence_95_0_high,
    ])


def staple_pulsetrain(pulses, periods):
    # Staple pulses according to given periods as if there is no gaps.
    # pulses    - pulse times in seconds
    # periods   - list of [t_start, t_end] periods to staple
    shift = 0
    adjusted_pulses = pulses.copy()
    for period in periods:
        t_start, t_end = period[0], period[1]

        idxs = np.where((pulses > t_start) & (pulses < t_end))[0]
        adjusted_pulses[idxs] -= t_start  # align to time 0
        adjusted_pulses[idxs] += shift
        shift += t_end - t_start

    return adjusted_pulses


def staple_spike_times(s_times, periods, mode='sequence'):
    # 'sequence' - periods follow each other in a sequence
    # 'overlay'  - all periods aligned to time zero
    # returns list of lists for each given period!
    all_spikes  = []  # collect as groups
    sil_dur = 0
    for period in periods:
        idxs_tl_l, idxs_tl_r = period[0], period[1]
        t_start, t_end = period[0], period[1]

        spikes = s_times[(s_times > t_start) & (s_times < t_end)]
        spikes -= t_start  # align to time 0
        if mode == 'sequence':
            spikes += sil_dur  # adjust to already processed silence periods
        all_spikes.append(spikes)

        sil_dur += t_end - t_start
    return all_spikes  #np.array([item for sublist in all_spikes for item in sublist])


def get_psth_matrix(psth_file, electrodes):
    # build response profile (PSTH) matrix for all units all sound conditions
    # electrodes - a list of electrodes to take into account, like [1, 2]
    psths_all = {}
    psth_bins = None
    with h5py.File(psth_file, 'r') as f:
        conditions = [x for x in f]  # should be like ['BGR', 'TGT', 'SIL', 'NOI']

    for k, cond in enumerate(conditions):
        psths = []
        with h5py.File(psth_file, 'r') as f:
            for unit_name in f[cond]:
                if int(unit_name[0]) not in electrodes:
                    continue
                psths.append(np.array(f[cond][unit_name]['profile_stats'][1]))
                if psth_bins is None:
                    psth_bins = np.array(f[cond][unit_name]['profile_stats'][0])
        psths = np.array(psths)

        # z-score
        psths_z = psths.copy()
        for i, psth in enumerate(psths):
            psths_z[i] = stats.zscore(psth)
            
        psths_all[cond] = psths_z
    return psth_bins, psths_all