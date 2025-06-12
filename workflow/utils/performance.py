import numpy as np


def when_successful(traj, x_isl, y_isl, r_isl, t_sit):
    """
    traj - t, x, y - a matrix Nx3 of position data, equally sampled!
    """
    splits = np.where( (traj[:, 1] - x_isl)**2 + (traj[:, 2] - y_isl)**2 < r_isl**2 )[0]
    df = np.where(np.diff(splits) > 5)[0]  # idxs of periods of starts  

    if len(splits) == 0:
        return None
    
    periods = [[0, df[0] if len(df) > 0 else len(splits)-1]]
    if len(df) > 1:
        for point in df[1:]:
            periods.append( [periods[-1][1] + 1, point] )

    if len(df) > 0:
        periods.append([periods[-1][1] + 1, len(splits)-1])

    for period in periods:
        if splits[period[0]] - 5 < 0:
            continue
        if traj[splits[period[1]]][0] - traj[splits[period[0]] - 5][0] > t_sit:  # -5 is a hack
            return traj[splits[period[0]]][0] + t_sit


def calculate_performance(tl, trial_idxs, cfg):
    """
    Returns a matrix of time_bins x metrics, usually (12 x 7) of 
    performance_median, performance_upper_CI, performance_lower_CI, chance_median, chance_upper_CI, chance_lower_CI, time
    """
    
    arena_r = cfg['position']['floor_r_in_meters']
    target_r = cfg['experiment']['target_radius']
    t_sit = cfg['experiment']['target_duration']
    timepoints = cfg['experiment']['timepoints']
    s_duration = cfg['experiment']['session_duration']

    trial_time = tl[trial_idxs[:, 1].astype(np.int32)][:, 0] - tl[trial_idxs[:, 0].astype(np.int32)][:, 0]
    correct_trial = (trial_idxs[:, 5] == 1)

    time_bin_length = 5  # in secs
    N_time_slot = int(cfg['experiment']['trial_duration'] / time_bin_length)  # 12 bins
    time_x_2plot = (np.arange(N_time_slot) + 1) * time_bin_length
    amount_trials = len(trial_time)

    amount_correct = np.zeros(N_time_slot, dtype=np.int32)
    for i, t_bin in enumerate(time_x_2plot):
        amount_correct[i] = len(np.where(trial_time < t_bin)[0])

    proportion_correct = amount_correct / amount_trials
    
    # bootstrapping real trials
    bs_count = 1000

    proportion_correct_bs = np.zeros((bs_count, N_time_slot))
    confidence_interval_real = np.zeros((2, N_time_slot))

    for i in range(N_time_slot):
        for bs in range(bs_count):
            temp_index = np.random.randint(0, amount_trials, amount_trials)
            temp_correct = np.zeros(amount_trials)
            temp_correct[:amount_correct[i]] = 1

            proportion_correct_bs[bs, i] = temp_correct[temp_index].sum() / float(amount_trials)
        confidence_interval_real[0, i] = np.percentile(proportion_correct_bs[:, i], 84.1) - np.median(proportion_correct_bs[:, i])
        confidence_interval_real[1, i] = np.percentile(proportion_correct_bs[:, i], 15.9) - np.median(proportion_correct_bs[:, i])

    # creating list of fake islands that will not overlap with target islands
    no_fake_islands = 1000
    fake_island_centers_x = np.empty((no_fake_islands, amount_trials))
    fake_island_centers_y = np.empty((no_fake_islands, amount_trials))
    fake_island_centers_x[:] = np.nan
    fake_island_centers_y[:] = np.nan

    for i in range(amount_trials):
        X_target, Y_target = trial_idxs[i][2], trial_idxs[i][3]
        count = 0

        d_target = np.hypot(X_target, Y_target)
        R1 = arena_r - target_r
        R2 = 2 * target_r
        # if R2 is larger than R1 + d_target, the two discs don't overlap
        if R2 > R1 + d_target:
            # no valid fake islands exist; skip this trial
            print(f"Skipping trial {i} due to no valid fake islands.")
            continue

        while np.isnan(fake_island_centers_x[:, i]).any():
            angle = 2 * np.pi * np.random.rand()
            r = arena_r * np.sqrt(np.random.rand())
            x_temp = r * np.cos(angle)  # add center of the arena if not centered
            y_temp = r * np.sin(angle)  # add center of the arena if not centered

            if np.sqrt((x_temp - X_target)**2 + (y_temp - Y_target)**2) > 2 * target_r and \
                np.sqrt(x_temp**2 + y_temp**2) < arena_r - target_r:
                fake_island_centers_x[count, i] = x_temp
                fake_island_centers_y[count, i] = y_temp
                count += 1        
                
    # surrogate islands work, now calculate the chance performance
    surrogate_correct = np.zeros((no_fake_islands, amount_trials))

    pos_downsample = 10  # think about reducing

    for trial in range(amount_trials):
        temp_index = np.arange(trial_idxs[trial][0], trial_idxs[trial][1], pos_downsample).astype(np.int32)
        temp_traj = tl[temp_index]
        temp_traj[:, 0] -= temp_traj[0][0]  # time relative to trial start

        for surr in range(no_fake_islands):
            x_fake, y_fake = fake_island_centers_x[surr, trial], fake_island_centers_y[surr, trial]
            fake_island_time_finish = when_successful(temp_traj, x_fake, y_fake, target_r, t_sit)

            if fake_island_time_finish is not None:
                surrogate_correct[surr, trial] = fake_island_time_finish

    # now i have to do the same curve as in the real correct, but for the matrix surrogate_correct
    surr_for_deleting = np.array(surrogate_correct)
    proportion_correct_bs_fake = np.zeros((bs_count, N_time_slot))
    confidence_interval_bs_fake = np.zeros((2, N_time_slot))

    for time_slot in range(N_time_slot):
        fake_trials_to_remove = np.where(trial_time < (time_slot + 1) * time_bin_length)[0]

        for trial in fake_trials_to_remove:
            #idxs = np.logical_or(surrogate_correct[:, trial] == 0, surrogate_correct[:, trial] > trial_time[trial])
            idxs = np.where( (surrogate_correct[:, trial] == 0) | (surrogate_correct[:, trial] > trial_time[trial]) )[0]
            for idx in idxs:
                surr_for_deleting[idx, trial] = np.nan

        scwr = surr_for_deleting.flatten()
        scwr = scwr[~np.isnan(scwr)]

        for bs in range(bs_count):
            temp_index = np.random.randint(0, len(scwr), amount_trials)

            temp_correct = np.logical_and(scwr[temp_index] < (time_slot + 1) * time_bin_length, scwr[temp_index] > 0)
            proportion_correct_bs_fake[bs, time_slot] = temp_correct.sum() / float(amount_trials)

        confidence_interval_bs_fake[0, time_slot] = np.percentile(proportion_correct_bs_fake[:, time_slot], 84.1) - np.median(proportion_correct_bs_fake[:, time_slot])
        confidence_interval_bs_fake[1, time_slot] = np.percentile(proportion_correct_bs_fake[:, time_slot], 15.9) - np.median(proportion_correct_bs_fake[:, time_slot])
        
    # compute performance metrics
    c_median = 100 * np.median(proportion_correct_bs_fake, axis=0)
    c_lower_CI = 100 * confidence_interval_bs_fake[1]
    c_upper_CI = 100 * confidence_interval_bs_fake[0]

    p_median = 100 * np.median(proportion_correct_bs, axis=0)
    p_lower_CI = 100 * confidence_interval_real[1]
    p_upper_CI = 100 * confidence_interval_real[0]
    
    return np.column_stack([p_median, p_lower_CI, p_upper_CI, c_median, c_lower_CI, c_upper_CI, time_x_2plot])
