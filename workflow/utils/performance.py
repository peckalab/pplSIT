import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


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


def calculate_performance(tl, trial_idxs, cfg, islands=None):
    """
    Returns a matrix of time_bins x metrics, usually (12 x 7) of 
    performance_median, performance_upper_CI, performance_lower_CI, chance_median, chance_upper_CI, chance_lower_CI, time
    """
    
    arena_r = cfg['position']['floor_r_in_meters']
    target_r = cfg['experiment']['target_radius']
    t_sit = cfg['experiment']['target_duration']
    timepoints = cfg['experiment']['timepoints']
    s_duration = cfg['experiment']['session_duration']
    distractor_fail = cfg['experiment']['distractor_fail']
    distractor_islands = cfg['experiment']['distractor_islands']

    trial_time = tl[trial_idxs[:, 1].astype(np.int32)][:, 0] - tl[trial_idxs[:, 0].astype(np.int32)][:, 0]
    correct_trial = (trial_idxs[:, 5] == 1)

    time_bin_length = 5  # in secs
    N_time_slot = int(cfg['experiment']['trial_duration'] / time_bin_length)  # 12 bins
    time_x_2plot = (np.arange(N_time_slot) + 1) * time_bin_length
    amount_trials = len(trial_time)

    amount_correct = np.zeros(N_time_slot, dtype=np.int32)
    if distractor_fail:
        amount_distractor_fail = np.zeros(N_time_slot, dtype=np.int32)
    for i, t_bin in enumerate(time_x_2plot):
        amount_correct[i] = len(np.where((trial_time < t_bin)&(correct_trial))[0])
        if distractor_fail:
            amount_distractor_fail[i] = len(np.where((trial_time < t_bin)&(~correct_trial))[0])

    proportion_correct = amount_correct / amount_trials
    
    # bootstrapping real trials
    bs_count = 1000

    proportion_correct_bs = np.zeros((bs_count, N_time_slot))
    confidence_interval_real = np.zeros((2, N_time_slot))
    if distractor_fail:
        proportion_distractor_fail_bs = np.zeros((bs_count, N_time_slot))
        confidence_interval_distractor_fail = np.zeros((2, N_time_slot))

    for i in range(N_time_slot):
        for bs in range(bs_count):
            temp_index = np.random.randint(0, amount_trials, amount_trials)
            temp_correct = np.zeros(amount_trials)
            temp_correct[:amount_correct[i]] = 1
            if distractor_fail:
                temp_distractor_fail = np.zeros(amount_trials)
                temp_distractor_fail[:amount_distractor_fail[i]] = 1

            proportion_correct_bs[bs, i] = temp_correct[temp_index].sum() / float(amount_trials)
            if distractor_fail:
                proportion_distractor_fail_bs[bs, i] = temp_distractor_fail[temp_index].sum() / float(amount_trials)
        confidence_interval_real[0, i] = np.percentile(proportion_correct_bs[:, i], 97.5) - np.median(proportion_correct_bs[:, i])
        confidence_interval_real[1, i] = np.percentile(proportion_correct_bs[:, i], 2.5) - np.median(proportion_correct_bs[:, i])
        if distractor_fail:
            confidence_interval_distractor_fail[0, i] = np.percentile(proportion_distractor_fail_bs[:, i], 97.5) - np.median(proportion_distractor_fail_bs[:, i])
            confidence_interval_distractor_fail[1, i] = np.percentile(proportion_distractor_fail_bs[:, i], 2.5) - np.median(proportion_distractor_fail_bs[:, i])

    # creating list of fake islands that will not overlap with target islands
    no_fake_islands = 1000
    fake_island_centers_x = np.empty((no_fake_islands, amount_trials))
    fake_island_centers_y = np.empty((no_fake_islands, amount_trials))
    fake_island_centers_x[:] = np.nan
    fake_island_centers_y[:] = np.nan

    # Check if islands array has data when distractors are configured
    if distractor_islands > 0 and len(islands) == 0:
        print("Warning: distractor_islands configured but islands array is empty. Skipping fake island generation.")
        # Return early or set distractor_islands to 0 to skip processing
        distractor_islands = 0

    for i in range(amount_trials):
        X_target, Y_target = trial_idxs[i][2], trial_idxs[i][3]
        if distractor_islands>0:
            # take x and y of distractors, from islands array, which is tgt_x, tgt_y, tgt_r, d1_x, d1_y, d1_r, d2_x, d2_y, d2_r, d3_x, d3_y, d3_r
            X_dist = []
            Y_dist = []
            R_dist = []

            for d in range(distractor_islands):
                X_dist.append(islands[i][3 + d*3])
                Y_dist.append(islands[i][4 + d*3])
                R_dist.append(islands[i][5 + d*3])

        count = 0

        # this block was replaced by empirical sampling with maximum number of attempts
        # d_target = np.hypot(X_target, Y_target)
        # R1 = arena_r - target_r
        # R2 = 2 * target_r
        # # if R2 is larger than R1 + d_target, the two discs don't overlap
        # if R2 > R1 + d_target:
        #     # no valid fake islands exist; skip this trial
        #     print(f"Skipping trial {i} due to no valid fake islands.")
        #     continue

        max_attempts = 100000
        attempts = 0

        while np.isnan(fake_island_centers_x[:, i]).any() and attempts < max_attempts:
            attempts += 1
            angle = 2 * np.pi * np.random.rand()
            r = arena_r * np.sqrt(np.random.rand())
            x_temp = r * np.cos(angle)  # add center of the arena if not centered
            y_temp = r * np.sin(angle)  # add center of the arena if not centered

            # check target overlap
            valid = np.sqrt((x_temp - X_target)**2 + (y_temp - Y_target)**2) > 2 * target_r

            # check arena boundary
            valid &= np.sqrt(x_temp**2 + y_temp**2) < arena_r - target_r

            # check distractor overlap (if they exist)
            if distractor_islands > 0:
                for xd, yd, rd in zip(X_dist, Y_dist, R_dist):
                    if np.sqrt((x_temp - xd)**2 + (y_temp - yd)**2) <= (target_r + rd):
                        valid = False
                        break

            if valid:
                fake_island_centers_x[count, i] = x_temp
                fake_island_centers_y[count, i] = y_temp
                count += 1
        if attempts == max_attempts:
            print(f"Skipping trial {i}: could not place all fake islands without overlap with target (or distractors).")
            continue
    
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

        confidence_interval_bs_fake[0, time_slot] = np.percentile(proportion_correct_bs_fake[:, time_slot], 97.5) - np.median(proportion_correct_bs_fake[:, time_slot])
        confidence_interval_bs_fake[1, time_slot] = np.percentile(proportion_correct_bs_fake[:, time_slot], 2.5) - np.median(proportion_correct_bs_fake[:, time_slot])
        
    # compute performance metrics
    c_median = 100 * np.median(proportion_correct_bs_fake, axis=0)
    c_lower_CI = 100 * confidence_interval_bs_fake[1]
    c_upper_CI = 100 * confidence_interval_bs_fake[0]

    p_median = 100 * np.median(proportion_correct_bs, axis=0)
    p_lower_CI = 100 * confidence_interval_real[1]
    p_upper_CI = 100 * confidence_interval_real[0]

    if distractor_fail:
        d_median = 100 * np.median(proportion_distractor_fail_bs, axis=0)
        d_lower_CI = 100 * confidence_interval_distractor_fail[1]
        d_upper_CI = 100 * confidence_interval_distractor_fail[0]
    
    if distractor_fail:
        return np.column_stack([time_x_2plot, c_median, c_lower_CI, c_upper_CI, p_median, p_lower_CI, p_upper_CI, d_median, d_lower_CI, d_upper_CI])
    return np.column_stack([time_x_2plot, c_median, c_lower_CI, c_upper_CI, p_median, p_lower_CI, p_upper_CI])


def plot_session_metrics(tl, trial_idxs, cfg, fig_path):
    arena_r = cfg['position']['floor_r_in_meters']

    fig = plt.figure(figsize=(12, 12))

    # trajectory and islands
    ax = fig.add_subplot(221)
    ax.scatter(tl[:, 1], tl[:, 2], s=1, alpha=0.1)  # positions
    scat = ax.scatter(trial_idxs[trial_idxs[:,-1]==0, 2], trial_idxs[trial_idxs[:,-1]==0, 3], s=1000, facecolors='none', edgecolors='r')  # incorrect islands, radius approx.
    scat = ax.scatter(trial_idxs[trial_idxs[:,-1]==1, 2], trial_idxs[trial_idxs[:,-1]==1, 3], s=1000, facecolors='none', edgecolors='g')  # incorrect islands, radius approx.
    ax.add_patch(plt.Circle((0, 0), arena_r, color='r', fill=False))
    ax.set_aspect('equal')
    ax.set_xlabel('X, m', fontsize=14)
    ax.set_ylabel('Y, m', fontsize=14)
    ax.set_title('Running', fontsize=14)
    ax.grid()

    # occupancy
    sigma = 0.1
    lin_profile = np.linspace(-15, 15, 20)
    bump = np.exp(-sigma * lin_profile**2)
    bump /= np.trapz(bump)  # normalize to 1
    kernel = bump[:, np.newaxis] * bump[np.newaxis, :]
    occupancy_map, _, _ = np.histogram2d(tl[:, 1], tl[:, 2], bins=[40, 40], range=np.array([[-0.5, 0.5], [-0.5, 0.5]]))
    occupancy_map = signal.convolve2d(occupancy_map, kernel, mode='same')

    ax = fig.add_subplot(222)
    ax.imshow(occupancy_map.T, origin='lower', extent=(-0.5, 0.5, -0.5, 0.5), cmap='Blues')
    ax.add_patch(plt.Circle((0, 0), arena_r, color='r', fill=False))
    ax.set_xlabel('X, m', fontsize=14)
    ax.set_title('Occupancy', fontsize=14)
    ax.grid()

    # trials
    durations = tl[trial_idxs[:, 1].astype(int)][:, 0] - tl[trial_idxs[:, 0].astype(int)][:, 0]
    colors = ['red' if x == 1 else 'grey' for x in trial_idxs[:, 5]]

    ax = fig.add_subplot(223)
    ax.barh(np.arange(len(trial_idxs)), durations, color=colors, align='center')
    ax.set_xlabel('Time, s', fontsize=14)
    ax.set_ylabel('Trial, #', fontsize=14)
    ax.set_title('Trials', fontsize=14)

    # speed
    ax = fig.add_subplot(224)

    s_rate = 100  # Hz
    window = 60   # secs
    step = 10     # secs
    duration = tl[-1][0]
    x_vals = np.arange(int(duration/step))

    inst_speed = [tl[x*step*s_rate:(x*step + window)*s_rate][:, 3].mean() for x in x_vals]
    ax.plot(x_vals*step, inst_speed)
    ax.set_ylabel('Speed, m/s', fontsize=14)
    ax.set_xlabel('Time, s', fontsize=14)
    ax.set_title('Speed', fontsize=14)

    fig.tight_layout()
    fig.savefig(fig_path)

def plot_performance(cfg, perf, fig_path):
    fig = plt.figure(figsize=(4, 4))

    # check distractor fail
    distractor_fail = cfg['experiment']['distractor_fail']
    x = perf[:, 0] # time points

    ax = fig.add_subplot(111)

    ax.plot(x, perf[:, 4], color='C0', label='Performance')  # performance
    ax.plot(x, perf[:, 1], color='C1', label='Chance')  # chance
    ax.fill_between(x, perf[:, 4] + perf[:, 5], perf[:, 4] + perf[:, 6], alpha=0.4, color='C0') # performance CI
    ax.fill_between(x, perf[:, 1] + perf[:, 2], perf[:, 1] + perf[:, 3], alpha=0.4, color='C1') # chance CI

    if distractor_fail:
        ax.plot(x, perf[:, 7], color='C2', label='Distractor Fail')  # distractor fail
        ax.fill_between(x, perf[:, 7] + perf[:, 8], perf[:, 7] + perf[:, 9], alpha=0.4, color='C2') # distractor fail CI
    ax.legend(fontsize=10, loc='upper left')
    ax.set_ylim(0, 110)
    ax.set_xlim(0, 65)
    ax.grid()
    # ax.set_title(experiment_id[-19:], fontsize=14)
    ax.set_xlabel('Time, s', fontsize=14)

    ax.set_ylabel('Successful trials, %', fontsize=14)
            
    fig.tight_layout()
    fig.savefig(fig_path)