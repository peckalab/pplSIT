import os, json
import numpy as np


def _normalize_bombcell_label(label):
    return str(label).strip().upper()


def load_bombcell_unit_labels(ks_path):
    """
    Read Bombcell unit type labels written for Phy compatibility.
    Returns {cluster_id: label}.
    """
    tsv_path = os.path.join(ks_path, "cluster_bc_unitType.tsv")
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(
            f"Bombcell labels requested, but cluster_bc_unitType.tsv is missing: {tsv_path}"
        )

    with open(tsv_path, "r") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        if "cluster_id" not in col_idx:
            raise ValueError(f"Missing cluster_id column in {tsv_path}")

        label_col = None
        for candidate in ("bc_unitType", "unitType", "bc_unit_type"):
            if candidate in col_idx:
                label_col = col_idx[candidate]
                break
        if label_col is None:
            raise ValueError(
                f"Missing Bombcell unit type column in {tsv_path}. "
                "Expected one of: bc_unitType, unitType, bc_unit_type."
            )

        labels = {}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            cluster_id = int(parts[col_idx["cluster_id"]])
            labels[cluster_id] = _normalize_bombcell_label(parts[label_col])

    return labels


def filter_cluster_ids_by_bombcell(ks_path, allowed_labels):
    allowed = {_normalize_bombcell_label(label) for label in allowed_labels}
    labels = load_bombcell_unit_labels(ks_path)
    return np.array(
        [cluster_id for cluster_id, label in labels.items() if label in allowed],
        dtype=int
    )


def load_ks_units_before(ks_path, label_source="ks", bombcell_unit_types=None):
    """
    Load units from Kilosort BEFORE manual curation.
    Uses KSLabel == 'good' to load clusters.
    """

    # ---- probe geometry ----
    with open(os.path.join(ks_path, 'probe.json'), 'r') as json_file:
        probe = json.load(json_file)

    ch_shank_map = np.array(probe['kcoords'], dtype=np.int16) + 1
    shanks = np.unique(ch_shank_map)

    # ---- load Kilosort outputs ----
    s_times   = np.load(os.path.join(ks_path, 'spike_times.npy'))
    s_clust   = np.load(os.path.join(ks_path, 'spike_clusters.npy'))
    templates = np.load(os.path.join(ks_path, 'templates.npy'))
    ch_pos    = np.load(os.path.join(ks_path, 'channel_positions.npy'))

    if label_source == "bombcell":
        good_idxs = filter_cluster_ids_by_bombcell(
            ks_path,
            bombcell_unit_types or ["GOOD", "NON-SOMA GOOD"]
        )
    else:
        # ---- read cluster_KSLabel.tsv without pandas ----
        good_idxs = []
        tsv_path = os.path.join(ks_path, 'cluster_KSLabel.tsv')

        with open(tsv_path, 'r') as f:
            header = f.readline().strip().split('\t')
            cluster_col = header.index('cluster_id') if 'cluster_id' in header else 0
            label_col   = header.index('KSLabel')

            for line in f:
                parts = line.strip().split('\t')
                cluster_id = int(parts[cluster_col])
                label = parts[label_col]
                if label == 'good':
                    good_idxs.append(cluster_id)

        good_idxs = np.array(good_idxs, dtype=int)

    # ---- template peak channel mapping ----
    template_maxchans = np.abs(templates).max(axis=1).argmax(axis=1)
    clu_ch_mapping = ch_shank_map[template_maxchans]
    clu_pos_mapping = ch_pos[template_maxchans]

    # ---- organize output by shank ----
    all_units = {}
    all_pos = {}

    for shank in shanks:
        clu_idxs = np.where(clu_ch_mapping == shank)[0]
        sel_clusters = np.intersect1d(good_idxs, clu_idxs)

        spiketrains = {}
        sel_pos = {}

        for clu_id in sel_clusters:
            spiketrains[clu_id] = s_times[s_clust == clu_id]
            sel_pos[clu_id] = clu_pos_mapping[clu_id]

        all_units[shank] = spiketrains
        all_pos[shank] = sel_pos

    return all_units, all_pos


def load_ks_units_after(ks_path, label_source="ks", bombcell_unit_types=None):
    """
    Load kilosorted units AFTER manual curation.
    Returns:
        all_units[shank][cluster_id] = spike_times
        unit_info[shank][cluster_id] = [Amplitude, ContamPct, amp, ch, depth, fr, label]
    """
    def clean_shank(sh):
        return int(float(sh))

    # ---- load spike data ----
    s_times = np.load(os.path.join(ks_path, 'spike_times.npy'))
    s_clust = np.load(os.path.join(ks_path, 'spike_clusters.npy'))

    # ---- read cluster_info.tsv manually ----
    tsv_path = os.path.join(ks_path, 'cluster_info.tsv')

    with open(tsv_path, 'r') as f:
        header = f.readline().strip().split('\t')
        col_idx = {name: i for i, name in enumerate(header)}

        records = []
        for line in f:
            parts = line.strip().split('\t')
            records.append(parts)

    # ---- collect unique shanks ----
    shank_idx = col_idx['sh']
    shanks = sorted(set(clean_shank(r[shank_idx]) for r in records))

    all_units = {}
    unit_info = {}
    bombcell_good_ids = None
    if label_source == "bombcell":
        bombcell_good_ids = set(filter_cluster_ids_by_bombcell(
            ks_path,
            bombcell_unit_types or ["GOOD", "NON-SOMA GOOD"]
        ).tolist())

    for shank in shanks:
        spiketrains = {}
        u_info_sh = {}

        for r in records:
            if clean_shank(r[shank_idx]) != shank:
                continue

            ks_label = r[col_idx['KSLabel']]
            group = r[col_idx['group']]
            clu_id = int(r[col_idx['cluster_id']])

            if label_source == "bombcell":
                if clu_id not in bombcell_good_ids:
                    continue
            else:
                # Apply same filtering logic as original
                if not ((ks_label == 'good') or (group == 'good')):
                    continue
                if group == 'noise':
                    continue

            # spike times
            spiketrains[clu_id] = s_times[s_clust == clu_id]

            # metadata
            u_info_sh[clu_id] = np.array([
                float(r[col_idx['Amplitude']]),
                float(r[col_idx['ContamPct']]),
                float(r[col_idx['amp']]),
                int(r[col_idx['ch']]),
                float(r[col_idx['depth']]),
                float(r[col_idx['fr']]),
                1 if ks_label == 'good' else 2,
            ])

        all_units[int(shank) + 1] = spiketrains
        unit_info[int(shank) + 1] = u_info_sh

    return all_units, unit_info


def infer_n_chan_bin(dat_path: str, base_n_chan: int = 384, dtype_bytes: int = 2) -> int:
    """
    Infer number of channels in a raw binary .dat by checking file size divisibility.

    We assume:
      - int16 => 2 bytes/sample/channel
      - real probe channels are base_n_chan (default 384)
      - optional extra sync channel => base_n_chan + 1
    """
    size = os.path.getsize(dat_path)
    for n in (base_n_chan, base_n_chan + 1):
        if size % (dtype_bytes * n) == 0:
            return n
    raise ValueError(
        f"Cannot infer n_chan_bin for {dat_path}. "
        f"File size {size} not divisible by 2*{base_n_chan} or 2*{base_n_chan+1}."
    )
