"""Legacy BGR/TGT micro-PSTH layout with saved Bayesian-block timing markers."""
import json
from pathlib import Path

import h5py
import numpy as np


def plot_psth_micro_overlay(meta_path, units_path, responses_path, output_pdf,
                            micro_config, phase_lock_path=None, preview_png=None):
    """Use the existing PSTH counts/bins and unit order; never refit responses.

    Bars retain all logged BGR/TGT events, as in psth_micro.py. Markers use
    singleton sound-ID conditions from the saved Bayesian-block analysis.
    Counts and timing availability are stated explicitly in the plot.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from ..psth import get_spike_counts
    from ..neurosuite import get_unit_names_sorted
    from ..maths import pval2text

    with h5py.File(meta_path) as f:
        parameters = json.loads(f['processed'].attrs['parameters'])
        events = f['processed/sound_events'][:]
    sounds = parameters['sound']['sounds']
    if 'background' not in sounds or 'target' not in sounds:
        raise ValueError('BGR/TGT overlay requires background and target sound definitions')
    hw = float(micro_config['latency'])
    bins_count = int(micro_config['bin_count'])
    xlim = micro_config.get('xlim', [-hw, hw])
    colors = {1: 'tab:blue', 2: 'tab:orange'}
    labels = {1: 'BGR', 2: 'TGT'}
    ids = (1, 2)
    times = {sid: events[events[:, 1] == sid, 0] for sid in ids}
    phase = {}
    if phase_lock_path and Path(phase_lock_path).is_file():
        with h5py.File(phase_lock_path) as f:
            for label in ('bgr', 'tgt'):
                if label in f:
                    phase[label] = {unit: (float(g['MRL_real'][()]), float(g['p_value'][()]))
                                    for unit, g in f[label].items()}

    with h5py.File(responses_path) as response, h5py.File(units_path) as units:
        if not np.array_equal(response['sound_events'][:], events):
            raise ValueError('PSTH metadata events differ from the Bayesian-block input events')
        conditions = {}
        for group in response['conditions'].values():
            sound_ids = np.asarray(group.attrs['sound_ids']).ravel().tolist()
            if len(sound_ids) == 1 and sound_ids[0] in ids:
                sid = sound_ids[0]
                if sid in conditions:
                    raise ValueError(f'Multiple Bayesian-block conditions for sound ID {sid}')
                conditions[sid] = group
        window = json.loads(response.attrs['effective_settings'])['window_s']
        names = get_unit_names_sorted(list(units))
        depth = {name: float(units[name]['kilosort_info'][4]) if 'kilosort_info' in units[name]
                 else float(units[name]['anatomical_position'][1]) for name in names}
        names.sort(key=lambda name: (int(name.split('-')[0]), depth[name]))
        cols = 3
        rows = max(1, int(np.ceil(len(names) / cols)))
        fig = plt.figure(figsize=(4 * cols, 4 * rows + 1.1))
        caption = []
        for sid in ids:
            used = int(conditions[sid].attrs['n_events_used']) if sid in conditions else 0
            caption.append(f'{labels[sid]}: bars {len(times[sid])}, BB {used}')
        fig.suptitle(f'{Path(meta_path).parent.name}\n'
                     'BGR / TGT micro-PSTH with Bayesian-block response timing\n'
                     + '; '.join(caption) + f' | BB window: {window * 1000:g} ms',
                     fontsize=12, y=1 - .1 / fig.get_figheight())
        header = [Line2D([0], [0], color='black', ls='--', label='Response onset'),
                  Line2D([0], [0], color='black', ls=':', label='Response offset')]
        fig.legend(handles=header, loc='upper center', ncol=2,
                   bbox_to_anchor=(.5, 1 - .65 / fig.get_figheight()), frameon=False)
        preview_axes = []
        for index, name in enumerate(names):
            ax = fig.add_subplot(rows, cols, index + 1)
            if index < 3:
                preview_axes.append(ax)
            spikes = units[name]['spike_times'][:]
            availability = []
            for sid in ids:
                edges, rates = get_spike_counts(spikes, times[sid], hw=hw, bin_count=bins_count)
                label = labels[sid].lower()
                if name in phase.get(label, {}):
                    mrl, pv = phase[label][name]
                    label += f' ({mrl:.2f}; {pval2text(pv)})'
                ax.hist(edges[:-1], bins=edges, weights=rates, edgecolor='black',
                        color=colors[sid], alpha=.7, label=label)
                group = conditions.get(sid)
                if group is None or name not in group['units']:
                    availability.append(f'{labels[sid]}: not analyzed')
                    continue
                attrs = group['units'][name].attrs
                timing = []
                for kind, style in [('onset', '--'), ('offset', ':')]:
                    value = float(attrs[f'{kind}_latency_s'])
                    if np.isfinite(value):
                        ax.axvline(value, color=colors[sid], ls=style, lw=1.8, zorder=5)
                        timing.append(f'{kind} {value * 1000:.1f}')
                    else:
                        timing.append(f'{kind}: {attrs[f"{kind}_status"]}')
                availability.append(f'{labels[sid]} ' + ', '.join(timing))
            ax.axvline(0, color='black', ls='--', lw=1)
            duration = float(sounds['background']['duration'])
            ax.axvspan(0, duration, alpha=.3, color='gray')
            ax.axvspan(-hw, -hw + duration, alpha=.3, color='gray')
            ax.set_title(name, fontsize=14)
            ax.set_xlim(*xlim)
            ax.set_xlabel('Time from sound onset (s)')
            if index % cols == 0:
                ax.set_ylabel('Firing Rate, Hz', fontsize=14)
            ax.legend(loc='lower right', prop={'size': 10})
            ax.text(.98, .98, '\n'.join(availability) + '\nTiming values in ms',
                    transform=ax.transAxes, ha='right', va='top', fontsize=7,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.8))
            ax.margins(y=.22)
        if not names:
            fig.text(.5, .5, 'No units available', ha='center')
        fig.tight_layout(rect=(0, 0, 1, 1 - 1.0 / fig.get_figheight()))
        Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_pdf)
        if preview_png and preview_axes:
            # Render just the first row for inspection; the PDF contains all units.
            from matplotlib.transforms import Bbox
            fig.canvas.draw()
            boxes = [ax.get_tightbbox(fig.canvas.get_renderer()) for ax in preview_axes]
            box = Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted()).expanded(1.02, 1.03)
            fig.savefig(preview_png, dpi=140, bbox_inches=box)
        elif preview_png:
            fig.savefig(preview_png, dpi=100)
        plt.close(fig)
