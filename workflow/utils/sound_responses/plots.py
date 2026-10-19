"""Diagnostic plotting from saved Bayesian-block results; no refitting."""
import json
from pathlib import Path

import h5py
import numpy as np


def plot_results(input_h5, output_pdf):
    """PSTHs with selected boundaries; plotting never changes the timing estimates."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(input_h5, 'r') as source, PdfPages(output_pdf) as pdf:
        settings = json.loads(source.attrs['effective_settings'])
        rows = [(cond, name) for cond in source['conditions'] for name in source[f'conditions/{cond}/units']]
        for start in range(0, max(1, len(rows)), 6):
            fig, axes = plt.subplots(3, 2, figsize=(11, 10), squeeze=False)
            for ax in axes.flat:
                ax.set_visible(False)
            for ax, (cond, name) in zip(axes.flat, rows[start:start + 6]):
                ax.set_visible(True)
                group = source[f'conditions/{cond}']
                unit = group[f'units/{name}']
                attrs = unit.attrs
                edges = unit['psth_edges_s'][:] * 1000
                ax.stairs(unit['psth_rate_hz'][:], edges, color='0.25')
                duration = group.attrs['duration_s']
                if np.isfinite(duration):
                    ax.axvspan(0, duration * 1000, color='0.6', alpha=0.15, label='Sound')
                for key, label, color in [('onset_latency_s', 'Response onset', 'tab:blue'),
                                          ('offset_latency_s', 'Response offset', 'tab:orange')]:
                    if np.isfinite(attrs[key]):
                        ax.axvline(attrs[key] * 1000, color=color, label=label)
                ax.set(xlim=(0, settings['window_s'] * 1000), xlabel='Time from sound onset (ms)',
                       ylabel='Firing rate (Hz)', title=f"{name} | {group.attrs['label']}\n"
                       f"{attrs['status']}; mod={attrs['is_modulated']}; n={attrs['n_events']}")
                ax.text(0.98, 0.97, f"onset: {attrs['onset_status']}\noffset: {attrs['offset_status']}",
                        transform=ax.transAxes, ha='right', va='top', fontsize=7)
                handles, labels = ax.get_legend_handles_labels()
                if handles:
                    ax.legend(handles, labels, fontsize=7, loc='upper left')
            if not rows:
                axes.flat[0].set_visible(True)
                axes.flat[0].text(0.05, 0.5, 'No positive sound conditions or selected units.')
                axes.flat[0].axis('off')
            fig.suptitle('Bayesian-block modulation and timing | -1 unavailable, 0 no change, 1 change')
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)
