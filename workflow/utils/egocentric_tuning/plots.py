"""Per-unit curves and condition contrasts, with eligibility shown explicitly."""
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from .core import PHASES


def plot_tuning(input_h5, output_pdf, output_preview):
    for p in (output_pdf, output_preview):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(input_h5, 'r') as f, PdfPages(output_pdf) as pdf:
        if not f.attrs.get('complete', False):
            raise ValueError('Cannot plot incomplete tuning artifact')
        conditions = json.loads(f.attrs['conditions'])
        pairs = json.loads(f.attrs['comparisons'])
        angles = f['angular_centers_rad'][:]
        close = lambda v: np.r_[v, v[0]]
        colors = plt.get_cmap('tab10')
        page = 0
        for unit in json.loads(f.attrs['unit_names']):
            for ref in json.loads(f.attrs['references']):
                fig = plt.figure(figsize=(13, 8))
                grid = fig.add_gridspec(2, 3, top=.84, bottom=.12, hspace=.65)
                for col, phase in enumerate(PHASES):
                    ax = fig.add_subplot(grid[0, col], projection='polar')
                    notes = []
                    for i, condition in enumerate(conditions):
                        cg = f['conditions'][condition]
                        pg = cg[f'units/{unit}/{phase}']
                        g = pg[ref]
                        label = str(cg.attrs.get('label', condition))
                        status = str(g.attrs['status'])
                        if 'rate_hz' in g:
                            ax.plot(close(angles), close(g['rate_hz'][:]), color=colors(i % 10), label=label)
                            ci = g['rate_ci95_hz'][:]
                            ax.fill_between(close(angles), close(ci[0]), close(ci[1]),
                                            color=colors(i % 10), alpha=.13)
                            label_info = (g.attrs['classification'] if status == 'eligible' else status)
                            notes.append(f"{label}: {label_info}, n={g.attrs['n_events']}, "
                                         f"spikes={g.attrs['n_spikes']}\n"
                                         f"VS={g.attrs['vector_strength']:.2f}; "
                                         f"VS2={g.attrs['orientation_strength']:.2f}")
                        else:
                            notes.append(f'{label}: {pg.attrs["timing_status"]}')
                    ax.set_title(phase.capitalize()+' (Hz)', pad=15)
                    ax.set_ylim(bottom=0)
                    if col == 0 and ax.lines:
                        ax.legend(loc='upper left', bbox_to_anchor=(-.55, 1.2), fontsize=7)
                    ax.text(.5, -.24, '\n'.join(notes), transform=ax.transAxes,
                            fontsize=7, ha='center', va='top')
                    bx = fig.add_subplot(grid[1, col])
                    if pairs:
                        for i, (a, b) in enumerate(pairs):
                            g = f[f'comparisons/{i}/{unit}/{ref}/{phase}']
                            label = f'{f["conditions"][b].attrs.get("label", b)} − {f["conditions"][a].attrs.get("label", a)}'
                            if 'rate_difference_hz' not in g:
                                bx.text(.5, .5, 'Timing unavailable', ha='center', transform=bx.transAxes)
                                continue
                            bx.plot(np.rad2deg(angles), g['rate_difference_hz'][:], 'o-', ms=3,
                                    color=colors(i % 10), label=label)
                            ci = g['rate_difference_ci95_hz'][:]
                            bx.fill_between(np.rad2deg(angles), ci[0], ci[1], color=colors(i % 10), alpha=.15)
                        bx.axhline(0, color='grey', lw=.8)
                        if bx.lines and bx.get_legend_handles_labels()[0]:
                            bx.legend(fontsize=7)
                        bx.set_ylabel('Rate difference (Hz)')
                    else:
                        for i, condition in enumerate(conditions):
                            g = f[f'conditions/{condition}/units/{unit}/{phase}/{ref}']
                            if 'presentations' in g:
                                bx.plot(np.rad2deg(angles), g['presentations'][:], 'o-', color=colors(i % 10))
                        bx.set_ylabel('Presentations')
                    bx.set_xlabel('Egocentric bearing (degrees)')
                    bx.set_xticks([0, 90, 180, 270])
                    bx.grid(alpha=.15)
                timing = json.loads(f.attrs['config'])['timing_condition'] or 'per condition'
                fig.suptitle(f'{f.attrs["session"]} | unit {unit} | {ref}\n'
                             f'Amaro-style egocentric tuning; timing: {timing}', fontsize=13)
                fig.text(.5, .035, 'Shading: pointwise 95% presentation-bootstrap intervals (no multiple-comparison correction).\n'
                         'Late includes pre-onset time. Classification uses coverage, spike-count and VS criteria; '
                         'condition differences are descriptive.', ha='center', fontsize=8)
                pdf.savefig(fig)
                if page == 0:
                    fig.savefig(output_preview, dpi=140)
                plt.close(fig)
                page += 1
