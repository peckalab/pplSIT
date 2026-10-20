"""Descriptive angle coverage; counts are not neural tuning curves."""
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .core import wrap


def plot_geometry(path, pdf_path, preview_path):
    Path(pdf_path).parent.mkdir(parents=True, exist_ok=True)
    Path(preview_path).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, 'r') as f, PdfPages(pdf_path) as pdf:
        if not f.attrs.get('complete', False):
            raise ValueError('Geometry file is incomplete')
        refs = json.loads(f.attrs['references'])
        conditions = json.loads(f.attrs['conditions'])
        edges = np.rad2deg(f['angular_bin_edges_rad'][:])
        # Paginate conditions to keep passive-session legends readable.
        pages = [conditions[i:i+6] for i in range(0, len(conditions), 6)] or [[]]
        for page_index, page in enumerate(pages):
            fig, axes = plt.subplots(len(refs), 2, figsize=(12, 3.1*len(refs)), squeeze=False,
                                     constrained_layout=True)
            for row, ref in enumerate(refs):
                for condition in page:
                    g = f[f"events/{ref}/coverage/{condition['sound_id']}"]
                    counts = g['counts'][:]
                    n = int(g.attrs['n_valid'])
                    label = f"{condition['name']} ({n}/{g.attrs['n_events']})"
                    axes[row, 0].stairs(counts, edges, label=label)
                    if n:
                        axes[row, 1].stairs(counts/n, edges, label=condition['name'])
                for col in range(2):
                    axes[row, col].set(xlim=(-180, 180), xlabel='Egocentric speaker bearing (degrees)',
                                       title=ref, ylabel='Sound events' if col == 0 else 'Fraction within condition')
                    if axes[row, col].get_legend_handles_labels()[0]:
                        axes[row, col].legend(fontsize=8)
            fig.suptitle(f"{f.attrs['session']} — sound-onset geometry\n"
                         "Coverage only; 0° in front; positive counterclockwise; legend = usable / total")
            pdf.savefig(fig)
            if page_index == 0:
                fig.savefig(preview_path, dpi=140)
            plt.close(fig)
        if len(refs) > 1:
            fig, axes = plt.subplots(len(refs)-1, 1, figsize=(10, 3*(len(refs)-1)),
                                     squeeze=False, constrained_layout=True)
            anchor = f[f'timeline/{refs[0]}/heading_rad'][:]
            for ax, ref in zip(axes[:, 0], refs[1:]):
                other = f[f'timeline/{ref}/heading_rad'][:]
                use = np.isfinite(anchor) & np.isfinite(other)
                delta = np.rad2deg(wrap(other[use]-anchor[use]))
                ax.hist(delta, bins=np.linspace(-180, 180, 73), color='0.4')
                ax.set(xlabel=f'{ref} minus {refs[0]} (degrees)', ylabel='Valid timeline samples',
                       title=f'{int(use.sum()):,} jointly valid samples')
            fig.suptitle('Direction-reference comparison; shared position origin')
            pdf.savefig(fig)
            plt.close(fig)
