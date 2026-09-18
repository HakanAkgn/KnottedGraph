#!/usr/bin/env python3
"""Plot witnessed continuation components, keeping unknown rectangles visible."""
from __future__ import annotations
import argparse
from collections import deque
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

FAMILIES = ('schwarz_p_to_diamond', 'gyroid_to_schwarz_p', 'gyroid_to_diamond')
TITLES = ('Schwarz-P to Diamond', 'Gyroid to Schwarz-P', 'Gyroid to Diamond')


def continuation_components(mask):
    """Closed certified rectangles may share an edge or a certified corner."""
    labels = np.zeros(mask.shape, dtype=int)
    sizes = []
    for start in zip(*np.nonzero(mask)):
        if labels[start]:
            continue
        label = len(sizes)+1
        labels[start] = label
        queue = deque([start])
        size = 0
        while queue:
            i, j = queue.popleft()
            size += 1
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    a, b = i+di, j+dj
                    if (0 <= a < mask.shape[0] and 0 <= b < mask.shape[1]
                            and mask[a, b] and labels[a, b] == 0):
                        labels[a, b] = label
                        queue.append((a, b))
        sizes.append(size)
    return labels, sizes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output exists; use a new directory')
    args.out.mkdir(parents=True)
    plan = json.loads((args.data/'plan.json').read_text())
    rows = json.loads((args.data/'records.json').read_text())
    n = plan['samples_per_axis']
    if plan['mode'] != 'grid' or len(rows) != plan['planned_cases']:
        raise ValueError('a complete grid attempt ledger is required')
    xs, ys = np.linspace(0., 1., n), np.linspace(0., .3, n)
    metadata = {'records_sha256': hashlib.sha256((args.data/'records.json').read_bytes()).hexdigest(),
                'quantity': 'connected components of witnessed analytic continuation',
                'unknown_is_inequivalent': False, 'different_colors_are_inequivalent': False,
                'cross_family_color_comparison': False, 'filtering': False,
                'closed_rectangle_adjacency': 'shared edge or corner', 'families': {}}
    for family, title in zip(FAMILIES, TITLES):
        grid = np.zeros((n-1, n-1), dtype=bool)
        seen = set()
        for row in rows:
            if row['family'] != family:
                continue
            ij = row['c_index'], row['lambda_index']
            if ij in seen:
                raise ValueError('duplicate parameter rectangle')
            seen.add(ij)
            if row['status'] == 'certified':
                if not row['replay_valid']:
                    raise ValueError('unreplayed certificate cannot be plotted as certified')
                grid[ij] = True
            elif row['status'] != 'unknown':
                raise ValueError('unexpected result status')
        if len(seen) != (n-1)**2:
            raise ValueError('missing parameter rectangles')
        labels, sizes = continuation_components(grid)
        fig, ax = plt.subplots(figsize=(5.1, 3.95))
        fig.subplots_adjust(left=.15, right=.97, bottom=.18, top=.87)
        ax.pcolormesh(xs, ys, np.ma.masked_where(~grid, labels),
                      vmin=.5, vmax=max(1.5, len(sizes)+.5), shading='flat', rasterized=False)
        for i, j in zip(*np.nonzero(~grid)):
            ax.add_patch(Rectangle((xs[j], ys[i]), xs[j+1]-xs[j], ys[i+1]-ys[i],
                                   fill=False, hatch='////', linewidth=0))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, .3)
        ax.set_xticks([0, .25, .5, .75, 1])
        ax.set_yticks([0, .1, .2, .3])
        ax.set_xlabel(r'Interpolation $\lambda$', fontsize=13)
        ax.set_ylabel(r'Level offset $c$', fontsize=13)
        ax.set_title(title, fontsize=15, pad=10)
        ax.tick_params(labelsize=11)
        fig.savefig(args.out/(family+'.pdf'))
        fig.savefig(args.out/(family+'.png'), dpi=180)
        plt.close(fig)
        metadata['families'][family] = {'certified_rectangles': int(grid.sum()),
            'unknown_rectangles': int((~grid).sum()), 'continuation_components': len(sizes),
            'component_sizes': sizes, 'labels': labels.tolist(), 'lambda_edges': xs.tolist(),
            'c_edges': ys.tolist()}
    (args.out/'figure_data.json').write_text(json.dumps(metadata, indent=2)+'\n')
    try:
        import fitz
    except ImportError:
        return
    output = fitz.open()
    width, height = 5.1*72, 3.95*72
    page = output.new_page(width=3*width, height=height)
    for i, family in enumerate(FAMILIES):
        with fitz.open(args.out/(family+'.pdf')) as source:
            page.show_pdf_page(fitz.Rect(i*width, 0, (i+1)*width, height), source, 0)
    output.save(args.out/'TPMS_certified_continuation.pdf', deflate=True)
    output[0].get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(args.out/'TPMS_certified_continuation.png')
    output.close()


if __name__ == '__main__':
    main()
