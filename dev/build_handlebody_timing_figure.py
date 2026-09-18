#!/usr/bin/env python3
"""Render historical timing records with explicit sample counts and matched engines."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, MaxNLocator
import numpy as np
import pandas as pd


def true_flag(series):
    return series.astype(str).str.lower().eq('true')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Output filename prefix, without extension')
    args = parser.parse_args()
    data = pd.read_csv(args.input)
    paired = (data['topoly_status'].eq('ok')
              & true_flag(data['topoly_matches_knottedgraph_up_to_unit'])
              & true_flag(data['recomputed_matches_stored_recovered_yamada'])
              & np.isfinite(data['topoly_yamada_seconds'])
              & np.isfinite(data['knottedgraph_yamada_seconds'])
              & data['topoly_yamada_seconds'].gt(0)
              & data['knottedgraph_yamada_seconds'].gt(0))
    fields = [
        ('Skeletonization', 'skeletonization_seconds'),
        ('Graph extraction', 'graph_extraction_seconds'),
        ('Edge contraction', 'contract_short_edges_seconds'),
        ('Leaf removal', 'remove_leaf_nodes_seconds'),
        ('Edge simplification', 'simplify_edges_seconds'),
        ('Edge smoothing', 'smooth_edges_seconds'),
        ('Projection attempts', 'projection_selection_seconds'),
        ('KnottedGraph Yamada', 'knottedgraph_yamada_seconds'),
        ('Topoly Yamada', 'topoly_yamada_seconds'),
    ]
    plt.rcParams.update({'font.family': 'DejaVu Serif', 'font.size': 9,
                         'axes.labelsize': 9, 'axes.titlesize': 10, 'xtick.labelsize': 8,
                         'ytick.labelsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'axes.linewidth': .65})
    fig, axes = plt.subplots(3, 3, figsize=(8.1, 7.2), constrained_layout=True)
    colors = ['#267E83', '#337C9B', '#4B78AE', '#7771AC', '#7A69A5', '#916EA7', '#A76591', '#BA6370', '#C7824E']
    summary, plotted = [], []
    for index, (ax, (title, field), color) in enumerate(zip(axes.flat, fields, colors)):
        valid = np.isfinite(data[field]) & data[field].gt(0)
        if index >= 7:
            valid &= paired
        values = data.loc[valid, field].to_numpy()
        mean = float(np.mean(values)); p95 = float(np.quantile(values, .95))
        bounds = np.geomspace(values.min() * (1 - 1e-12), values.max() * (1 + 1e-12), 41)
        hist, _, _ = ax.hist(values, bins=bounds, color=color, edgecolor='white', linewidth=.25)
        assert int(hist.sum()) == len(values), 'Histogram omitted a timing record'
        ax.set_xscale('log')
        ax.axvline(mean, color='#202020', linewidth=1.05, label='Mean')
        ax.axvline(p95, color='#202020', linewidth=1.05, linestyle=(0, (4, 3)), label='95th percentile')
        ax.set_title(f'({chr(97 + index)}) {title}', loc='left', pad=25, fontweight='bold')
        ax.text(.98, .96, f'n = {len(values):,}', transform=ax.transAxes, va='top', ha='right', fontsize=8, bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .95, 'pad': 1})
        ax.text(0, 1.045, f'mean {mean:.3g} s; P95 {p95:.3g} s', transform=ax.transAxes, ha='left', va='bottom', fontsize=7.5)
        ax.set_xlabel('Time (s)'); ax.set_ylabel('Count')
        ax.xaxis.set_major_locator(LogLocator(base=10, numticks=4))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
        ax.set_axisbelow(True); ax.grid(axis='y', color='#DADDE0', linewidth=.45)
        ax.spines[['top', 'right']].set_visible(False)
        ax.margins(x=.035, y=.16)
        cohort = 'same completed and matched diagrams' if index >= 7 else 'all available historical records'
        summary.append({'panel': chr(97 + index), 'stage': field, 'cohort': cohort, 'n': int(len(values)),
                        'mean_seconds': mean, 'p95_seconds': p95, 'minimum_seconds': float(values.min()), 'maximum_seconds': float(values.max())})
        plotted.extend({'panel': chr(97 + index), 'index': int(row['index']), 'case': row['case'], 'time_seconds': float(row[field]), 'cohort': cohort} for _, row in data.loc[valid].iterrows())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output.with_suffix('.pdf'), metadata={'Title': 'Historical handlebody stage timings with matched Yamada comparisons', 'Creator': 'build_timing_figure.py'})
    plt.close(fig)
    pd.DataFrame(summary).to_csv(args.output.with_name(args.output.name + '_summary.csv'), index=False)
    pd.DataFrame(plotted).to_csv(args.output.with_name(args.output.name + '_plotted_records.csv'), index=False)
    shutil.copyfile(args.input, args.output.with_name(args.output.name + '_source.csv'))
    metadata = {'source_path': str(args.input.resolve()), 'source_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
                'source_rows': len(data), 'paired_engine_rows': int(paired.sum()),
                'topoly_status_counts': data['topoly_status'].value_counts(dropna=False).to_dict(),
                'selection': 'Panels a-g use all finite positive recorded stage timings. Panels h-i both require Topoly status ok, recorded agreement up to a unit, recorded KG/stored-recovered agreement, and finite positive durations for both engines.',
                'records_are_historical': True, 'current_projection_code_rerun': False,
                'topoly_completed_above_reported_timeout': int((paired & data['topoly_yamada_seconds'].gt(data['topoly_timeout_s'])).sum()),
                'timeout_rows_with_stale_match_flag': int((data['topoly_status'].eq('timeout') & true_flag(data['topoly_matches_knottedgraph_up_to_unit'])).sum()),
                'panels': summary}
    args.output.with_name(args.output.name + '_provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
