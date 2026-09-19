#!/usr/bin/env python3
"""Data-only summary of complete map records; no manuscript or figure edits.

All counts derive from the collected exact grids. Polynomial signatures are
explicitly distinguished from continuum isotopy classes. Source arrays, grid
identities, data hashes and every unavailable outcome are retained in CSV/JSON.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path

from run_full_resolution_maps import dump, plan


def analyze(maps, out, adaptive=None):
    maps, out = Path(maps), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((maps / 'full_map_manifest.json').read_text())
    results = {'experiment': 'original full grids at 120 cubed', 'modes': {},
               'input_manifest_sha256': sha256((maps / 'full_map_manifest.json').read_bytes()).hexdigest(),
               'source_isotopy_classification_proved': False,
               'source_to_voxel_correspondence_proved': False,
               'figures_regenerated': False, 'manuscript_rebuilt': False}
    text = ['# Complete-grid numerical execution', '',
            'These are newly executed finite-domain voxel-complex and selected-spine results. '
            'A polynomial signature is not an equivalence certificate for the analytic source solid.', '']
    for mode, expected_count in (('hamiltonian', 7140), ('tpms', 1323)):
        path = maps / f'{mode}_records.json'
        rows = json.loads(path.read_text())
        expected = {r['id']: r for r in plan(mode)}
        observed = {r['id']: r for r in rows}
        if len(rows) != expected_count or len(observed) != expected_count or observed.keys() != expected.keys():
            raise ValueError('summary input is not the exact complete requested grid')
        if any(any(row.get(k) != value for k, value in expected[row['id']].items()) for row in rows):
            raise ValueError('summary input changed parameter coordinates')
        available = [r for r in rows if r['status'] == 'evaluated' and r.get('yamada') is not None]
        full = [r for r in rows if r.get('reconstruction', {}).get('kind') == 'graph_retract']
        family_summary = {}
        for family in dict.fromkeys(r['family'] for r in rows):
            subset = [r for r in rows if r['family'] == family]
            successful = [r for r in subset if r['status'] == 'evaluated' and r.get('yamada') is not None]
            changes = 0
            lookup = {(r['level_index'], r['lambda_index']): r for r in subset}
            for row in subset:
                right = lookup.get((row['level_index'], row['lambda_index'] + 1))
                if right is None or row not in successful or right['status'] != 'evaluated':
                    continue
                left_scope = row.get('evaluation', {})
                right_scope = right.get('evaluation', {})
                compatible = all(left_scope.get(k) == right_scope.get(k)
                                 for k in ('evaluation_kind', 'normalization', 'is_subcubic'))
                if compatible and row['signature'] != right['signature']:
                    changes += 1
            family_summary[family] = {
                'cells': len(subset), 'statuses': dict(Counter(r['status'] for r in subset)),
                'lambda_samples': len({r['lambda'] for r in subset}),
                'level_samples': len({r['level'] for r in subset}),
                'level_range': [min(r['level'] for r in subset), max(r['level'] for r in subset)],
                'boundary_touching_source_masks': sum(r.get('volume', {}).get('touches_boundary', False) for r in subset),
                'distinct_successful_polynomial_values': len({r['yamada'] for r in successful}),
                'distinct_recorded_successful_signatures': len({r['signature'] for r in successful}),
                'adjacent_compatible_spine_signature_changes': changes,
                'these_counts_are_not_source_topology_classes_or_transitions': True,
            }
        mode_summary = {
            'source_commit': manifest[mode]['commit'], 'requested_cells': expected_count,
            'recorded_cells': len(rows), 'dimension': 120,
            'statuses': dict(Counter(r['status'] for r in rows)),
            'successful_spatial_or_isolated_vertex_cells': len(available),
            'graph_retract_cells': len(full),
            'nontrivial_subcubic_cells_with_two_projection_checks': sum(len(r.get('projection_checks', [])) >= 2 for r in available),
            'distinct_polynomial_values_not_isotopy_classes': len({r['yamada'] for r in available}),
            'data_sha256': sha256(path.read_bytes()).hexdigest(), 'by_family': family_summary,
        }
        results['modes'][mode] = mode_summary
        print('FULL_GRID_RESULT ' + json.dumps({mode: mode_summary}, sort_keys=True), flush=True)
        dump(out / f'{mode}_summary.json', mode_summary)
        (out / path.name).write_bytes(path.read_bytes())
        columns = ['id', 'family', 'lambda', 'level', 'lambda_hex', 'level_hex', 'dimension', 'status',
                   'yamada', 'source_mask_sha256', 'artifact_relative_path', 'error']
        with (out / f'{mode}_all_cells.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows({k: r.get(k) for k in columns} for r in rows)
        text.extend([f'## {mode.title()}', '', f'Requested and recorded cells: **{expected_count} / {len(rows)}** at **120³**.', '',
                     '| Family | Cells | Successful polynomial evaluations | Distinct values |',
                     '|---|---:|---:|---:|'])
        for family, result in family_summary.items():
            text.append(f"| {family} | {result['cells']} | {result['statuses'].get('evaluated', 0)} | {result['distinct_successful_polynomial_values']} |")
        text.extend(['', 'Outcome counts: `' + json.dumps(mode_summary['statuses'], sort_keys=True) + '`.', '',
                     'The last column counts polynomial values, not proved source-solid classes. '
                     'All fixed-diagram, unavailable and existing-cavity-route outcomes remain separate.', ''])
        if mode == 'tpms':
            selected = [r for r in rows if r['family'] == 'gyroid_to_diamond' and r['level_index'] == 7
                        and r['lambda_index'] in (0, 1, 5, 11, 17)]
            results['original_five_tpms_representatives'] = [
                {k: r.get(k) for k in ('id', 'lambda', 'level', 'status', 'graph', 'yamada',
                                       'source_mask_sha256', 'artifact_relative_path')}
                for r in selected]
    if adaptive is not None:
        adaptive = Path(adaptive)
        results['adaptive_continuation'] = json.loads((adaptive / 'adaptive_summary.json').read_text())
        text.extend(['## Analytic TPMS continuation', '',
                     'Accepted rectangles have the recorded dual interval replay. '
                     'Their connected cover components are positive continuation groups, not counts of distinct isotopy classes.', '',
                     '| Family | Certified area fraction | Unresolved area fraction |',
                     '|---|---:|---:|'])
        for family, result in results['adaptive_continuation']['by_family'].items():
            text.append(f"| {family} | {result['certified_parameter_area_fraction']:.6f} | {result['unresolved_parameter_area_fraction']:.6f} |")
    text.extend(['', '## Scope', '',
                 'The full numerical runs do not establish analytic-source-to-voxel correspondence, '
                 'a complete source-solid isotopy classification, regenerated main figures, or a rebuilt manuscript. '
                 'The all-word formulas and existing cavity implementation were not changed.', ''])
    dump(out / 'numerical_execution_summary.json', results)
    (out / 'NUMERICAL_RESULTS.md').write_text('\n'.join(text))
    hashes = {p.name: sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    dump(out / 'data_manifest.json', hashes)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maps', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--adaptive', type=Path)
    args = parser.parse_args()
    analyze(args.maps, args.out, args.adaptive)
