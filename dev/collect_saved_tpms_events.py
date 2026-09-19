#!/usr/bin/env python3
"""Recover complete event results after the historical empty-shard write bug.

Read retained artifacts without modifying them. Re-audit the original adaptive
partition, require each shard's exact planned IDs and input hashes, and freshly
replay every local event certificate against its trusted defining field with
both interval backends. Missing data are accepted only for an explicitly empty
historical shard with a matching zero-row plan and summary.

Known exact P-to-D stationary branches can supply additional positive event
witnesses, using rational arithmetic for full parameter-cell membership. These
are event-containing rectangles, not proved topology-changing transitions or
counts of distinct isotopy classes. No figures or manuscript are produced.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from time import monotonic

from audit_refined_tpms_records import audit
from certify_tpms_event_candidates import event_inside_cell, problem_functions
from certify_tpms_continuation import write_json
from verify_exact_tpms_branches import verify_identities
from knotted_graph.core.critical_events import verify_critical_point

POSITIVE = 'certified_stationary_event_in_cell'
EXACT = 'exact_stationary_branch_event_in_cell'


def read_saved_rows(directory, expected_ids):
    directory = Path(directory)
    plan = json.loads((directory / 'plan.json').read_text())
    summary = json.loads((directory / 'event_summary.json').read_text())
    if (plan['selected_ids'] != expected_ids or summary['selected_cells'] != len(expected_ids)
            or summary['recorded_cells'] != len(expected_ids)):
        raise ValueError('event selection or completion count mismatch')
    path = directory / 'event_records.json'
    legacy_empty = not path.exists()
    if legacy_empty:
        if expected_ids or summary['status_counts'] != {}:
            raise ValueError('missing records in a nonempty or inconsistent event shard')
        rows = []
    else:
        rows = json.loads(path.read_text())
    if not isinstance(rows, list) or [r['id'] for r in rows] != expected_ids:
        raise ValueError('event records do not exactly match the planned IDs')
    if dict(Counter(r['status'] for r in rows)) != summary['status_counts']:
        raise ValueError('event status summary does not match retained records')
    return rows, legacy_empty, plan


def exact_branch_witness(row):
    """Algebraic membership only; caller must first verify field identities."""
    if row['family'] != 'schwarz_p_to_diamond':
        return None
    left, right = map(Fraction, row['lambda_bounds'])
    bottom, top = map(Fraction, row['c_bounds'])
    if left > right or bottom > top:
        raise ValueError('unordered parameter cell')
    candidates = [
        ('1-2*lambda', (1-top)/2, (1-bottom)/2, ('pi', '0', '0'), -2, 1),
        ('2*lambda-1', (bottom+1)/2, (top+1)/2, ('pi', 'pi', '0'), 2, -1),
    ]
    for branch, low, high, point, slope, offset in candidates:
        low, high = max(left, low), min(right, high)
        if low > high:
            continue
        lam = (low + high) / 2
        level = slope * lam + offset
        if not left <= lam <= right or not bottom <= level <= top:
            raise AssertionError('exact event membership arithmetic failed')
        return {'critical_branch': branch, 'point': list(point),
                'lambda_fraction': [lam.numerator, lam.denominator],
                'c_fraction': [level.numerator, level.denominator],
                'membership_arithmetic': 'exact_rationals_from_binary_parameter_endpoints',
                'isolated_point_asserted': False, 'topology_change_certified': False}
    return None


def collect(saved, continuation, out):
    saved, continuation, out = Path(saved), Path(continuation), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    base = audit(continuation, out / 'partition', 20)
    leaf_list = json.loads((out / 'partition/adaptive_leaves.json').read_text())
    targets = {r['id']: r for r in leaf_list if r['status'] != 'certified'}
    identities = verify_identities()
    if not identities['all_symbolic_residuals_zero'] or not identities['rational_ball_inclusion_proved']:
        raise RuntimeError('exact stationary-branch identities were not verified')
    write_json(out / 'exact_branch_identities.json', identities)
    results, input_hashes, empty_legacy = {}, {}, []
    source_commits = set()
    started = monotonic()
    for shard in range(20):
        input_dir = continuation / f'refined-tpms-shard-{shard}'
        directory = saved / f'event-shard-{shard}'
        source_rows = json.loads((input_dir / 'records.json').read_text())
        source_plan = json.loads((input_dir / 'plan.json').read_text())
        expected = [r['id'] for r in source_rows if not r['children'] and r['status'] != 'certified']
        rows, legacy_empty, plan = read_saved_rows(directory, expected)
        if (plan['input_records_sha256'] != sha256((input_dir / 'records.json').read_bytes()).hexdigest()
                or plan['source_continuation_plan'] != source_plan):
            raise ValueError('event search refers to different adaptive source records')
        source_commits.add(plan['commit'])
        if legacy_empty:
            empty_legacy.append(shard)
        for name in ('plan.json', 'event_summary.json', 'event_records.json'):
            path = directory / name
            if path.exists():
                input_hashes[str(path.relative_to(saved))] = sha256(path.read_bytes()).hexdigest()
        for row in rows:
            key = row['id']
            if key in results or key not in targets:
                raise ValueError('duplicate or unexpected event record')
            target = targets[key]
            result = dict(row)
            result.update(family=target['family'], lambda_bounds=target['lambda_bounds'],
                          c_bounds=target['c_bounds'], original_search_status=row['status'],
                          parameter_depth=target['depth'], cell_weight=4 ** (-target['depth']))
            if row['status'] == POSITIVE:
                if any(row.get(k) != target[k] for k in ('family', 'lambda_bounds', 'c_bounds')):
                    raise ValueError('certified event query metadata changed')
                lam = float.fromhex(row['lambda_hex'])
                problem, *_ = problem_functions(target['family'], lam, row['kind'])
                replays = {}
                for backend in ('primary', 'rational'):
                    checked = verify_critical_point(problem, row['certificate'], backend=backend)
                    if not checked['valid'] or not event_inside_cell(target, lam, checked):
                        raise ValueError('saved event failed fresh existence or query-membership verification')
                    replays[backend] = checked
                result['fresh_collector_replays'] = replays
            elif row['status'] == 'unknown':
                witness = exact_branch_witness(target)
                if witness is not None:
                    result['status'] = EXACT
                    result['exact_branch_witness'] = witness
            elif row['status'] != 'not_executed_time_budget':
                raise ValueError('unexpected event-search status')
            result['topology_change_certified'] = False
            result['all_critical_points_enumerated'] = False
            results[key] = result
        write_json(out / 'event_records.json', sorted(results.values(), key=lambda r: r['id']))
        print('REPLAYED_EVENT_SHARD ' + json.dumps({'shard': shard, 'accounted': len(results),
              'expected': len(targets), 'statuses': dict(Counter(r['status'] for r in results.values()))}), flush=True)
    if results.keys() != targets.keys() or len(source_commits) != 1:
        raise ValueError('incomplete event coverage or mixed search source commits')
    statuses = Counter(r['status'] for r in results.values())
    original_statuses = Counter(r['original_search_status'] for r in results.values())
    summary = {'event_search_run_id': 35414418156, 'retained_artifact_id': 10574974757,
               'event_search_source_commit': source_commits.pop(),
               'collector_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
               'original_unresolved_leaf_count': len(targets), 'accounted_leaf_count': len(results),
               'original_search_outcomes': dict(original_statuses), 'replayed_and_exact_outcomes': dict(statuses),
               'every_retained_local_positive_replayed_with_both_backends': True,
               'exact_branch_identities_reverified': True, 'legacy_empty_shards': empty_legacy,
               'original_artifacts_modified': False, 'by_family': {}, 'seconds': monotonic()-started,
               'distinct_stationary_points_counted': False, 'source_topology_changes_proved': False,
               'all_source_isotopy_classes_distinguished': False, 'manuscript_rebuilt': False,
               'figures_regenerated': False, 'analytic_to_voxel_correspondence_proved': False}
    for family, regular in base['by_family'].items():
        subset = [r for r in results.values() if r['family'] == family]
        event_area = sum(r['cell_weight'] for r in subset if r['status'] in (POSITIVE, EXACT)) / 400
        unknown_area = sum(r['cell_weight'] for r in subset if r['status'] not in (POSITIVE, EXACT)) / 400
        regular_area = regular['certified_parameter_area_fraction']
        if regular_area + event_area + unknown_area != 1:
            raise ValueError('regular/event/unknown area accounting does not partition the domain')
        summary['by_family'][family] = {
            'outcomes': dict(Counter(r['status'] for r in subset)),
            'certified_regular_area_fraction': regular_area,
            'event_containing_cell_area_fraction': event_area,
            'still_unknown_cell_area_fraction': unknown_area,
            'event_cell_area_is_not_area_of_the_critical_set': True,
        }
    write_json(out / 'event_summary.json', summary)
    write_json(out / 'input_artifact_file_hashes.json', input_hashes)
    write_json(out / 'result_file_hashes.json', {p.name: sha256(p.read_bytes()).hexdigest()
                                               for p in out.glob('*.json')})
    print('COMPLETE_EVENT_RECOVERY ' + json.dumps(summary, sort_keys=True), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--saved', type=Path, required=True)
    parser.add_argument('--continuation', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    collect(args.saved, args.continuation, args.out)
