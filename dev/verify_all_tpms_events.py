#!/usr/bin/env python3
"""Run the validated stationary-event search on all retained unresolved cells.

The base continuation partition is audited first. Every successful local
certificate is freshly replayed with both backends against the original field
and required to lie inside the same queried parameter cell. This is not a
complete isotopy classification, not a proof every event changes topology,
and not a request to alter formulas, cavity processing, figures or manuscript.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from time import monotonic

from audit_refined_tpms_records import audit
from certify_tpms_event_candidates import event_inside_cell, problem_functions
from certify_tpms_continuation import write_json
from knotted_graph.core.critical_events import verify_critical_point


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('use a new output directory')
    args.out.mkdir(parents=True)
    base = audit(args.root, args.out / 'partition', 20)
    leaves = json.loads((args.out / 'partition/adaptive_leaves.json').read_text())
    expected = {r['id']: r for r in leaves if r['status'] != 'certified'}
    outcomes = {}
    started = monotonic()
    for shard in range(20):
        source = args.root / f'refined-tpms-shard-{shard}'
        directory = args.out / f'event-shard-{shard}'
        subprocess.run([sys.executable, str(Path(__file__).with_name('certify_tpms_event_candidates.py')),
                        '--input', str(source), '--out', str(directory), '--seconds', '600'], check=True)
        results = json.loads((directory / 'event_records.json').read_text())
        for result in results:
            key = result['id']
            if key in outcomes or key not in expected:
                raise ValueError('duplicate event result or query outside unresolved partition')
            target = expected[key]
            if result['status'] == 'certified_stationary_event_in_cell':
                lam = float.fromhex(result['lambda_hex'])
                problem, *_ = problem_functions(target['family'], lam, result['kind'])
                certificate = result['certificate']
                for backend in ('primary', 'rational'):
                    checked = verify_critical_point(problem, certificate, backend=backend)
                    if not checked['valid'] or not event_inside_cell(target, lam, checked):
                        raise ValueError('event existence or exact parameter-cell inclusion failed replay')
                result['collector_replayed_both_backends'] = True
            result['parameter_area_in_root_cell_units'] = 4 ** (-target['depth'])
            result['input_family'] = target['family']
            outcomes[key] = result
        print('ALL_EVENT_SHARD ' + json.dumps({'shard': shard, 'processed': len(outcomes),
              'expected': len(expected), 'statuses': dict(Counter(r['status'] for r in outcomes.values()))}), flush=True)
    if outcomes.keys() != expected.keys():
        raise ValueError('not all unresolved leaves have an event-search outcome')
    summary = {'base_partition': base, 'unresolved_leaf_cells': len(expected),
               'all_unresolved_cells_accounted': True,
               'status_counts': dict(Counter(r['status'] for r in outcomes.values())),
               'by_family': {}, 'seconds': monotonic() - started,
               'all_isotopy_classes_distinguished': False, 'topology_changes_certified': False,
               'one_event_may_be_represented_in_multiple_parameter_cells': True,
               'new_result_is_stationary_event_existence_not_phase_count': True}
    for family in base['by_family']:
        subset = [r for r in outcomes.values() if r['input_family'] == family]
        summary['by_family'][family] = {
            'statuses': dict(Counter(r['status'] for r in subset)),
            'event_cell_area_fraction': sum(r['parameter_area_in_root_cell_units'] for r in subset
                                            if r['status'] == 'certified_stationary_event_in_cell') / 400,
            'still_unknown_cell_area_fraction': sum(r['parameter_area_in_root_cell_units'] for r in subset
                                            if r['status'] != 'certified_stationary_event_in_cell') / 400,
            'certified_regular_area_fraction': base['by_family'][family]['certified_parameter_area_fraction'],
        }
    write_json(args.out / 'all_event_records.json', sorted(outcomes.values(), key=lambda r: r['id']))
    write_json(args.out / 'all_event_summary.json', summary)
    write_json(args.out / 'result_hashes.json', {p.name: sha256(p.read_bytes()).hexdigest()
                                               for p in args.out.glob('*.json')})
    print('ALL_TPMS_EVENT_SUMMARY ' + json.dumps(summary, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
