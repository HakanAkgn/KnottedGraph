#!/usr/bin/env python3
"""Seek certified stationary events in unresolved TPMS parameter rectangles.

Approximate roots are candidates only. For each accepted event, both interval
backends prove a unique stationary point in a spatial box and enclose its
critical value entirely inside the queried level interval at an exact binary
lambda. No missing root, failed search or event is relabelled as inequivalence
or as a proved topology-changing transition. No plots or manuscript edits.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from time import monotonic
import warnings

import numpy as np
from scipy.optimize import root
import sympy as sp

from certify_tpms_continuation import write_json
from knotted_graph.core.field_isotopy import tpms_problem
from knotted_graph.core.critical_events import (
    CriticalPointProblem, certify_critical_point, verify_critical_point,
)


@lru_cache(maxsize=256)
def problem_functions(family, lam, kind):
    source = tpms_problem(family, (lam, lam), (0., 0.))
    coordinates = source.variables[:3]
    field = source.expression.subs({source.variables[3]: sp.Rational(lam), source.variables[4]: 0})
    problem = CriticalPointProblem(field, coordinates, source.radius, kind)
    variables, system = problem.system()
    matrix = sp.Matrix(system)
    jacobian = matrix.jacobian(variables)
    return (problem, sp.lambdify(variables, matrix, 'numpy'),
            sp.lambdify(variables, jacobian, 'numpy'),
            sp.lambdify(coordinates, field, 'numpy'),
            sp.lambdify(coordinates, [sp.diff(field, v) for v in coordinates], 'numpy'))


def event_inside_cell(row, lam, result):
    low, high = map(float.fromhex, result['critical_value_hex'])
    return (row['lambda_bounds'][0] <= lam <= row['lambda_bounds'][1]
            and row['c_bounds'][0] <= low <= high <= row['c_bounds'][1])


def attempt(row):
    spatial = row.get('unresolved_spatial_box')
    if not spatial or 'box_hex' not in spatial:
        return {'id': row['id'], 'status': 'unknown', 'reason': 'no_retained_spatial_seed', 'attempts': []}
    bounds = np.asarray([[float.fromhex(a), float.fromhex(b)] for a, b in spatial['box_hex']])
    seed = bounds[:3].mean(axis=1)
    lambda_candidates = [float(bounds[3].mean()),
                         row['lambda_bounds'][0]/2 + row['lambda_bounds'][1]/2,
                         row['lambda_bounds'][0], row['lambda_bounds'][1]]
    lambdas = list(dict.fromkeys(x for x in lambda_candidates
                                if row['lambda_bounds'][0] <= x <= row['lambda_bounds'][1]))
    attempts = []
    radius = tpms_problem(row['family'], row['lambda_bounds'], row['c_bounds']).radius
    kinds = ('wall', 'bulk') if abs(np.linalg.norm(seed) - radius) < .25 else ('bulk', 'wall')
    for lam in lambdas:
        for kind in kinds:
            problem, function, jacobian, scalar, gradient = problem_functions(row['family'], lam, kind)
            initial = seed.copy()
            if kind == 'wall':
                norm = np.linalg.norm(initial)
                if norm < 1e-12:
                    initial = np.array([radius, 0., 0.])
                else:
                    initial *= radius / norm
                multiplier = float(np.dot(initial, np.asarray(gradient(*initial)).ravel()) / radius**2)
                initial = np.r_[initial, multiplier]
            def fun(value):
                return np.asarray(function(*value), dtype=float).ravel()
            def jac(value):
                return np.asarray(jacobian(*value), dtype=float)
            record = {'lambda_hex': float(lam).hex(), 'kind': kind}
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    numerical = root(fun, initial, jac=jac, method='hybr', options={'xtol': 1e-10, 'maxfev': 120})
                center = np.asarray(numerical.x, dtype=float)
                if not np.isfinite(center).all() or np.max(np.abs(fun(center))) > 1e-7:
                    record['status'] = 'candidate_not_converged'
                    attempts.append(record)
                    continue
                approximate_level = float(scalar(*center[:3]))
                if (not np.isfinite(approximate_level) or approximate_level < row['c_bounds'][0] - 1e-7
                        or approximate_level > row['c_bounds'][1] + 1e-7):
                    record['status'] = 'candidate_outside_level_interval'
                    attempts.append(record)
                    continue
                if kind == 'bulk' and np.linalg.norm(center[:3]) >= radius:
                    record['status'] = 'candidate_outside_ball'
                    attempts.append(record)
                    continue
                record['certificate_attempts'] = []
                for width in (1e-5, 1e-7, 1e-9):
                    box = [(float(x-width), float(x+width)) for x in center]
                    certificate = certify_critical_point(problem, box, center=center)
                    record['certificate_attempts'].append({'width': width, 'status': certificate['status'],
                                                           'reason': certificate['reason']})
                    if certificate['status'] != 'certified':
                        continue
                    primary = verify_critical_point(problem, certificate, backend='primary')
                    rational = verify_critical_point(problem, certificate, backend='rational')
                    if not primary['valid'] or not rational['valid']:
                        raise RuntimeError('local event certificate failed fresh replay')
                    if not event_inside_cell(row, lam, primary) or not event_inside_cell(row, lam, rational):
                        continue
                    record['status'] = 'certified_stationary_event_in_cell'
                    attempts.append(record)
                    return {'id': row['id'], 'family': row['family'], 'lambda_bounds': row['lambda_bounds'],
                            'c_bounds': row['c_bounds'], 'lambda_hex': float(lam).hex(), 'kind': kind,
                            'status': 'certified_stationary_event_in_cell', 'certificate': certificate,
                            'primary_replay': primary, 'rational_replay': rational, 'attempts': attempts,
                            'topology_change_certified': False, 'all_critical_points_enumerated': False}
                record['status'] = 'local_root_or_level_membership_not_certified'
            except (ValueError, FloatingPointError, np.linalg.LinAlgError, OverflowError) as exc:
                record['status'] = 'candidate_error'
                record['error'] = f'{type(exc).__name__}: {exc}'
            attempts.append(record)
    return {'id': row['id'], 'family': row['family'], 'lambda_bounds': row['lambda_bounds'],
            'c_bounds': row['c_bounds'], 'status': 'unknown', 'attempts': attempts,
            'topology_change_certified': False, 'all_critical_points_enumerated': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=int, default=600)
    args = parser.parse_args()
    if args.out.exists() or args.seconds < 1:
        parser.error('new output directory and positive time budget required')
    args.out.mkdir(parents=True)
    rows = json.loads((args.input / 'records.json').read_text())
    selected = [r for r in rows if not r['children'] and r['status'] != 'certified']
    write_json(args.out / 'plan.json', {
        'input_records_sha256': sha256((args.input / 'records.json').read_bytes()).hexdigest(),
        'selection': 'every unresolved leaf in supplied complete shard', 'selected_ids': [r['id'] for r in selected],
        'seconds_budget': args.seconds, 'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'driver_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_continuation_plan': json.loads((args.input / 'plan.json').read_text()),
    })
    results = []
    started = monotonic()
    for number, row in enumerate(selected):
        if monotonic() - started >= args.seconds:
            result = {'id': row['id'], 'status': 'not_executed_time_budget'}
        else:
            archive = args.input / row['certificate']
            raw = archive.read_bytes()
            if sha256(raw).hexdigest() != row['certificate_sha256']:
                raise RuntimeError('input continuation certificate hash mismatch')
            source = tpms_problem(row['family'], row['lambda_bounds'], row['c_bounds'])
            original = json.loads(gzip.decompress(raw))
            if original['problem'] != source.specification() or original['problem_sha256'] != source.fingerprint:
                raise RuntimeError('input continuation problem mismatch')
            result = attempt(row)
        results.append(result)
        write_json(args.out / 'event_records.json', results)
        if number % 10 == 0 or number + 1 == len(selected):
            print('CRITICAL_EVENT_PROGRESS ' + json.dumps({'processed': len(results), 'selected': len(selected),
                  'statuses': dict(Counter(r['status'] for r in results))}), flush=True)
    summary = {'selected_cells': len(selected), 'recorded_cells': len(results),
               'status_counts': dict(Counter(r['status'] for r in results)), 'seconds': monotonic()-started,
               'all_isotopy_regions_classified': False, 'topology_changing_events_proved': False}
    write_json(args.out / 'event_summary.json', summary)
    print('CERTIFIED_EVENT_SUMMARY ' + json.dumps(summary, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
