#!/usr/bin/env python3
"""Validate the complete adaptive numerical partition; no plotting or editing.

The per-leaf runs have already executed both interval replays. This collector
checks immutable certificates, exact source identities, dyadic tree coverage,
and same-revision provenance. It does not turn disconnected cover components
into distinct source-topology classes or make any graph-spine correspondence.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label

from certify_tpms_continuation import build_plan, write_json
from refine_tpms_continuation import subdivide
from knotted_graph.core.field_isotopy import tpms_problem


def audit(root, out, shards=20):
    root, out = Path(root), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    original = {r['id']: r for r in build_plan('grid', 21)}
    records, revisions, depths = {}, set(), set()
    input_hashes = {}
    for shard in range(shards):
        directory = root / f'refined-tpms-shard-{shard}'
        plan = json.loads((directory / 'plan.json').read_text())
        summary = json.loads((directory / 'summary.json').read_text())
        rows = json.loads((directory / 'records.json').read_text())
        if (plan['shard'] != shard or plan['shards'] != shards or not summary['complete_shard']
                or [r['id'] for r in plan['roots']] != list(original)[shard::shards]):
            raise ValueError('incomplete or changed adaptive shard')
        revisions.add(plan['commit'])
        depths.add(plan['max_parameter_depth'])
        input_hashes[directory.name] = {name: sha256((directory / name).read_bytes()).hexdigest()
                                       for name in ('plan.json', 'records.json', 'summary.json')}
        for row in rows:
            key = row['id']
            if key in records:
                raise ValueError('duplicate adaptive record')
            certificate_path = directory / row['certificate']
            raw = certificate_path.read_bytes()
            if sha256(raw).hexdigest() != row['certificate_sha256']:
                raise ValueError('altered certificate archive')
            certificate = json.loads(gzip.decompress(raw))
            problem = tpms_problem(row['family'], row['lambda_bounds'], row['c_bounds'])
            if (certificate['problem'] != problem.specification()
                    or certificate['problem_sha256'] != problem.fingerprint
                    or row['problem_sha256'] != problem.fingerprint
                    or certificate['status'] != row['status']):
                raise ValueError('certificate refers to a different analytic source or region')
            if row['status'] == 'certified':
                if (not row['primary_replay']['valid'] or not row['rational_replay']['valid']
                        or row['children']):
                    raise ValueError('accepted leaf lacks the recorded dual replay')
            row['certificate_artifact_path'] = str(certificate_path.relative_to(root))
            records[key] = row
    if len(revisions) != 1 or len(depths) != 1:
        raise ValueError('mixed source revisions or refinement depths')
    depth = depths.pop()
    if not 0 <= depth <= 4:
        raise ValueError('unsupported recorded depth')
    visited = set()
    for key, source in original.items():
        expected_root = {**source, 'depth': 0, 'tile_x': source['lambda_index'], 'tile_y': source['c_index']}
        stack = [(key, expected_root)]
        while stack:
            name, expected = stack.pop()
            row = records[name]
            if name in visited or any(row.get(k) != v for k, v in expected.items()):
                raise ValueError('tree repeats or changes a parameter region')
            visited.add(name)
            if row['children']:
                children = subdivide(expected)
                if row['children'] != [r['id'] for r in children] or row['depth'] >= depth:
                    raise ValueError('incomplete or over-deep subdivision')
                stack.extend((r['id'], r) for r in children)
            elif row['status'] != 'certified' and row['depth'] != depth:
                raise ValueError('unresolved leaf ended before the declared refinement depth')
    if visited != set(records):
        raise ValueError('orphan parameter records')
    leaves = [r for r in records.values() if not r['children']]
    families = sorted({r['family'] for r in original.values()})
    by_family = {}
    for family in families:
        cover = np.zeros((20 * 2**depth,) * 2, dtype=np.uint8)
        accepted = np.zeros_like(cover, dtype=bool)
        for row in leaves:
            if row['family'] != family:
                continue
            width = 2 ** (depth - row['depth'])
            x, y = row['tile_x'] * width, row['tile_y'] * width
            region = np.s_[y:y+width, x:x+width]
            cover[region] += 1
            accepted[region] = row['status'] == 'certified'
        if not np.all(cover == 1):
            raise ValueError('leaf rectangles do not cover the original square exactly once')
        # Closed verified blocks touching at a corner share that exact parameter.
        _, components = label(accepted, structure=np.ones((3, 3), dtype=int))
        subset = [r for r in leaves if r['family'] == family]
        root_rows = [r for r in records.values() if r['family'] == family and r['depth'] == 0]
        by_family[family] = {
            'root_cells': len(root_rows),
            'certified_root_cells': sum(r['status'] == 'certified' for r in root_rows),
            'leaf_status_counts': dict(Counter(r['status'] for r in subset)),
            'certified_parameter_area_fraction': float(accepted.mean()),
            'unresolved_parameter_area_fraction': float((~accepted).mean()),
            'connected_certified_cover_components': int(components),
            'distinct_isotopy_classes_proved': False,
        }
    summary = {'complete_partition': True, 'root_cells': 1200, 'parameter_depth': depth,
               'source_commit': revisions.pop(), 'total_records': len(records),
               'leaf_counts': dict(Counter(r['status'] for r in leaves)), 'by_family': by_family,
               'replays_executed_in_per_leaf_runs': True, 'replays_reexecuted_by_this_collector': False,
               'all_regions_classified_up_to_isotopy': False,
               'input_files': input_hashes}
    write_json(out / 'adaptive_summary.json', summary)
    write_json(out / 'adaptive_leaves.json', sorted(leaves, key=lambda r: r['id']))
    print('ADAPTIVE_COMPLETE_AUDIT ' + json.dumps(summary, sort_keys=True), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--shards', type=int, default=20)
    args = parser.parse_args()
    audit(args.root, args.out, args.shards)
