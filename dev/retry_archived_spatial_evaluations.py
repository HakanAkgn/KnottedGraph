#!/usr/bin/env python3
"""Re-evaluate unavailable polynomials on unchanged archived embedded graphs.

Only timed-out cells are selected. Each retained witness/graph hash is checked.
All candidate views are generic projections of that same graph. Two completed
subcubic values must agree; higher-valence results remain fixed-diagram values.
No original records, geometry, model, cavity route, or formulas are changed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from hashlib import sha256
import json
import multiprocessing as mp
import os
from pathlib import Path
import signal
import subprocess
from time import monotonic
import warnings

import networkx as nx
import numpy as np
import sympy as sp

from run_full_resolution_maps import dump
from knotted_graph.projection.pd_code import sample_projections


def load_graph(directory, record):
    directory = Path(directory)
    raw = (directory / 'graph.json.gz').read_bytes()
    if sha256(raw).hexdigest() != record['graph_sha256']:
        raise ValueError('retained graph hash differs')
    witness = directory / 'collapse.npz'
    if sha256(witness.read_bytes()).hexdigest() != record['reconstruction']['archive_sha256']:
        raise ValueError('retained collapse hash differs')
    data = json.loads(gzip.decompress(raw))
    graph = nx.MultiGraph()
    for node in data['nodes']:
        p = np.asarray(node['pos'], dtype=float)
        if node['id'] in graph or p.shape != (3,) or not np.isfinite(p).all():
            raise ValueError('invalid retained node')
        graph.add_node(node['id'], pos=p)
    for edge in data['edges']:
        p = np.asarray(edge['pts'], dtype=float)
        u, v, k = edge['u'], edge['v'], edge['key']
        if (u not in graph or v not in graph or graph.has_edge(u, v, k)
                or p.ndim != 2 or p.shape[1] != 3 or len(p) < 2 or not np.isfinite(p).all()):
            raise ValueError('invalid retained edge')
        graph.add_edge(u, v, key=k, pts=p)
    components = nx.number_connected_components(graph)
    actual = {'vertices': len(graph), 'edges': graph.number_of_edges(), 'components': components,
              'cycle_rank': graph.number_of_edges()-len(graph)+components,
              'max_degree': max(dict(graph.degree()).values(), default=0)}
    if actual != record['graph']:
        raise ValueError('retained graph counts differ')
    return graph


def bounded(function, arguments, seconds, output):
    """File-based child result avoids pipe-capacity deadlocks for large values."""
    def target():
        try:
            answer = function(*arguments)
        except Exception as exc:
            answer = {'status': 'unavailable', 'error': f'{type(exc).__name__}: {exc}'}
        dump(output, answer)
    child = mp.get_context('fork').Process(target=target)
    child.start()
    child.join(seconds)
    if child.is_alive():
        child.kill()
        child.join()
        answer = {'status': 'time_budget', 'seconds_limit': seconds}
        dump(output, answer)
        return answer
    if not output.exists():
        return {'status': 'worker_error', 'exitcode': child.exitcode}
    return json.loads(output.read_text())


def polynomial(view):
    started = monotonic()
    value = sp.expand(view.processor.compute_yamada(sp.Symbol('A'), normalize=True, n_jobs=1))
    terms = []
    if value != 0:
        for term in sp.Add.make_args(value):
            power = term.as_powers_dict().get(sp.Symbol('A'), sp.Integer(0))
            coefficient = sp.simplify(term / sp.Symbol('A')**power)
            if power.is_Integer is not True or coefficient.is_Integer is not True:
                raise ValueError('result is not an integer Laurent polynomial')
            terms.append([int(power), int(coefficient)])
        if min(power for power, coefficient in terms if coefficient) != 0:
            raise ValueError('normalization did not give minimum degree zero')
    return {'status': 'evaluated', 'value': str(value), 'terms': sorted(terms),
            'seconds': monotonic()-started,
            'projection': {'rotation_angles': list(view.rotation_angles),
                           'rotation_order': view.rotation_order,
                           'num_crossings': view.num_crossings, 'pd_code': view.pd_code}}


def retry(record, root, directory, view_seconds):
    directory = Path(directory)
    os.setsid()
    started = monotonic()
    result = {k: record[k] for k in ('id', 'family', 'source', 'lambda', 'level', 'lambda_hex', 'level_hex')}
    result.update(status='started', original_status=record['status'], yamada=None,
                  graph_sha256=record.get('graph_sha256'), source_mask_sha256=record.get('source_mask_sha256'),
                  geometry_changed=False, attempts=[], source_isotopy_classification_proved=False)
    dump(directory / 'record.json', result)
    try:
        graph = load_graph(Path(root) / record['artifact_relative_path'], record)
        subcubic = max(dict(graph.degree()).values(), default=0) <= 3
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            views = sample_projections(graph, num_rotation_samples=40)
        views = [v for _, v in sorted(enumerate(views), key=lambda p: (p[1].num_crossings, p[0]))]
        result['candidate_projections'] = [{'angles': list(v.rotation_angles), 'crossings': v.num_crossings} for v in views]
        result['warnings'] = [str(w.message) for w in captured]
        result['is_subcubic'] = subcubic
        dump(directory / 'record.json', result)
        completed = []
        for index, view in enumerate(views[:6]):
            answer = bounded(polynomial, (view,), view_seconds, directory / f'view-{index}.json')
            result['attempts'].append(answer)
            if answer['status'] == 'evaluated':
                completed.append(answer)
            if len(completed) >= (2 if subcubic else 1):
                if subcubic and completed[0]['terms'] != completed[1]['terms']:
                    result.update(status='projection_disagreement', error='two normalized generic projections disagree')
                    break
                first = completed[0]
                result.update(status='evaluated' if subcubic else 'fixed_diagram_evaluated',
                              yamada=first['value'], terms=first['terms'],
                              evaluation={'evaluation_kind': 'spatial-yamada' if subcubic else 'diagram-yamada',
                                          'normalization': 'signed-minimum-degree-zero',
                                          'projection': first['projection'], 'is_subcubic': subcubic},
                              projection_checks=[{'value': r['value'], 'projection': r['projection']} for r in completed])
                break
            dump(directory / 'record.json', result)
        if result['status'] == 'started':
            result.update(status='unavailable', error='insufficient completed generic projections')
    except Exception as exc:
        result.update(status='unavailable', error=f'{type(exc).__name__}: {exc}', yamada=None)
    finally:
        result['seconds'] = monotonic()-started
        dump(directory / 'record.json', result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--records', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--shard', type=int, required=True)
    parser.add_argument('--shards', type=int, default=10)
    parser.add_argument('--view-seconds', type=int, default=30)
    parser.add_argument('--case-seconds', type=int, default=210)
    args = parser.parse_args()
    if args.out.exists() or not 0 <= args.shard < args.shards:
        parser.error('new output directory and valid shard required')
    args.out.mkdir(parents=True)
    all_rows = json.loads(args.records.read_text())
    pending = [r for r in all_rows if r['status'] == 'time_budget']
    selected = pending[args.shard::args.shards]
    dump(args.out / 'plan.json', {'input_records_sha256': sha256(args.records.read_bytes()).hexdigest(),
         'all_timeout_ids': [r['id'] for r in pending], 'selected_ids': [r['id'] for r in selected],
         'shard': args.shard, 'shards': args.shards, 'views': 40, 'maximum_evaluations': 6,
         'view_seconds': args.view_seconds, 'case_seconds': args.case_seconds,
         'commit': subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip()})
    results = []
    dump(args.out / 'records.json', results)
    for record in selected:
        directory = args.out / record['id']
        directory.mkdir()
        child = mp.get_context('fork').Process(target=retry, args=(record, str(args.root), str(directory), args.view_seconds))
        child.start()
        child.join(args.case_seconds)
        if child.is_alive():
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                child.kill()
            child.join()
            result = json.loads((directory/'record.json').read_text()) if (directory/'record.json').exists() else {'id':record['id']}
            result.update(status='time_budget', yamada=None, error='ranked-view case budget reached')
            dump(directory/'record.json', result)
        elif (directory/'record.json').exists():
            result = json.loads((directory/'record.json').read_text())
        else:
            result = {'id':record['id'], 'status':'worker_error', 'yamada':None}
        results.append(result)
        dump(args.out / 'records.json', results)
        print('ARCHIVED_RETRY ' + json.dumps({'id':result['id'],'status':result['status'],
              'crossings':[v['crossings'] for v in result.get('candidate_projections',[])[:6]]}), flush=True)
    dump(args.out / 'summary.json', {'selected':len(selected),'recorded':len(results),
         'statuses':dict(Counter(r['status'] for r in results)), 'original_geometry_changed':False})


if __name__ == '__main__':
    main()
