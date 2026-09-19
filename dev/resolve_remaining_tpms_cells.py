#!/usr/bin/env python3
"""Try new spatial seeds without weakening event-certificate conditions.

Previously verified nearby critical points propose starts for unresolved cells.
Only a fresh square-system interval proof supplies a positive result. A separate
exact G-to-P degenerate stationary point is checked symbolically and by rational
domain membership. Counts refer to event-containing rectangles, not transitions.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
import subprocess
from time import monotonic

import sympy as sp

from level_resolved_tpms_events import attempt, event_system, replay, domain_membership
from run_full_resolution_maps import dump
from knotted_graph.core.field_isotopy import tpms_problem
from knotted_graph.core.field_isotopy_rational import _pi_bounds


@lru_cache(maxsize=1)
def exact_gyroid_schwarz_identity():
    source = tpms_problem('gyroid_to_schwarz_p',(0.,1.),(0.,.3))
    x,y,z,lam,c = source.variables
    field = source.expression+c
    replacements = {x:-sp.pi/2,y:-sp.pi/2,z:-sp.pi/2,lam:sp.Rational(1,2)}
    residuals = [sp.simplify(e.subs(replacements)) for e in (field,*[sp.diff(field,v) for v in (x,y,z)])]
    if any(v!=0 for v in residuals):
        raise ValueError('proposed exact G-to-P event is not stationary on the zero level')
    _,upper = _pi_bounds()
    if not 3*upper*upper/4 < Fraction(source.radius)**2:
        raise ValueError('exact G-to-P point is not proved inside the original ball')
    return {'field_sha256':source.fingerprint,'lambda_fraction':[1,2],'c_fraction':[0,1],
            'point':['-pi/2']*3,'all_four_symbolic_residuals_zero':True,
            'strict_ball_inclusion_by_rational_pi':True,'topology_change_proved':False}


def exact_member(leaf):
    if leaf['family']!='gyroid_to_schwarz_p':
        return False
    left,right = map(Fraction,leaf['lambda_bounds'])
    bottom,top = map(Fraction,leaf['c_bounds'])
    return left <= Fraction(1,2) <= right and bottom <= 0 <= top


def seed_pool(events):
    result = []
    for row in events:
        if row['status']=='unknown':
            continue
        if 'level_resolved_result' in row:
            certificate = row['level_resolved_result']['certificate']
            point = [float.fromhex(v) for v in certificate['center_hex'][:3]]
        elif 'certificate' in row:
            point = [float.fromhex(v) for v in row['certificate']['center_hex'][:3]]
        elif 'exact_branch_witness' in row:
            point = [math.pi if v=='pi' else 0. for v in row['exact_branch_witness']['point']]
        else:
            continue
        result.append({'family':row['family'],'point':point,'lambda':sum(row['lambda_bounds'])/2,
                       'level':sum(row['c_bounds'])/2,'source_event_id':row['id']})
    return result


def candidates(leaf,pool):
    lam,level = sum(leaf['lambda_bounds'])/2,sum(leaf['c_bounds'])/2
    nearby = sorted((r for r in pool if r['family']==leaf['family']),
                    key=lambda r:(r['lambda']-lam)**2+((r['level']-level)/.3)**2)
    points,seen = [],set()
    for row in nearby:
        key = tuple(round(v,4) for v in row['point'])
        if key in seen:
            continue
        seen.add(key)
        points.append((row['point'],row['source_event_id']))
        if len(points)>=20:
            break
    for value in (-3*math.pi/4,-math.pi/2,-math.pi/4,0.,math.pi/4,math.pi/2,3*math.pi/4):
        points.append(([value]*3,'deterministic_diagonal_start'))
    radius = tpms_problem(leaf['family'],leaf['lambda_bounds'],leaf['c_bounds']).radius
    for axis in range(3):
        for sign in (-1,1):
            point = [0.,0.,0.]
            point[axis] = sign*radius
            points.append((point,'deterministic_wall_axis_start'))
    return points


def search_one(leaf,pool,seconds):
    if exact_member(leaf):
        return {'id':leaf['id'],'family':leaf['family'],'status':'exact_degenerate_gyroid_schwarz_event',
                'exact_identity':exact_gyroid_schwarz_identity()}
    started = monotonic()
    tried = []
    for point,source_id in candidates(leaf,pool):
        if monotonic()-started > seconds:
            break
        proposal = deepcopy(leaf)
        bounds = [[float(v-1e-7).hex(),float(v+1e-7).hex()] for v in point]
        bounds += [[float(v).hex() for v in leaf[key]] for key in ('lambda_bounds','c_bounds')]
        proposal['unresolved_spatial_box'] = {'box_hex':bounds}
        answer = attempt(proposal)
        tried.append({'seed_source':source_id,'seed':point,'status':answer['status']})
        if answer['status']=='certified_level_resolved_event':
            answer['status'] = 'certified_multiseed_level_event'
            answer['seed_attempts'] = tried
            answer['seconds'] = monotonic()-started
            return answer
    return {'id':leaf['id'],'family':leaf['family'],'status':'unknown','seed_attempts':tried,
            'seconds':monotonic()-started,'seconds_budget':seconds}


def search(data,out,shard,shards):
    data,out = Path(data),Path(out)
    out.mkdir(parents=True,exist_ok=False)
    events = json.loads((data/'event_records.json').read_text())
    leaves = {r['id']:r for r in json.loads((data/'adaptive_leaves.json').read_text())}
    pending = [r for r in events if r['status']=='unknown']
    chosen = pending[shard::shards]
    pool = seed_pool(events)
    dump(out/'plan.json',{'input_sha256':sha256((data/'event_records.json').read_bytes()).hexdigest(),
         'selected_ids':[r['id'] for r in chosen],'shard':shard,'shards':shards,
         'case_seconds_budget':60,'maximum_neighbor_seeds':20,'additional_fixed_seeds':13,
         'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()})
    results = []
    dump(out/'records.json',results)
    for row in chosen:
        answer = search_one(leaves[row['id']],pool,60)
        results.append(answer)
        dump(out/'records.json',results)
        print('MULTISEED_EVENT '+json.dumps({'id':answer['id'],'status':answer['status'],
              'tried_starts':len(answer.get('seed_attempts',[]))}),flush=True)
    dump(out/'summary.json',{'selected':len(chosen),'completed':len(results),
         'statuses':dict(Counter(r['status'] for r in results))})


def collect(data,root,out):
    data,root,out = Path(data),Path(root),Path(out)
    shutil.copytree(data,out)
    events = json.loads((data/'event_records.json').read_text())
    leaves = json.loads((data/'adaptive_leaves.json').read_text())
    lookup = {r['id']:r for r in leaves}
    pending = [r for r in events if r['status']=='unknown']
    merged = {r['id']:r for r in events}
    old_hash = sha256((data/'event_records.json').read_bytes()).hexdigest()
    seen = set()
    commits = set()
    extra_success = 0
    exact_success = 0
    for shard in range(10):
        directory = root/f'final-event-retry-{shard}'
        plan = json.loads((directory/'plan.json').read_text())
        rows = json.loads((directory/'records.json').read_text())
        wanted = [r['id'] for r in pending[shard::10]]
        if plan['input_sha256']!=old_hash or plan['selected_ids']!=wanted or [r['id'] for r in rows]!=wanted:
            raise ValueError('multiseed record identity or coverage differs')
        commits.add(plan['commit'])
        for row in rows:
            key = row['id']
            if key in seen:
                raise ValueError('duplicate multiseed case')
            seen.add(key)
            target = lookup[key]
            if row['status']=='certified_multiseed_level_event':
                level = float.fromhex(row['fixed_level_hex'])
                system,radius,_ = event_system(target['family'],level,row['kind'])
                checks = {b:replay(system,row['certificate'],b) for b in ('primary','rational')}
                if not all(c['valid'] for c in checks.values()) or not domain_membership(target,level,row['kind'],row['certificate'],radius):
                    raise ValueError('multiseed event failed fresh interval/domain verification')
                merged[key].update(status=row['status'],multiseed_result=row,multiseed_replays=checks)
                extra_success += 1
            elif row['status']=='exact_degenerate_gyroid_schwarz_event':
                if not exact_member(target) or row['exact_identity']!=exact_gyroid_schwarz_identity():
                    raise ValueError('exact degenerate event identity or membership invalid')
                merged[key].update(status=row['status'],multiseed_result=row)
                exact_success += 1
            elif row['status']!='unknown':
                raise ValueError('unexpected multiseed status')
            else:
                merged[key]['multiseed_result'] = row
    if len(seen)!=len(pending) or len(commits)!=1:
        raise ValueError('incomplete or mixed-protocol multiseed results')
    summary = json.loads((data/'completion_summary.json').read_text())
    summary['multiseed_completion'] = {'input_unknown':len(pending),'additional_local_events':extra_success,
        'additional_exact_events':exact_success,'remaining_unknown':sum(r['status']=='unknown' for r in merged.values()),
        'search_commit':commits.pop(),'freshly_replayed':True}
    event_summary = summary['events']
    event_summary['event_statuses'] = dict(Counter(r['status'] for r in merged.values()))
    event_summary['new_event_certificates_replayed'] += extra_success
    for family in event_summary['by_family']:
        areas = {'regular':Fraction(0),'event':Fraction(0),'unknown':Fraction(0)}
        counts = Counter()
        for leaf in leaves:
            if leaf['family']!=family:
                continue
            left,right = map(Fraction,leaf['lambda_bounds'])
            bottom,top = map(Fraction,leaf['c_bounds'])
            key = 'regular' if leaf['status']=='certified' else ('unknown' if merged[leaf['id']]['status']=='unknown' else 'event')
            areas[key] += (right-left)*(top-bottom)
            counts[key] += 1
        total = sum(areas.values())
        event_summary['by_family'][family].update(leaf_counts=dict(counts),
            area_fractions={k:float(v/total) for k,v in areas.items()},
            exact_area_fractions={k:[(v/total).numerator,(v/total).denominator] for k,v in areas.items()})
    dump(out/'event_records.json',[merged[r['id']] for r in events])
    dump(out/'event_summary.json',event_summary)
    dump(out/'completion_summary.json',summary)
    dump(out/'data_sha256.json',{p.name:sha256(p.read_bytes()).hexdigest() for p in out.glob('*.json') if p.name!='data_sha256.json'})
    print('FINAL_COMPLETION_EVIDENCE '+json.dumps(summary,sort_keys=True),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--root',type=Path)
    parser.add_argument('--shard',type=int,default=0)
    arguments = parser.parse_args()
    if arguments.root is not None:
        collect(arguments.data,arguments.root,arguments.out)
    else:
        if not 0 <= arguments.shard < 10:
            parser.error('shard must lie in 0..9')
        search(arguments.data,arguments.out,arguments.shard,10)
