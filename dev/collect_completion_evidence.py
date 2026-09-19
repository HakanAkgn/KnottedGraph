#!/usr/bin/env python3
"""Collect exact-grid retries, replay event certificates and audit voxel links.

All original records remain intact. Successful retries are merged only into a
new derivative record array after checking source/geometry identities. Every
new level-resolved event is independently replayed against its original field.
The geometric and continuum checks are separate; an unproved bridge is never
filled by polynomial equality or by a digital-homology match.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import shutil
import subprocess
from time import monotonic

import numpy as np
import sympy as sp

from run_full_resolution_maps import dump, plan
from level_resolved_tpms_events import event_system, replay, domain_membership
from voxel_link_audit import audit_mask, pattern_table
from knotted_graph.core.field_isotopy import tpms_problem, verify
from knotted_graph.core.field_isotopy_rational import verify_rational


def terms(value):
    variable = sp.Symbol('A')
    expression = sp.expand(sp.sympify(value, locals={'A': variable}))
    if expression == 0:
        return []
    answer = []
    for term in sp.Add.make_args(expression):
        power = term.as_powers_dict().get(variable, sp.Integer(0))
        coefficient = sp.simplify(term / variable**power)
        if power.is_Integer is not True or coefficient.is_Integer is not True:
            raise ValueError('noninteger polynomial output')
        answer.append([int(power), int(coefficient)])
    if min(p for p,c in answer if c) != 0:
        raise ValueError('nonzero normalized output has nonzero minimum degree')
    return sorted(answer)


def merge_retries(original, directory, mode):
    original = Path(original)
    rows = json.loads(original.read_text())
    expected = plan(mode)
    if len(rows) != len(expected) or [r['id'] for r in rows] != [r['id'] for r in expected]:
        raise ValueError('original grid coverage differs')
    if any(any(r.get(k) != v for k,v in e.items()) for r,e in zip(rows,expected)):
        raise ValueError('original grid coordinates differ')
    pending = [r for r in rows if r['status'] == 'time_budget']
    known = {r['id']: r for r in rows}
    updated = {r['id']: dict(r, baseline_status=r['status']) for r in rows}
    observed = set()
    commits = set()
    original_hash = sha256(original.read_bytes()).hexdigest()
    for shard in range(10):
        folder = Path(directory)/f'archived-spatial-retry-{shard}'/mode
        spec = json.loads((folder/'plan.json').read_text())
        results = json.loads((folder/'records.json').read_text())
        wanted = [r['id'] for r in pending[shard::10]]
        if (spec['selected_ids'] != wanted or spec['input_records_sha256'] != original_hash
                or spec['shard'] != shard or spec['shards'] != 10
                or [r['id'] for r in results] != wanted):
            raise ValueError('retry shard coverage or source hash differs')
        commits.add(spec['commit'])
        for result in results:
            key = result['id']
            if key in observed:
                raise ValueError('duplicate retry')
            observed.add(key)
            base = known[key]
            target = updated[key]
            target['retry_result'] = result
            if result['status'] not in ('evaluated','fixed_diagram_evaluated'):
                target['last_attempt_status'] = result['status']
                continue
            if (result.get('geometry_changed') is not False or result['graph_sha256'] != base['graph_sha256']
                    or result['source_mask_sha256'] != base['source_mask_sha256']
                    or any(result.get(k) != base[k] for k in ('lambda_hex','level_hex','family','source'))):
                raise ValueError('successful retry is not the same geometric input')
            polynomial = terms(result['yamada'])
            if polynomial != result['terms']:
                raise ValueError('retry polynomial coefficients disagree')
            checks = result['projection_checks']
            if any(terms(c['value']) != polynomial for c in checks):
                raise ValueError('retry projection values disagree')
            if result['status'] == 'evaluated':
                if (not result['evaluation']['is_subcubic'] or len(checks) < 2
                        or len({tuple(c['projection']['rotation_angles']) for c in checks}) < 2):
                    raise ValueError('subcubic retry lacks distinct agreeing projections')
            target.update(status=result['status'], yamada=result['yamada'], evaluation=result['evaluation'],
                          projection_checks=checks, polynomial_terms=polynomial)
    if observed != {r['id'] for r in pending} or len(commits) != 1:
        raise ValueError('incomplete or mixed-protocol retry')
    final = [updated[r['id']] for r in rows]
    for row in final:
        if row['status'] in ('evaluated','fixed_diagram_evaluated'):
            row['polynomial_terms'] = terms(row['yamada'])
    return final, {'original_record_sha256': original_hash, 'retry_commit': commits.pop(),
                   'baseline_statuses':dict(Counter(r['status'] for r in rows)),
                   'final_statuses':dict(Counter(r['status'] for r in final)),
                   'retried_cells':len(pending), 'complete_exact_grid':True,
                   'by_family':{f:dict(Counter(r['status'] for r in final if r['family']==f))
                                for f in dict.fromkeys(r['family'] for r in final)}}


def collect_events(previous, additional):
    previous, additional = Path(previous), Path(additional)
    leaves = json.loads((previous/'partition/adaptive_leaves.json').read_text())
    events = json.loads((previous/'event_records.json').read_text())
    extra = json.loads((additional/'records.json').read_text())
    lookup = {r['id']:r for r in leaves}
    unknown = {r['id'] for r in events if r['status']=='unknown'}
    if len(extra) != len(unknown) or {r['id'] for r in extra} != unknown:
        raise ValueError('new event calculation does not cover every previously unknown cell')
    merged = {r['id']:dict(r) for r in events}
    for row in extra:
        target = lookup[row['id']]
        if row['status']=='certified_level_resolved_event':
            level = float.fromhex(row['fixed_level_hex'])
            system, radius, _ = event_system(target['family'],level,row['kind'])
            checks = {b: replay(system,row['certificate'],b) for b in ('primary','rational')}
            if (not all(c['valid'] for c in checks.values())
                    or not domain_membership(target,level,row['kind'],row['certificate'],radius)):
                raise ValueError('new event certificate failed fresh replay/domain membership')
            merged[row['id']].update(status=row['status'], level_resolved_result=row,
                                     completion_replays=checks)
        elif row['status'] != 'unknown':
            raise ValueError('unexpected additional event outcome')
    families = sorted({r['family'] for r in leaves})
    summary = {'regular_leaf_count':sum(r['status']=='certified' for r in leaves),
               'event_statuses':dict(Counter(r['status'] for r in merged.values())),
               'new_event_certificates_replayed':sum(r['status']=='certified_level_resolved_event' for r in extra),
               'by_family':{}, 'distinct_isotopy_classes_proved':False}
    for family in families:
        selected = [r for r in leaves if r['family']==family]
        total_area = Fraction(0)
        areas = {'regular':Fraction(0),'event':Fraction(0),'unknown':Fraction(0)}
        counts = Counter()
        for leaf in selected:
            l,r = map(Fraction,leaf['lambda_bounds'])
            b,t = map(Fraction,leaf['c_bounds'])
            area = (r-l)*(t-b)
            total_area += area
            state = 'regular' if leaf['status']=='certified' else (
                'unknown' if merged[leaf['id']]['status']=='unknown' else 'event')
            areas[state] += area
            counts[state] += 1
        if total_area != Fraction(float(.3)):
            raise ValueError('exact binary-endpoint cell areas do not cover the original domain')
        summary['by_family'][family] = {
            'leaf_counts':dict(counts),
            'area_fractions':{k:float(a/total_area) for k,a in areas.items()},
            'exact_area_fractions':{k:[(a/total_area).numerator,(a/total_area).denominator] for k,a in areas.items()},
            'event_area_refers_to_cells_not_the_critical_set':True}
    return leaves, list(merged.values()), summary


def local_audits(records, artifacts, mode, out):
    answers = []
    started = monotonic()
    for index, row in enumerate(records):
        directory = Path(artifacts)/row['artifact_relative_path']
        source = directory/'source.npz'
        with np.load(source,allow_pickle=False) as saved:
            mask = saved['mask']
            if sha256(mask.astype(np.uint8).tobytes()).hexdigest() != row['source_mask_sha256']:
                raise ValueError('source mask hash mismatch in link audit')
            result = audit_mask(mask)
        answers.append({'id':row['id'],'family':row['family'],
                        'source_mask_sha256':row['source_mask_sha256'],**result})
        if index % 200 == 0 or index+1==len(records):
            print('VOXEL_LINKS '+json.dumps({'mode':mode,'completed':len(answers),'total':len(records),
                  'nonmanifold_cases':sum(r['nonmanifold_lattice_vertices']>0 for r in answers)}),flush=True)
    dump(out/f'{mode}_voxel_links.json',answers)
    return {'cases':len(answers),'seconds':monotonic()-started,
            'pl_manifold_cases':sum(r['is_pl_3_manifold_with_boundary'] for r in answers),
            'nonmanifold_cases':sum(r['nonmanifold_lattice_vertices']>0 for r in answers),
            'by_family':{f:{'cases':sum(r['family']==f for r in answers),
                           'nonmanifold_cases':sum(r['family']==f and r['nonmanifold_lattice_vertices']>0 for r in answers)}
                         for f in dict.fromkeys(r['family'] for r in answers)}}


def continuum_digital_pairs(leaves, records, continuation, artifacts, out):
    conflicts, variations = [], []
    copied = set()
    for leaf in leaves:
        if leaf['status'] != 'certified':
            continue
        selected = [r for r in records if r['family']==leaf['family']
                    and leaf['lambda_bounds'][0] <= r['lambda'] <= leaf['lambda_bounds'][1]
                    and leaf['c_bounds'][0] <= r['level'] <= leaf['c_bounds'][1]]
        for first, second in combinations(selected,2):
            b1 = tuple(first['volume'][k] for k in ('connected_components','handle_rank','enclosed_voids'))
            b2 = tuple(second['volume'][k] for k in ('connected_components','handle_rank','enclosed_voids'))
            type_ = 'digital_homology_conflict' if b1 != b2 else (
                'spine_signature_variation' if all(r['status']=='evaluated' for r in (first,second))
                and first['polynomial_terms'] != second['polynomial_terms'] else None)
            if type_ is None:
                continue
            destination = conflicts if type_=='digital_homology_conflict' else variations
            if len(destination)>=12:
                continue
            problem = tpms_problem(leaf['family'],leaf['lambda_bounds'],leaf['c_bounds'])
            proof_path = Path(continuation)/leaf['certificate_artifact_path']
            raw = proof_path.read_bytes()
            if sha256(raw).hexdigest() != leaf['certificate_sha256']:
                raise ValueError('regularity proof hash mismatch')
            certificate = json.loads(gzip.decompress(raw))
            checks = {'primary':verify(problem,certificate),'rational':verify_rational(problem,certificate)}
            if not all(c['valid'] for c in checks.values()):
                raise ValueError('common analytic rectangle failed fresh proof replay')
            folder = out/'comparison_witnesses'
            folder.mkdir(exist_ok=True)
            target = folder/(leaf['id']+'.json.gz')
            if leaf['id'] not in copied:
                target.write_bytes(raw)
                copied.add(leaf['id'])
            pair = {'kind':type_,'family':leaf['family'],'regular_leaf_id':leaf['id'],
                    'lambda_bounds':leaf['lambda_bounds'],'c_bounds':leaf['c_bounds'],
                    'regularity_certificate':str(target.relative_to(out)),
                    'regularity_certificate_sha256':sha256(raw).hexdigest(),'fresh_replays':checks,
                    'first':{k:first.get(k) for k in ('id','lambda','level','source_mask_sha256','graph_sha256','volume','yamada')},
                    'second':{k:second.get(k) for k in ('id','lambda','level','source_mask_sha256','graph_sha256','volume','yamada')},
                    'first_betti':list(b1),'second_betti':list(b2)}
            for r in (first,second):
                source = Path(artifacts)/r['artifact_relative_path']
                for name in ('source.npz','graph.json.gz','collapse.npz'):
                    if (source/name).exists():
                        dest = folder/r['id']/name
                        dest.parent.mkdir(exist_ok=True)
                        if not dest.exists():
                            shutil.copyfile(source/name,dest)
            destination.append(pair)
    answer = {'digital_homology_conflict_examples':conflicts,'spine_signature_variation_examples':variations,
              'maximum_retained_examples_per_type':12,
              'scope':'Each pair lies in one freshly verified regular analytic rectangle. A different voxel Betti tuple obstructs simultaneous correct continuum reconstruction; different spine polynomials alone do not identify the source of disagreement.'}
    dump(out/'continuum_vs_digital_witnesses.json',answer)
    return {'retained_digital_homology_conflicts':len(conflicts),'retained_spine_signature_variations':len(variations)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('maps','retries','events','additional','artifacts','continuation','out'):
        parser.add_argument('--'+option,type=Path,required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    summary = {'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
               'original_sources_modified':False,'all_word_formulas_modified':False,'cavity_implementation_modified':False,
               'analytic_to_voxel_correspondence_proved':False,'complete_source_isotopy_classification_proved':False,
               'maps':{},'voxel_links':{}}
    modes = {}
    for mode in ('hamiltonian','tpms'):
        rows, stats = merge_retries(args.maps/f'{mode}_records.json',args.retries,mode)
        modes[mode] = rows
        summary['maps'][mode] = stats
        dump(args.out/f'{mode}_records.json',rows)
        print('MERGED_MAP_RESULT '+json.dumps({mode:stats}),flush=True)
    leaves, events, stats = collect_events(args.events,args.additional)
    summary['events'] = stats
    dump(args.out/'adaptive_leaves.json',leaves)
    dump(args.out/'event_records.json',events)
    dump(args.out/'event_summary.json',stats)
    print('MERGED_EVENT_RESULT '+json.dumps(stats),flush=True)
    dump(args.out/'octant_link_lookup.json',[{'pattern':i,'kind':v} for i,v in enumerate(pattern_table())])
    for mode in ('tpms','hamiltonian'):
        summary['voxel_links'][mode] = local_audits(modes[mode],args.artifacts,mode,args.out)
    summary['continuum_vs_digital'] = continuum_digital_pairs(
        leaves,modes['tpms'],args.continuation,args.artifacts,args.out)
    dump(args.out/'completion_summary.json',summary)
    dump(args.out/'data_sha256.json',{p.name:sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.json')})
    print('COMPLETION_EVIDENCE '+json.dumps(summary,sort_keys=True),flush=True)


if __name__ == '__main__':
    main()
