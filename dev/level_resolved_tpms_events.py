#!/usr/bin/env python3
"""Certified event search with the interpolation parameter as an unknown.

The older search fixes lambda and encloses the resulting critical level. Here
we fix a level in the requested cell and solve for (x,y,z,lambda), or
(x,y,z,mu,lambda) on the wall. A strict interval Krawczyk inclusion and norm <1
prove existence/uniqueness of that square-system root. Query membership and
bulk-domain inclusion are independently checked. This proves event existence,
not source isotopy classification or a global topology change.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
from time import monotonic

import numpy as np
from scipy.optimize import root
import sympy as sp

from knotted_graph.core.field_isotopy import (
    FieldProblem, _Evaluator, _ZERO, _add, _sub, _mul, _up, _square, tpms_problem,
)
from knotted_graph.core.field_isotopy_rational import _RationalEvaluator
from run_full_resolution_maps import dump


class SquareSystem:
    def __init__(self, expressions, variables):
        self.expressions = tuple(map(sp.sympify, expressions))
        self.variables = tuple(variables)
        self.n = len(self.variables)
        if (not 3 <= self.n <= 8 or len(self.expressions) != self.n
                or len(set(self.variables)) != self.n
                or any(not e.free_symbols.issubset(set(self.variables)) for e in self.expressions)):
            raise ValueError('invalid square analytic system')
        self.spec = {'equations': [sp.srepr(e) for e in self.expressions],
                     'variables': [sp.srepr(v) for v in self.variables]}
        self.digest = sha256(json.dumps(self.spec, sort_keys=True).encode()).hexdigest()
        matrix = sp.Matrix(self.expressions)
        self.numeric = sp.lambdify(self.variables, matrix, 'numpy')
        self.jacobian = sp.lambdify(self.variables, matrix.jacobian(self.variables), 'numpy')

    def evaluators(self, box, backend):
        cls = _Evaluator if backend == 'primary' else _RationalEvaluator
        evaluators = []
        for expression in self.expressions:
            evaluator = cls(FieldProblem(expression, self.variables, tuple(box), domain='box'))
            # Extend the existing safe expression DAG to all unknowns, including
            # lambda and the Lagrange multiplier. No serialized source is eval'd.
            evaluator.outputs = [evaluator._compile(expression)] + [
                evaluator._compile(sp.diff(expression, v)) for v in self.variables]
            evaluators.append(evaluator)
        return evaluators


def replay(system, certificate, backend='rational'):
    if backend not in ('primary', 'rational'):
        raise ValueError('unknown interval backend')
    if (certificate.get('schema') != 'knottedgraph.square_system.v1'
            or certificate.get('system') != system.spec
            or certificate.get('system_sha256') != system.digest):
        return {'valid': False, 'reason': 'system_mismatch'}
    try:
        box = tuple(tuple(map(float.fromhex, b)) for b in certificate['box_hex'])
        center = tuple(map(float.fromhex, certificate['center_hex']))
        inverse = np.array([[float.fromhex(x) for x in r] for r in certificate['inverse_hex']])
        if (len(box) != system.n or len(center) != system.n or inverse.shape != (system.n, system.n)
                or not np.isfinite(box).all() or not np.isfinite(center).all() or not np.isfinite(inverse).all()
                or any(not a < x < b for (a,b),x in zip(box,center))):
            raise ValueError('malformed box or matrix')
        evaluators = system.evaluators(box, backend)
        point = tuple((x,x) for x in center)
        f0 = [ev(point)[0] for ev in evaluators]
        jac = [ev(box)[1] for ev in evaluators]
        remainder = []
        for i in range(system.n):
            row = []
            for j in range(system.n):
                term = _ZERO
                for k in range(system.n):
                    c = float(inverse[i,k])
                    term = _add(term, _mul((c,c),jac[k][j]))
                row.append(_sub((1.,1.) if i == j else _ZERO, term))
            remainder.append(row)
        norm = 0.
        for row in remainder:
            total = 0.
            for a,b in row:
                total = _up(total+max(abs(a),abs(b)))
            norm = max(norm,total)
        if not math.isfinite(norm) or norm >= 1.:
            return {'valid':False,'reason':'contraction_not_established'}
        differences = [_sub(b,(x,x)) for b,x in zip(box,center)]
        image = []
        for i in range(system.n):
            correction = _ZERO
            for j in range(system.n):
                c = float(inverse[i,j])
                correction = _add(correction,_mul((c,c),f0[j]))
            value = _sub((center[i],center[i]),correction)
            for j in range(system.n):
                value = _add(value,_mul(remainder[i][j],differences[j]))
            image.append(value)
        if not all(a < lo <= hi < b for (a,b),(lo,hi) in zip(box,image)):
            return {'valid':False,'reason':'strict_inclusion_not_established'}
        return {'valid':True,'reason':'strict_inclusion_and_contraction','backend':backend,
                'norm_upper':norm,'root_enclosure_hex':[[a.hex(),b.hex()] for a,b in image]}
    except (ValueError, KeyError, TypeError, OverflowError):
        return {'valid':False,'reason':'malformed_or_unsupported_certificate'}


def certify(system, center, box):
    center=np.asarray(center,dtype=float)
    if center.shape != (system.n,) or not np.isfinite(center).all():
        return None
    try:
        inverse=np.linalg.inv(np.asarray(system.jacobian(*center),dtype=float))
    except (np.linalg.LinAlgError,ValueError,OverflowError):
        return None
    certificate={'schema':'knottedgraph.square_system.v1','system':system.spec,'system_sha256':system.digest,
                 'box_hex':[[float(a).hex(),float(b).hex()] for a,b in box],
                 'center_hex':[float(x).hex() for x in center],
                 'inverse_hex':[[float(x).hex() for x in r] for r in inverse]}
    checks={backend:replay(system,certificate,backend) for backend in ('primary','rational')}
    if not all(c['valid'] for c in checks.values()):
        return None
    certificate['replays']=checks
    return certificate


@lru_cache(maxsize=512)
def event_system(family, level, kind):
    source=tpms_problem(family,(0.,1.),(0.,.3))
    x,y,z,lam,c=source.variables
    field=sp.expand(source.expression+c)
    gradient=[sp.diff(field,v) for v in (x,y,z)]
    if kind == 'bulk':
        variables=(x,y,z,lam)
        equations=gradient+[field-sp.Rational(level)]
    else:
        mu=sp.Symbol('_event_mu',real=True)
        variables=(x,y,z,mu,lam)
        equations=[g-mu*v for g,v in zip(gradient,(x,y,z))]
        equations += [x*x+y*y+z*z-sp.Rational(source.radius)**2,field-sp.Rational(level)]
    return SquareSystem(equations,variables),source.radius,sp.lambdify((x,y,z,lam),gradient,'numpy')


def domain_membership(row, level, kind, certificate, radius):
    box=[tuple(map(float.fromhex,b)) for b in certificate['box_hex']]
    if not (row['c_bounds'][0] <= level <= row['c_bounds'][1]
            and row['lambda_bounds'][0] <= box[-1][0] < box[-1][1] <= row['lambda_bounds'][1]):
        return False
    if kind == 'bulk':
        squared=_ZERO
        for interval in box[:3]:
            squared=_add(squared,_square(interval))
        return _sub(squared,_square((radius,radius)))[1] < 0
    return kind == 'wall'


def attempt(row):
    retained=row.get('unresolved_spatial_box') or {}
    if 'box_hex' not in retained:
        return {'id':row['id'],'status':'unknown','reason':'missing_retained_seed'}
    bounds=np.array([[float.fromhex(a),float.fromhex(b)] for a,b in retained['box_hex']])
    seed=bounds[:3].mean(axis=1)
    left,right=row['lambda_bounds']
    bottom,top=row['c_bounds']
    levels=list(dict.fromkeys((bottom/2+top/2, (3*bottom+top)/4, (bottom+3*top)/4, bottom, top)))
    lam0=left/2+right/2
    attempts=[]
    for level in levels:
        for kind in ('wall','bulk'):
            system,radius,gradient=event_system(row['family'],level,kind)
            initial=seed.copy()
            if kind == 'wall':
                norm=np.linalg.norm(initial)
                initial=initial*(radius/norm) if norm > 1e-12 else np.array([radius,0.,0.])
                mu=float(np.dot(initial,np.asarray(gradient(*initial,lam0)).ravel())/radius**2)
                initial=np.r_[initial,mu,lam0]
            else:
                initial=np.r_[initial,lam0]
            tried={'level_hex':float(level).hex(),'kind':kind}
            try:
                numerical=root(lambda v:np.asarray(system.numeric(*v),dtype=float).ravel(),initial,
                               jac=lambda v:np.asarray(system.jacobian(*v),dtype=float),
                               method='hybr',options={'xtol':1e-10,'maxfev':160})
                center=np.asarray(numerical.x,dtype=float)
                residual=np.asarray(system.numeric(*center),dtype=float).ravel()
                if not np.isfinite(center).all() or not np.isfinite(residual).all() or np.max(np.abs(residual))>1e-7:
                    tried['status']='candidate_not_converged'
                elif not left < center[-1] < right:
                    tried['status']='candidate_outside_lambda_cell'
                else:
                    for width in (1e-5,1e-7,1e-9):
                        widths=np.full(system.n,width)
                        widths[-1]=min(width,(center[-1]-left)/2,(right-center[-1])/2)
                        box=[(float(x-w),float(x+w)) for x,w in zip(center,widths)]
                        certificate=certify(system,center,box)
                        if certificate is not None and domain_membership(row,level,kind,certificate,radius):
                            # Verification never relies on saved success flags.
                            fresh={b:replay(system,certificate,b) for b in ('primary','rational')}
                            if not all(v['valid'] for v in fresh.values()):
                                raise RuntimeError('fresh square-system replay failed')
                            return {'id':row['id'],'family':row['family'],'lambda_bounds':row['lambda_bounds'],
                                    'c_bounds':row['c_bounds'],'status':'certified_level_resolved_event',
                                    'kind':kind,'fixed_level_hex':float(level).hex(),'certificate':certificate,
                                    'fresh_replays':fresh,'domain_membership_verified':True,
                                    'source_topology_change_proved':False,'attempts':attempts}
                    tried['status']='interval_certificate_not_obtained'
            except (ValueError,OverflowError,FloatingPointError,np.linalg.LinAlgError) as exc:
                tried.update(status='candidate_error',error=f'{type(exc).__name__}: {exc}')
            attempts.append(tried)
    return {'id':row['id'],'family':row['family'],'status':'unknown','attempts':attempts}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists(): parser.error('use a new output directory')
    args.out.mkdir(parents=True)
    original=json.loads((args.input/'event_records.json').read_text())
    leaves={r['id']:r for r in json.loads((args.input/'partition/adaptive_leaves.json').read_text())}
    selected=[r for r in original if r['status']=='unknown']
    dump(args.out/'plan.json',{'selected_ids':[r['id'] for r in selected],
         'source_records_sha256':sha256((args.input/'event_records.json').read_bytes()).hexdigest(),
         'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()})
    results=[]
    started=monotonic()
    dump(args.out/'records.json',results)
    for row in selected:
        results.append(attempt(leaves[row['id']]))
        dump(args.out/'records.json',results)
        if len(results)%10==0 or len(results)==len(selected):
            print('LEVEL_RESOLVED_EVENTS '+json.dumps({'processed':len(results),'planned':len(selected),
                  'statuses':dict(Counter(r['status'] for r in results))}),flush=True)
    summary={'requested':len(selected),'recorded':len(results),'seconds':monotonic()-started,
             'statuses':dict(Counter(r['status'] for r in results)),
             'by_family':{f:dict(Counter(r['status'] for r in results if leaves[r['id']]['family']==f))
                          for f in sorted({r['family'] for r in selected})},
             'all_source_isotopy_classes_classified':False,'voxel_correspondence_proved':False}
    dump(args.out/'summary.json',summary)
    print('LEVEL_RESOLVED_SUMMARY '+json.dumps(summary),flush=True)


if __name__=='__main__': main()
