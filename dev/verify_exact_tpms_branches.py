#!/usr/bin/env python3
"""Exact stationary branches already contained in the original P-to-D family.

This is symbolic verification of explicit candidate sets, not an enumeration
of all stationary points or an alteration of the defining field. It explains
why the discriminant cannot be removed by increasing a regularity budget.
No graph-spine, full isotopy classification, cavity-treatment or figure claim.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from knotted_graph.core.field_isotopy import tpms_problem
from knotted_graph.core.field_isotopy_rational import _pi_bounds


def verify_identities():
    source = tpms_problem('schwarz_p_to_diamond', (0., 1.), (0., .3))
    x, y, z, lam, level = source.variables
    variables = (x, y, z)
    field = sp.expand(source.expression + level)
    gradient = sp.Matrix([sp.diff(field, v) for v in variables])
    hessian = sp.hessian(field, variables)
    first, second = [], []
    for axis in range(3):
        for sign in (-1, 1):
            point = [sp.Integer(0)] * 3
            point[axis] = sign * sp.pi
            first.append(tuple(point))
        for sign1 in (-1, 1):
            for sign2 in (-1, 1):
                point = [sp.Integer(0)] * 3
                other = [a for a in range(3) if a != axis]
                point[other[0]], point[other[1]] = sign1 * sp.pi, sign2 * sp.pi
                second.append(tuple(point))
    for points, target in ((first, 1 - 2*lam), (second, 2*lam - 1)):
        for point in points:
            substitutions = dict(zip(variables, point))
            if any(sp.simplify(v) != 0 for v in gradient.subs(substitutions)):
                raise AssertionError('candidate is not stationary for every lambda')
            if sp.simplify(field.subs(substitutions) - target) != 0:
                raise AssertionError('critical-value identity failed')
    first_hessian = hessian.subs({x: sp.pi, y: 0, z: 0}).applyfunc(sp.simplify)
    second_hessian = hessian.subs({x: sp.pi, y: sp.pi, z: 0}).applyfunc(sp.simplify)
    if first_hessian != sp.diag(1, 2*lam - 1, 2*lam - 1):
        raise AssertionError('first critical Hessian identity failed')
    if second_hessian != sp.diag(1 - 2*lam, 1 - 2*lam, -1):
        raise AssertionError('second critical Hessian identity failed')
    restriction = {x: sp.pi, z: 0, lam: sp.Rational(1, 2)}
    if sp.simplify(field.subs(restriction)) != 0:
        raise AssertionError('critical line does not lie on the zero level')
    if any(sp.simplify(v) != 0 for v in gradient.subs(restriction)):
        raise AssertionError('entire claimed line is not stationary')
    normal = hessian.extract([0, 2], [0, 2]).subs(restriction).applyfunc(sp.simplify)
    determinant = sp.trigsimp(normal.det())
    if sp.trigsimp(determinant + sp.sin(y)**2 / 2) != 0:
        raise AssertionError('normal critical-line Hessian determinant failed')
    # Exact rational Machin-series enclosure proves all discrete representatives
    # lie strictly inside the supplied binary-radius ball, without decimal pi.
    pi_low, pi_high = _pi_bounds()
    radius = Fraction(source.radius)
    if not 2*pi_high*pi_high < radius*radius < 3*pi_low*pi_low:
        raise AssertionError('retained radius does not have the stated strict bounds')
    return {
        'source_problem_sha256': source.fingerprint, 'source_field': sp.srepr(field),
        'radius_hex': source.radius.hex(),
        'first_branch': {'critical_value': '1-2*lambda', 'representative': ['pi', '0', '0'],
                         'symmetry_related_points': len(set(first)),
                         'hessian_diagonal': ['1', '2*lambda-1', '2*lambda-1'],
                         'morse_index_for_lambda_below_half': 2},
        'second_branch': {'critical_value': '2*lambda-1', 'representative': ['pi', 'pi', '0'],
                          'symmetry_related_points': len(set(second)),
                          'hessian_diagonal': ['1-2*lambda', '1-2*lambda', '-1'],
                          'morse_index_for_lambda_above_half': 3},
        'degenerate_parameter': {'lambda': '1/2', 'c': '0',
                                'critical_line': '(pi,t,0) with pi^2+t^2<R^2, and its coordinate/sign images',
                                'all_gradient_components_identically_zero': True,
                                'normal_hessian': [[str(v) for v in normal.row(i)] for i in range(2)],
                                'normal_determinant': '-sin(t)^2/2',
                                'isolated_root_uniqueness_inapplicable': True},
        'rational_ball_inclusion_proved': True,
        'all_symbolic_residuals_zero': True,
        'all_critical_points_enumerated': False,
        'distinct_source_isotopy_classes_proved': False,
        'all_word_formulas_changed': False, 'existing_cavity_treatment_changed': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('use a new result path')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result = verify_identities()
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print('EXACT_TPMS_BRANCHES ' + json.dumps(result, sort_keys=True), flush=True)
