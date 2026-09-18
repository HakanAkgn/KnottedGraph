"""Adversarial checks for whole-domain, boundary-aware continuation."""
from copy import deepcopy
from fractions import Fraction
import math
import random

import mpmath as mp
import pytest
import sympy as sp

from knotted_graph.core.field_isotopy import (
    FieldProblem, _Evaluator, _add, _mul, _power, _rational,
    certify, tpms_problem, verify,
)

x, y, z, p = sp.symbols("x y z p", real=True)
VARIABLES = (x, y, z, p)


def ball(expr, pb=(0.0, 1.0), radius=3.0):
    return FieldProblem(expr, VARIABLES, ((-radius, radius),)*3+(pb,), "ball", radius)


def box(expr, pb=(0.0, 1.0)):
    return FieldProblem(expr, VARIABLES, ((-1.0, 1.0),)*3+(pb,), "box")


def test_rational_conversion_encloses_exact_values():
    for numerator in range(-31, 32):
        for denominator in (1, 3, 7, 11, 2**80):
            lo, hi = _rational(numerator, denominator)
            q = Fraction(numerator, denominator)
            assert Fraction(lo) <= q <= Fraction(hi)


@pytest.mark.parametrize("operation", ["add", "multiply"])
def test_arithmetic_encloses_exact_rational_endpoint_operations(operation):
    rng = random.Random(1417)
    for _ in range(300):
        a = tuple(sorted(rng.uniform(-10, 10) for _ in range(2)))
        b = tuple(sorted(rng.uniform(-10, 10) for _ in range(2)))
        lo, hi = (_add if operation == "add" else _mul)(a, b)
        for av in a:
            for bv in b:
                value = (Fraction(av)+Fraction(bv) if operation == "add"
                         else Fraction(av)*Fraction(bv))
                assert Fraction(lo) <= value <= Fraction(hi)


@pytest.mark.parametrize("power", range(9))
def test_powers_enclose_samples(power):
    for a in ((-2.1, 1.3), (-3., -1.1), (0., 0.), (.1, .8)):
        lo, hi = _power(a, power)
        for val in (a[0], (a[0]+a[1])/2, a[1]):
            exact = Fraction(val)**power
            assert Fraction(lo) <= exact <= Fraction(hi)


@pytest.mark.parametrize("function", [sp.sin, sp.cos])
def test_interval_trigonometry_encloses_high_precision_point_values(function):
    problem = box(function(x+p), (-0.2, 0.3))
    ev = _Evaluator(problem)
    rng = random.Random(21)
    with mp.workdps(90):
        for _ in range(50):
            lo = rng.uniform(-8, 8)
            hi = lo + rng.uniform(0, 4)
            value, _ = ev(((lo, hi), (-1., 1.), (-1., 1.), (-.2, .3)))
            for t in (0., .23, .75, 1.):
                v = mp.mpf(lo)+(mp.mpf(hi)-mp.mpf(lo))*t
                q = mp.mpf(-.2)+(mp.mpf(.3)-mp.mpf(-.2))*t
                actual = (mp.sin if function == sp.sin else mp.cos)(v+q)
                assert mp.mpf(value[0]) <= actual <= mp.mpf(value[1])


def test_box_plane_has_all_27_strata_and_verified_cover():
    problem = box(x-p, (-.25, .25))
    cert = certify(problem)
    assert cert["status"] == "certified"
    assert len(cert["trees"]) == 27
    assert verify(problem, cert)["valid"]
    assert not cert["periodic_identification"]
    assert not cert["graph_spine_certified"]


def test_ball_sphere_radius_continuation():
    problem = ball(x*x+y*y+z*z-p, (1., 2.))
    cert = certify(problem, max_boxes=30000)
    assert cert["status"] == "certified", cert
    assert verify(problem, cert)["valid"]


def test_cavity_bearing_shell_continuation_needs_no_graph_spine():
    r2 = x*x+y*y+z*z
    problem = ball((r2-p)*(r2-4), (0.8, 1.2))
    cert = certify(problem, max_boxes=100000)
    assert cert["status"] == "certified", cert
    assert verify(problem, cert)["valid"]


def test_interior_birth_is_not_certified():
    problem = ball(x*x+y*y+z*z-p, (-.1, .1))
    cert = certify(problem, max_boxes=4000, max_depth=36)
    assert cert["status"] == "unknown"
    assert not verify(problem, cert)["valid"]


def test_spherical_wall_tangency_is_not_hidden_by_nonzero_bulk_gradient():
    problem = ball(x-p, (.9, 1.1), radius=1.)
    cert = certify(problem, max_boxes=6000, max_depth=36)
    assert cert["status"] == "unknown"


def test_box_corner_event_requires_zero_dimensional_strata():
    problem = box(x+y+z-p, (2.9, 3.1))
    cert = certify(problem, max_boxes=10000, max_depth=28)
    assert cert["status"] == "unknown"
    assert cert["unresolved"]["stratum"] == "+++"


def test_identical_endpoint_topology_does_not_certify_a_singular_path():
    problem = ball(x*x+y*y+z*z-(p-sp.Rational(1, 2))**2+sp.Rational(1, 100))
    cert = certify(problem, max_boxes=6000, max_depth=36)
    assert cert["status"] == "unknown"


def test_constant_full_and_empty_domains():
    for expr in (sp.Integer(-1), sp.Integer(1)):
        problem = ball(expr)
        cert = certify(problem)
        assert cert["status"] == "certified"
        assert verify(problem, cert)["valid"]


def test_budget_exhaustion_cannot_pass():
    problem = ball(x*x+y*y+z*z-p, (1., 2.))
    cert = certify(problem, max_boxes=1)
    assert cert["status"] == "unknown"
    assert cert["reason"] == "box_budget"


@pytest.mark.parametrize("change", ["problem", "drop_stratum", "truncate", "extra", "token", "false_leaf"])
def test_replay_rejects_modified_certificates(change):
    problem = (ball(x*x+y*y+z*z-p, (1., 2.)) if change in ("truncate", "false_leaf")
               else box(x-p, (-.25, .25)))
    cert = deepcopy(certify(problem))
    assert cert["status"] == "certified"
    if change == "problem":
        cert["problem"]["domain"] = "other"
    elif change == "drop_stratum":
        cert["trees"].pop()
    elif change == "truncate":
        cert["trees"][0]["preorder"] = cert["trees"][0]["preorder"][:-1]
    elif change == "extra":
        cert["trees"][0]["preorder"] += "."
    elif change == "token":
        cert["trees"][0]["preorder"] = "?"
    else:
        cert["trees"][0]["preorder"] = "."
    assert not verify(problem, cert)["valid"]


def test_replay_matches_caller_supplied_problem_not_serialized_code():
    original = box(x-p, (-.25, .25))
    other = box(x*x-p, (-.25, .25))
    assert not verify(other, certify(original))["valid"]


@pytest.mark.parametrize("expr", [sp.exp(x), sp.Abs(x), sp.sqrt(x), 1/x, x**257])
def test_unsupported_operations_fail_explicitly(expr):
    with pytest.raises(ValueError):
        certify(box(expr))


@pytest.mark.parametrize("budget", [0, -1, True, 1.5])
def test_invalid_budgets(budget):
    with pytest.raises(ValueError):
        certify(box(x), max_boxes=budget)


def test_invalid_domain_and_symbols():
    with pytest.raises(ValueError):
        FieldProblem(x, VARIABLES, ((-1., 1.),)*4, "periodic")
    with pytest.raises(ValueError):
        FieldProblem(x, VARIABLES, ((-1., 1.),)*4, "ball", 2.)
    with pytest.raises(ValueError):
        FieldProblem(x+sp.Symbol("unknown"), VARIABLES, ((-1., 1.),)*4, "box")
    with pytest.raises(ValueError):
        FieldProblem(x, VARIABLES, ((-math.inf, 1.),)*4, "box")


def test_tpms_models_match_source_formulas():
    rng = random.Random(18)
    functions = {
        "gyroid": lambda a,b,c: math.sin(a)*math.cos(b)+math.sin(b)*math.cos(c)+math.sin(c)*math.cos(a),
        "schwarz_p": lambda a,b,c: math.cos(a)+math.cos(b)+math.cos(c),
        "diamond": lambda a,b,c: math.cos(a)*math.cos(b)*math.cos(c)-math.sin(a)*math.sin(b)*math.sin(c),
    }
    for family, ends in (("gyroid_to_diamond", ("gyroid", "diamond")),
                         ("gyroid_to_schwarz_p", ("gyroid", "schwarz_p")),
                         ("schwarz_p_to_diamond", ("schwarz_p", "diamond"))):
        problem = tpms_problem(family, (0., 1.), (0., .3))
        function = sp.lambdify(problem.variables, problem.expression, "math")
        for _ in range(20):
            xyz = tuple(rng.uniform(-5, 5) for _ in range(3))
            lam, c = rng.random(), .3*rng.random()
            reference = (1-lam)*functions[ends[0]](*xyz)+lam*functions[ends[1]](*xyz)-c
            assert function(*xyz, lam, c) == pytest.approx(reference, abs=1e-14)
        assert problem.radius.hex() == (0.72*(2.25*math.pi)).hex()


def test_explicit_primitive_diamond_critical_loci():
    problem = tpms_problem('schwarz_p_to_diamond', (0., 1.), (0., .3))
    a, b, d, lam, level = problem.variables
    for point, value in (((sp.pi, 0, 0), 1-2*lam),
                         ((sp.pi, sp.pi, 0), 2*lam-1)):
        substitutions = dict(zip((a, b, d), point))
        assert all(sp.simplify(sp.diff(problem.expression, v).subs(substitutions)) == 0
                   for v in (a, b, d))
        assert sp.simplify(problem.expression.subs(substitutions) - (value-level)) == 0
        assert float(sum(v*v for v in point)) < problem.radius**2
