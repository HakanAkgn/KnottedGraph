from copy import deepcopy
import math

import pytest
import sympy as sp

from knotted_graph.core.critical_events import (
    CriticalPointProblem, certify_critical_point, verify_critical_point,
)

X, Y, Z = sp.symbols("x y z", real=True)
VARIABLES = (X, Y, Z)


def box(center, radius=1e-5):
    return [(float(x-radius), float(x+radius)) for x in center]


@pytest.mark.parametrize("backend", ["primary", "rational"])
@pytest.mark.parametrize("expression,center", [
    (X*X + Y*Y + Z*Z, (0., 0., 0.)),
    ((X-1)**2 + (Y+1)**2 - Z*Z, (1., -1., 0.)),
    (sp.sin(X) + Y*Y + Z*Z, (math.pi/2, 0., 0.)),
    (sp.cos(X) + sp.cos(Y) + sp.cos(Z), (0., 0., 0.)),
])
def test_bulk_local_root_exists_and_replays(expression, center, backend):
    problem = CriticalPointProblem(expression, VARIABLES, 5.)
    certificate = certify_critical_point(problem, box(center))
    assert certificate["status"] == "certified"
    result = verify_critical_point(problem, certificate, backend=backend)
    assert result["valid"]
    assert result["contraction_norm_upper"] < 1
    assert not result["source_topology_change_certified"]


@pytest.mark.parametrize("backend", ["primary", "rational"])
def test_spherical_wall_maximum(backend):
    problem = CriticalPointProblem(Z, VARIABLES, 2., "wall")
    certificate = certify_critical_point(problem, box((0., 0., 2., .5)))
    assert certificate["status"] == "certified"
    result = verify_critical_point(problem, certificate, backend=backend)
    assert result["valid"]
    low, high = map(float.fromhex, result["critical_value_hex"])
    assert low <= 2 <= high


def test_nonroot_box_is_not_accepted():
    problem = CriticalPointProblem((X-1)**2 + Y*Y + Z*Z, VARIABLES, 5.)
    assert certify_critical_point(problem, box((0., 0., 0.)))["status"] == "unknown"


def test_degenerate_critical_point_does_not_get_unique_root_claim():
    problem = CriticalPointProblem(X**4 + Y*Y + Z*Z, VARIABLES, 5.)
    assert certify_critical_point(problem, box((0., 0., 0.)))["status"] == "unknown"


def test_entire_sphere_critical_is_not_an_isolated_wall_root():
    problem = CriticalPointProblem(X*X + Y*Y + Z*Z, VARIABLES, 2., "wall")
    assert certify_critical_point(problem, box((0., 0., 2., 2.)))["status"] == "unknown"


def test_outside_ball_root_is_rejected():
    problem = CriticalPointProblem((X-3)**2 + Y*Y + Z*Z, VARIABLES, 2.)
    assert certify_critical_point(problem, box((3., 0., 0.)))["status"] == "unknown"


@pytest.mark.parametrize("change", ["matrix", "center", "box", "field", "scope", "status"])
def test_tampered_local_certificate(change):
    problem = CriticalPointProblem(X*X + Y*Y + Z*Z, VARIABLES, 5.)
    certificate = deepcopy(certify_critical_point(problem, box((0., 0., 0.))))
    if change == "matrix":
        certificate["preconditioner_hex"] = [[0.0.hex()] * 3 for _ in range(3)]
    elif change == "center":
        certificate["center_hex"][0] = 2.0.hex()
    elif change == "box":
        certificate["box_hex"][0] = [1.0.hex(), 2.0.hex()]
    elif change == "field":
        certificate["problem"]["field"] = "forged"
    elif change == "scope":
        certificate["source_topology_change_certified"] = True
    else:
        certificate["status"] = "unknown"
    assert not verify_critical_point(problem, certificate)["valid"]


def test_stored_success_flags_are_not_the_verification():
    problem = CriticalPointProblem(X*X + Y*Y + Z*Z, VARIABLES, 5.)
    certificate = certify_critical_point(problem, box((0., 0., 0.)))
    certificate.pop("primary_replay")
    certificate.pop("rational_replay")
    assert verify_critical_point(problem, certificate)["valid"]


@pytest.mark.parametrize("bounds", [[], [(0, 0)] * 3, [(1, -1)] * 3, [(0, float("inf"))] * 3])
def test_invalid_root_boxes(bounds):
    problem = CriticalPointProblem(X*X + Y*Y + Z*Z, VARIABLES, 5.)
    with pytest.raises(ValueError):
        certify_critical_point(problem, bounds)
