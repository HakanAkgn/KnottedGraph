"""Replayable local certificates for bulk and spherical-wall critical points.

A numerical root is only a candidate. Acceptance requires a Krawczyk inclusion
strictly inside the proposed box and a rigorously bounded contraction norm < 1.
Both existing interval backends replay the same candidate. No result here says
that all critical points were found or that a critical event changes topology.

For a field h(x,y,z), the bulk system is grad(h)=0. The wall system is
(grad(h)-mu*x, |x|^2-R^2)=0. Constants and supplied box/preconditioner floats
represent their exact binary values. Serialized SymPy text is never evaluated.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math

import numpy as np
import sympy as sp

from .field_isotopy import (
    FieldProblem, _Evaluator, _ZERO, _add, _mul, _neg, _square, _sub, _up,
)
from .field_isotopy_rational import _RationalEvaluator

__all__ = ["CriticalPointProblem", "certify_critical_point", "verify_critical_point"]
_SCHEMA = "knottedgraph.critical_point.v1"


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class CriticalPointProblem:
    field: sp.Expr
    variables: tuple[sp.Symbol, sp.Symbol, sp.Symbol]
    radius: float
    kind: str = "bulk"

    def __post_init__(self):
        variables = tuple(self.variables)
        expression = sp.sympify(self.field)
        radius = float(self.radius)
        if (len(variables) != 3 or len(set(variables)) != 3
                or any(not isinstance(v, sp.Symbol) or v.is_real is not True for v in variables)):
            raise ValueError("three distinct real coordinate symbols are required")
        if (not expression.free_symbols.issubset(set(variables))
                or expression.is_real is not True):
            raise ValueError("field must be real and have no undeclared variables")
        if any(v.name == "_kg_mu" for v in variables):
            raise ValueError("_kg_mu is reserved for the wall multiplier")
        if not math.isfinite(radius) or radius <= 0 or self.kind not in ("bulk", "wall"):
            raise ValueError("positive finite radius and bulk/wall kind required")
        object.__setattr__(self, "variables", variables)
        object.__setattr__(self, "field", expression)
        object.__setattr__(self, "radius", radius)

    @property
    def dimension(self):
        return 3 if self.kind == "bulk" else 4

    def specification(self):
        return {"field": sp.srepr(self.field), "variables": [sp.srepr(v) for v in self.variables],
                "radius_hex": self.radius.hex(), "kind": self.kind}

    @property
    def fingerprint(self):
        return sha256(_json(self.specification()).encode()).hexdigest()

    def system(self):
        gradient = tuple(sp.diff(self.field, v) for v in self.variables)
        if self.kind == "bulk":
            return self.variables, gradient
        multiplier = sp.Symbol("_kg_mu", real=True)
        radius = sp.Rational(self.radius)
        equations = tuple(g - multiplier * v for g, v in zip(gradient, self.variables))
        return self.variables + (multiplier,), equations + (sum(v*v for v in self.variables) - radius*radius,)


def _box(problem, bounds):
    values = np.asarray(bounds, dtype=float)
    if (values.shape != (problem.dimension, 2) or not np.isfinite(values).all()
            or np.any(values[:, 0] >= values[:, 1])):
        raise ValueError("root box must have finite strictly ordered bounds")
    return tuple((float(lo), float(hi)) for lo, hi in values)


class _System:
    def __init__(self, problem, bounds, backend):
        evaluator = _Evaluator if backend == "primary" else _RationalEvaluator
        variables, expressions = problem.system()
        self.problem = problem
        self.evaluators = [evaluator(FieldProblem(f, variables, bounds, domain="box")) for f in expressions]
        self.field = evaluator(FieldProblem(problem.field, variables, bounds, domain="box"))

    def evaluate(self, bounds):
        values, jacobian = [], []
        for index, evaluator in enumerate(self.evaluators):
            value, gradient = evaluator(bounds)
            row = list(gradient)
            if self.problem.kind == "wall":
                row.append(_neg(bounds[index]) if index < 3 else _ZERO)
            values.append(value)
            jacobian.append(row)
        return values, jacobian


def _enclosure(problem, bounds, center, preconditioner, backend):
    system = _System(problem, bounds, backend)
    point = tuple((x, x) for x in center)
    value, _ = system.evaluate(point)
    _, jacobian = system.evaluate(bounds)
    size = problem.dimension
    remainder = []
    for i in range(size):
        row = []
        for j in range(size):
            entry = _ZERO
            for k in range(size):
                coefficient = float(preconditioner[i, k])
                entry = _add(entry, _mul((coefficient, coefficient), jacobian[k][j]))
            identity = (1., 1.) if i == j else _ZERO
            row.append(_sub(identity, entry))
        remainder.append(row)
    norms = []
    for row in remainder:
        norm = 0.
        for lo, hi in row:
            norm = _up(norm + max(abs(lo), abs(hi)))
        norms.append(norm)
    contraction = max(norms)
    if not math.isfinite(contraction) or contraction >= 1:
        return {"valid": False, "reason": "contraction_not_established", "backend": backend}
    differences = [_sub(interval, (x, x)) for interval, x in zip(bounds, center)]
    krawczyk = []
    for i in range(size):
        correction = _ZERO
        for j in range(size):
            c = float(preconditioner[i, j])
            correction = _add(correction, _mul((c, c), value[j]))
        entry = _sub((center[i], center[i]), correction)
        for j in range(size):
            entry = _add(entry, _mul(remainder[i][j], differences[j]))
        krawczyk.append(entry)
    if not all(a < lo <= hi < b for (a, b), (lo, hi) in zip(bounds, krawczyk)):
        return {"valid": False, "reason": "strict_inclusion_not_established", "backend": backend}
    if problem.kind == "bulk":
        norm2 = _ZERO
        for interval in bounds[:3]:
            norm2 = _add(norm2, _square(interval))
        radial = _sub(norm2, _square((problem.radius, problem.radius)))
        if radial[1] >= 0:
            return {"valid": False, "reason": "bulk_box_not_strictly_inside_ball", "backend": backend}
    critical_value, _ = system.field(bounds)
    if not all(math.isfinite(x) for x in critical_value):
        return {"valid": False, "reason": "unbounded_critical_value", "backend": backend}
    return {"valid": True, "reason": "strict_Krawczyk_inclusion_and_contraction", "backend": backend,
            "contraction_norm_upper": contraction,
            "root_enclosure_hex": [[a.hex(), b.hex()] for a, b in krawczyk],
            "critical_value_hex": [x.hex() for x in critical_value],
            "unique_critical_point_in_root_box": True, "source_topology_change_certified": False}


def certify_critical_point(problem, bounds, *, center=None):
    """Attempt a local existence/uniqueness proof, not a global root search."""
    bounds = _box(problem, bounds)
    center = (np.asarray([lo/2 + hi/2 for lo, hi in bounds]) if center is None else np.asarray(center, dtype=float))
    if (center.shape != (problem.dimension,) or not np.isfinite(center).all()
            or any(not lo < x < hi for (lo, hi), x in zip(bounds, center))):
        raise ValueError("center must lie strictly inside the root box")
    center = tuple(map(float, center))
    certificate = {"schema": _SCHEMA, "status": "unknown", "problem": problem.specification(),
                   "problem_sha256": problem.fingerprint,
                   "box_hex": [[lo.hex(), hi.hex()] for lo, hi in bounds],
                   "center_hex": [x.hex() for x in center], "source_topology_change_certified": False}
    try:
        _, jacobian = _System(problem, bounds, "primary").evaluate(tuple((x, x) for x in center))
        midpoint = np.array([[lo/2 + hi/2 for lo, hi in row] for row in jacobian])
        inverse = np.linalg.inv(midpoint)
        if not np.isfinite(inverse).all():
            raise ValueError("nonfinite preconditioner")
    except (np.linalg.LinAlgError, ValueError, OverflowError) as exc:
        certificate["reason"] = f"no_valid_preconditioner: {type(exc).__name__}"
        return certificate
    certificate["preconditioner_hex"] = [[float(x).hex() for x in row] for row in inverse]
    primary = _enclosure(problem, bounds, center, inverse, "primary")
    rational = _enclosure(problem, bounds, center, inverse, "rational")
    certificate["primary_replay"] = primary
    certificate["rational_replay"] = rational
    if primary["valid"] and rational["valid"]:
        certificate["status"] = "certified"
        certificate["reason"] = "both_interval_backends_verified"
    else:
        certificate["reason"] = "local_root_not_certified"
    return certificate


def verify_critical_point(problem, certificate, *, backend="rational"):
    """Recompute every inequality for the caller's own trusted problem."""
    if backend not in ("primary", "rational"):
        raise ValueError("backend must be primary or rational")
    if (not isinstance(certificate, dict) or certificate.get("schema") != _SCHEMA
            or certificate.get("status") != "certified"
            or certificate.get("problem") != problem.specification()
            or certificate.get("problem_sha256") != problem.fingerprint
            or certificate.get("source_topology_change_certified") is not False):
        return {"valid": False, "reason": "status_problem_or_scope_mismatch"}
    try:
        bounds = _box(problem, [[float.fromhex(a), float.fromhex(b)] for a, b in certificate["box_hex"]])
        center = tuple(float.fromhex(x) for x in certificate["center_hex"])
        inverse = np.asarray([[float.fromhex(x) for x in row] for row in certificate["preconditioner_hex"]])
        if (len(center) != problem.dimension or inverse.shape != (problem.dimension, problem.dimension)
                or not np.isfinite(center).all() or not np.isfinite(inverse).all()
                or any(not lo < x < hi for (lo, hi), x in zip(bounds, center))):
            raise ValueError("invalid center or preconditioner")
        return _enclosure(problem, bounds, center, inverse, backend)
    except (KeyError, ValueError, TypeError, OverflowError):
        return {"valid": False, "reason": "malformed_local_certificate"}
