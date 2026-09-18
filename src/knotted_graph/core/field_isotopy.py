"""Replayable regular-level continuation certificates for analytic solids.

Experimental: this is a sufficient test, not a complete isotopy classifier.
It compares sublevel sets in a fixed ball or box, never their extracted graphs.
Supported expressions use rational constants, pi, +, *, nonnegative integer
powers, sin and cos. No sampled-gradient or equal-invariant shortcut is used.

Basic arithmetic is outward-rounded IEEE-754 binary64. Transcendental ranges
come from a private 80-bit mpmath interval context, with outward conversion.
A certificate is a full, replayable binary partition of every required stratum.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
import math
from numbers import Integral
import sys
from time import perf_counter
from typing import Any, Iterable

import sympy as sp
from mpmath.ctx_iv import MPIntervalContext

__all__ = ["FieldProblem", "certify", "verify", "tpms_problem"]
Interval = tuple[float, float]
Box = tuple[Interval, ...]
_FULL = (-math.inf, math.inf)
_ZERO = (0.0, 0.0)
_ONE = (1.0, 1.0)


def _down(x: float) -> float:
    return math.nextafter(x, -math.inf)


def _up(x: float) -> float:
    return math.nextafter(x, math.inf)


def _add(a: Interval, b: Interval) -> Interval:
    if a == _ZERO:
        return b
    if b == _ZERO:
        return a
    lo, hi = a[0] + b[0], a[1] + b[1]
    if not math.isfinite(lo) or not math.isfinite(hi):
        return _FULL
    return _down(lo), _up(hi)


def _neg(a: Interval) -> Interval:
    return -a[1], -a[0]


def _sub(a: Interval, b: Interval) -> Interval:
    return _add(a, _neg(b))


def _mul(a: Interval, b: Interval) -> Interval:
    if a == _ZERO or b == _ZERO:
        return _ZERO
    if a == _ONE:
        return b
    if b == _ONE:
        return a
    values = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    if not all(math.isfinite(v) for v in values):
        return _FULL
    return _down(min(values)), _up(max(values))


def _square(a: Interval) -> Interval:
    values = (a[0]*a[0], a[1]*a[1])
    if not all(math.isfinite(v) for v in values):
        return _FULL
    lo = 0.0 if a[0] <= 0.0 <= a[1] else max(0.0, _down(min(values)))
    return lo, _up(max(values))


def _power(a: Interval, n: int) -> Interval:
    if n == 0:
        return _ONE
    if n == 1:
        return a
    if n == 2:
        return _square(a)
    result = _ONE
    while n:
        if n & 1:
            result = _mul(result, a)
        n //= 2
        if n:
            a = _square(a)
    return result


def _rational(p: int, q: int = 1) -> Interval:
    value = Fraction(p, q)
    rounded = float(value)
    if not math.isfinite(rounded):
        raise ValueError("constant is outside the finite binary64 range")
    exact = Fraction.from_float(rounded)
    return (rounded if exact <= value else _down(rounded),
            rounded if exact >= value else _up(rounded))


def _nonzero(a: Interval) -> bool:
    return a[0] > 0.0 or a[1] < 0.0


def _valid_bound(bound: Iterable[float]) -> Interval:
    values = tuple(bound)
    if len(values) != 2:
        raise ValueError("each bound must have two endpoints")
    a, b = (float(v) for v in values)
    if not math.isfinite(a) or not math.isfinite(b) or a > b:
        raise ValueError("bounds must be ordered and finite")
    return a, b


def _canonical(data: Any) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class FieldProblem:
    """An analytic f(x,y,z,parameters) <= 0 in a fixed compact domain.

    ``variables[:3]`` are spatial; the remainder are independent parameters.
    Bounds use the exact binary values of the supplied floats. For a ball the
    first three bounds must contain that entire ball, centered at the origin.
    On a box all face/edge/corner restrictions are checked, not just the bulk.
    """
    expression: sp.Expr
    variables: tuple[sp.Symbol, ...]
    bounds: Box
    domain: str = "ball"
    radius: float | None = None

    def __post_init__(self) -> None:
        expression = sp.sympify(self.expression)
        variables = tuple(self.variables)
        bounds = tuple(_valid_bound(b) for b in self.bounds)
        if len(variables) < 3 or len(variables) != len(bounds):
            raise ValueError("need three spatial variables and matching bounds")
        if len(variables) > 10:
            raise ValueError("at most ten variables are supported")
        if any(not isinstance(s, sp.Symbol) for s in variables):
            raise ValueError("variables must be SymPy symbols")
        if len(set(variables)) != len(variables):
            raise ValueError("variables must be distinct")
        if not expression.free_symbols.issubset(set(variables)):
            raise ValueError("expression has undeclared variables")
        if any(lo == hi for lo, hi in bounds[:3]):
            raise ValueError("spatial box must have positive extent")
        if self.domain not in ("ball", "box"):
            raise ValueError("domain must be 'ball' or 'box'")
        radius = None if self.radius is None else float(self.radius)
        if self.domain == "ball":
            if radius is None or not math.isfinite(radius) or radius <= 0:
                raise ValueError("ball radius must be finite and positive")
            if any(lo > -radius or hi < radius for lo, hi in bounds[:3]):
                raise ValueError("spatial bounds must contain the entire ball")
        elif radius is not None:
            raise ValueError("a box does not take a radius")
        object.__setattr__(self, "expression", expression)
        object.__setattr__(self, "variables", variables)
        object.__setattr__(self, "bounds", bounds)
        object.__setattr__(self, "radius", radius)

    def specification(self) -> dict[str, Any]:
        return {"expression": sp.srepr(self.expression),
                "variables": [sp.srepr(v) for v in self.variables],
                "bounds_hex": [[lo.hex(), hi.hex()] for lo, hi in self.bounds],
                "domain": self.domain,
                "radius_hex": None if self.radius is None else self.radius.hex()}

    @property
    def fingerprint(self) -> str:
        return sha256(_canonical(self.specification()).encode()).hexdigest()


class _Evaluator:
    """A shared expression DAG; never evaluates generated Python source."""
    def __init__(self, problem: FieldProblem):
        if sys.float_info.radix != 2 or sys.float_info.mant_dig != 53:
            raise RuntimeError("the interval backend requires binary64 floats")
        self.iv = MPIntervalContext()
        self.iv.prec = 80
        self.nodes: list[tuple[str, Any]] = []
        self.known: dict[sp.Expr, int] = {}
        self.symbols = {v: i for i, v in enumerate(problem.variables)}
        self.outputs = [self._compile(problem.expression)]
        self.outputs += [self._compile(sp.diff(problem.expression, v))
                         for v in problem.variables[:3]]
        self.trig = lru_cache(maxsize=65536)(self._trig)

    def _enclosure(self, value: Any) -> Interval:
        return _down(float(value.a)), _up(float(value.b))

    def _compile(self, expr: sp.Expr) -> int:
        if expr in self.known:
            return self.known[expr]
        if expr in self.symbols:
            operation = ("var", self.symbols[expr])
        elif isinstance(expr, (sp.Rational, sp.Float)):
            value = sp.Rational(expr)
            operation = ("constant", _rational(int(value.p), int(value.q)))
        elif expr == sp.pi:
            operation = ("constant", self._enclosure(self.iv.pi))
        elif expr.func in (sp.Add, sp.Mul):
            operation = ("add" if expr.func == sp.Add else "mul",
                         tuple(self._compile(a) for a in expr.args))
        elif expr.func == sp.Pow and expr.args[1].is_Integer and expr.args[1] >= 0:
            if expr.args[1] > 256:
                raise ValueError("integer powers above 256 are unsupported")
            operation = ("pow", (self._compile(expr.args[0]), int(expr.args[1])))
        elif expr.func in (sp.sin, sp.cos):
            operation = ("sin" if expr.func == sp.sin else "cos",
                         self._compile(expr.args[0]))
        else:
            raise ValueError(f"unsupported analytic operation: {expr.func}")
        index = len(self.nodes)
        self.nodes.append(operation)
        self.known[expr] = index
        return index

    def _trig(self, kind: str, lo: float, hi: float) -> Interval:
        if not math.isfinite(lo) or not math.isfinite(hi):
            return -1.0, 1.0
        interval = self.iv.mpf([lo, hi])
        value = self.iv.sin(interval) if kind == "sin" else self.iv.cos(interval)
        a, b = self._enclosure(value)
        return max(-1.0, a), min(1.0, b)

    def __call__(self, box: Box) -> tuple[Interval, tuple[Interval, ...]]:
        values = []
        for kind, data in self.nodes:
            if kind == "var":
                result = box[data]
            elif kind == "constant":
                result = data
            elif kind == "add":
                result = _ZERO
                for index in data:
                    result = _add(result, values[index])
            elif kind == "mul":
                result = _ONE
                for index in data:
                    result = _mul(result, values[index])
            elif kind == "pow":
                result = _power(values[data[0]], data[1])
            else:
                result = self.trig(kind, *values[data])
            values.append(result)
        return values[self.outputs[0]], tuple(values[i] for i in self.outputs[1:])


def _strata(problem: FieldProblem) -> list[tuple[str, Box, tuple[int, ...]]]:
    if problem.domain == "ball":
        return [("ball_with_wall", problem.bounds, (0, 1, 2))]
    result = []
    # 0 means free, -/+ mean fixed at the corresponding box endpoint.
    for roles in product((0, -1, 1), repeat=3):
        bounds = list(problem.bounds)
        for axis, role in enumerate(roles):
            if role:
                value = bounds[axis][0 if role == -1 else 1]
                bounds[axis] = (value, value)
        result.append(("".join("0" if r == 0 else ("-" if r == -1 else "+")
                               for r in roles), tuple(bounds),
                       tuple(i for i, r in enumerate(roles) if r == 0)))
    return result


def _criterion(problem: FieldProblem, evaluate: _Evaluator, box: Box,
               free: tuple[int, ...]) -> str | None:
    phi = None
    if problem.domain == "ball":
        norm2 = _ZERO
        for interval in box[:3]:
            norm2 = _add(norm2, _square(interval))
        phi = _sub(norm2, _square((problem.radius, problem.radius)))
        if phi[0] > 0:
            return "outside_ball"
    value, gradient = evaluate(box)
    if _nonzero(value):
        return "no_level"
    if problem.domain == "box":
        return "stratum_regular" if any(_nonzero(gradient[i]) for i in free) else None
    if phi[1] < 0:
        return "bulk_regular" if any(_nonzero(g) for g in gradient) else None
    # This box may meet the spherical wall. A bulk-gradient test is not enough.
    cross = tuple(_sub(_mul(box[j], gradient[k]), _mul(box[k], gradient[j]))
                  for j, k in ((1, 2), (2, 0), (0, 1)))
    return "wall_transverse" if any(_nonzero(c) for c in cross) else None


def _split(box: Box, axis: int) -> tuple[Box, Box] | None:
    lo, hi = box[axis]
    mid = lo / 2.0 + hi / 2.0
    if not lo < mid < hi:
        return None
    left, right = list(box), list(box)
    left[axis], right[axis] = (lo, mid), (mid, hi)
    return tuple(left), tuple(right)


def _axis(box: Box) -> int | None:
    ordered = sorted(range(len(box)), key=lambda i: box[i][1]-box[i][0], reverse=True)
    return next((i for i in ordered if _split(box, i) is not None), None)


def _positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def certify(problem: FieldProblem, *, max_boxes: int = 250000,
            max_depth: int = 72) -> dict[str, Any]:
    """Cover every stratum and the whole parameter rectangle, or return unknown.

    No negative equivalence verdict follows from exhaustion. A successful result
    supplies the regularity hypotheses for every parameter path in the rectangle.
    A failed path test does not mean its endpoints are inequivalent.
    """
    max_boxes = _positive_int(max_boxes, "max_boxes")
    max_depth = _positive_int(max_depth, "max_depth")
    evaluate = _Evaluator(problem)
    started = perf_counter()
    certificate: dict[str, Any] = {
        "schema": "knottedgraph.field_isotopy.v1", "status": "unknown",
        "problem": problem.specification(), "problem_sha256": problem.fingerprint,
        "arithmetic": "outward_binary64_with_mpmath_iv_80bit_trigonometry",
        "max_boxes": max_boxes, "max_depth": max_depth, "trees": [],
        "counts": {}, "checked_boxes": 0,
        "scope": "analytic_sublevel_ambient_isotopy_in_fixed_domain",
        "periodic_identification": False, "graph_spine_certified": False,
    }
    for label, root, free in _strata(problem):
        stack = [(root, 0)]
        tree = []
        while stack:
            box, depth = stack.pop()
            if certificate["checked_boxes"] >= max_boxes:
                reason = "box_budget"
            else:
                certificate["checked_boxes"] += 1
                accepted = _criterion(problem, evaluate, box, free)
                if accepted is not None:
                    tree.append(".")
                    counts = certificate["counts"]
                    counts[accepted] = counts.get(accepted, 0) + 1
                    continue
                axis = _axis(box)
                if depth < max_depth and axis is not None:
                    tree.append(str(axis))
                    left, right = _split(box, axis)
                    stack.append((right, depth + 1))
                    stack.append((left, depth + 1))
                    continue
                reason = "depth_or_precision_limit"
            certificate["reason"] = reason
            certificate["unresolved"] = {
                "stratum": label, "box_hex": [[a.hex(), b.hex()] for a, b in box],
                "depth": depth, "remaining_stack_boxes": len(stack),
            }
            certificate["elapsed_seconds"] = perf_counter() - started
            # Partial trees are intentionally not mistaken for complete covers.
            certificate["partial_tree"] = "".join(tree) + "?"
            return certificate
        certificate["trees"].append({"stratum": label, "preorder": "".join(tree)})
    certificate["status"] = "certified"
    certificate["reason"] = "all_bulk_and_boundary_strata_regular"
    certificate["elapsed_seconds"] = perf_counter() - started
    certificate["cover_sha256"] = sha256(_canonical(certificate["trees"]).encode()).hexdigest()
    return certificate


def verify(problem: FieldProblem, certificate: dict[str, Any], *,
           max_nodes: int = 2000000) -> dict[str, Any]:
    """Independently replay a complete cover; never trust stored leaf verdicts.

    The caller supplies its own trusted FieldProblem. Serialized symbolic text
    is hashed for matching only; it is never eval'd or sympified by the verifier.
    """
    max_nodes = _positive_int(max_nodes, "max_nodes")
    if not isinstance(certificate, dict):
        return {"valid": False, "reason": "not_a_certificate"}
    if (certificate.get("schema") != "knottedgraph.field_isotopy.v1"
            or certificate.get("status") != "certified"
            or certificate.get("problem_sha256") != problem.fingerprint
            or certificate.get("problem") != problem.specification()):
        return {"valid": False, "reason": "status_or_problem_mismatch"}
    strata = _strata(problem)
    trees = certificate.get("trees")
    if not isinstance(trees, list) or len(trees) != len(strata):
        return {"valid": False, "reason": "missing_strata"}
    evaluate = _Evaluator(problem)
    nodes = 0
    counts: dict[str, int] = {}
    for record, (label, root, free) in zip(trees, strata):
        if (not isinstance(record, dict) or record.get("stratum") != label
                or not isinstance(record.get("preorder"), str)):
            return {"valid": False, "reason": "stratum_format"}
        stack = [root]
        for token in record["preorder"]:
            nodes += 1
            if nodes > max_nodes:
                return {"valid": False, "reason": "verification_budget"}
            if not stack:
                return {"valid": False, "reason": "trailing_tree_data"}
            box = stack.pop()
            if token == ".":
                accepted = _criterion(problem, evaluate, box, free)
                if accepted is None:
                    return {"valid": False, "reason": "invalid_leaf", "stratum": label}
                counts[accepted] = counts.get(accepted, 0) + 1
            elif token in "0123456789" and int(token) < len(box):
                children = _split(box, int(token))
                if children is None:
                    return {"valid": False, "reason": "invalid_split"}
                left, right = children
                stack.extend((right, left))
            else:
                return {"valid": False, "reason": "invalid_tree_token"}
        if stack:
            return {"valid": False, "reason": "incomplete_cover"}
    return {"valid": True, "reason": "cover_recomputed", "nodes": nodes, "counts": counts}


def tpms_problem(family: str, lambda_bounds: Interval, c_bounds: Interval,
                 *, radius: float = 0.72 * (2.25 * math.pi)) -> FieldProblem:
    """Exact analytic family used by _tpms.py at the audited revision.

    Radius defaults to the same binary64 constant, not the printed 5.08938.
    This is the continuous field inside that ball, not a voxel interpolation.
    """
    x, y, z, lam, c = sp.symbols("x y z lam c", real=True)
    fields = {
        "gyroid": sp.sin(x)*sp.cos(y) + sp.sin(y)*sp.cos(z) + sp.sin(z)*sp.cos(x),
        "schwarz_p": sp.cos(x) + sp.cos(y) + sp.cos(z),
        "diamond": sp.cos(x)*sp.cos(y)*sp.cos(z) - sp.sin(x)*sp.sin(y)*sp.sin(z),
    }
    families = {"schwarz_p_to_diamond": ("schwarz_p", "diamond"),
                "gyroid_to_schwarz_p": ("gyroid", "schwarz_p"),
                "gyroid_to_diamond": ("gyroid", "diamond")}
    if family not in families:
        raise ValueError("unknown TPMS family")
    lb, cb = _valid_bound(lambda_bounds), _valid_bound(c_bounds)
    if lb[0] < 0 or lb[1] > 1:
        raise ValueError("TPMS lambda must lie in [0,1]")
    first, last = families[family]
    expr = (1-lam)*fields[first] + lam*fields[last] - c
    return FieldProblem(expr, (x, y, z, lam, c),
                        ((-radius, radius),)*3 + (lb, cb), "ball", radius)
