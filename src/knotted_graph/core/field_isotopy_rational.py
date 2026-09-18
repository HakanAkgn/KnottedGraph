"""Independent rational Taylor trigonometry for replaying field certificates.

No mpmath transcendental values are used by this backend. Endpoint Taylor
remainders and the Machin pi enclosure use Fraction arithmetic. The enclosing
DAG arithmetic remains the core outward-rounded binary64 implementation.
Arguments outside [-8,8] conservatively return [-1,1], never a guessed range.
This is an independent trigonometric check, not formal verification of Python.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from typing import Any

import sympy as sp

from .field_isotopy import (
    FieldProblem, _Evaluator, _criterion, _positive_int, _rational, _split, _strata,
)

__all__ = ["verify_rational"]


@lru_cache(maxsize=1)
def _pi_bounds() -> tuple[Fraction, Fraction]:
    """pi = 16 atan(1/5) - 4 atan(1/239), with alternating-series bounds."""
    def atan(q: int) -> tuple[Fraction, Fraction]:
        total = sum((Fraction((-1)**k, (2*k+1)*q**(2*k+1)) for k in range(64)), Fraction(0))
        error = Fraction(1, 129*q**129)
        return total-error, total+error
    a, b = atan(5), atan(239)
    return 16*a[0]-4*b[1], 16*a[1]-4*b[0]


@lru_cache(maxsize=65536)
def _endpoint(kind: str, value: float) -> tuple[Fraction, Fraction]:
    x = Fraction(value)
    square = x*x
    term = x if kind == "sin" else Fraction(1)
    total = term
    for n in range(1, 48):
        divisor = (2*n)*(2*n+1) if kind == "sin" else (2*n-1)*(2*n)
        term *= -square / divisor
        total += term
    divisor = 96*97 if kind == "sin" else 95*96
    error = abs(term)*square/divisor
    return total-error, total+error


def _trig_interval(kind: str, lo: float, hi: float) -> tuple[float, float]:
    if not math.isfinite(lo) or not math.isfinite(hi) or lo < -8 or hi > 8:
        return -1., 1.
    ends = (_endpoint(kind, lo), _endpoint(kind, hi))
    lower, upper = min(e[0] for e in ends), max(e[1] for e in ends)
    pi_lo, pi_hi = _pi_bounds()
    # All stationary points in [-8,8] occur in this finite set (pi > 3).
    for k in range(-4, 5):
        for positive in (False, True):
            coefficient = (Fraction(2*k) + (Fraction(1, 2) if positive else Fraction(-1, 2))
                           if kind == "sin" else Fraction(2*k + (0 if positive else 1)))
            a, b = sorted((coefficient*pi_lo, coefficient*pi_hi))
            if Fraction(lo) <= b and a <= Fraction(hi):
                if positive:
                    upper = max(upper, Fraction(1))
                else:
                    lower = min(lower, Fraction(-1))
    lower, upper = max(Fraction(-1), lower), min(Fraction(1), upper)
    return _rational(lower.numerator, lower.denominator)[0], _rational(upper.numerator, upper.denominator)[1]


class _RationalEvaluator(_Evaluator):
    def _compile(self, expr: sp.Expr) -> int:
        if expr == sp.pi and expr not in self.known:
            a, b = _pi_bounds()
            value = (_rational(a.numerator, a.denominator)[0],
                     _rational(b.numerator, b.denominator)[1])
            index = len(self.nodes)
            self.nodes.append(("constant", value))
            self.known[expr] = index
            return index
        return super()._compile(expr)

    def _trig(self, kind: str, lo: float, hi: float) -> tuple[float, float]:
        return _trig_interval(kind, lo, hi)


def verify_rational(problem: FieldProblem, certificate: dict[str, Any], *,
                    max_nodes: int = 2000000) -> dict[str, Any]:
    """Reconstruct all strata and replay using rational trigonometric bounds."""
    max_nodes = _positive_int(max_nodes, "max_nodes")
    if (not isinstance(certificate, dict)
            or certificate.get("schema") != "knottedgraph.field_isotopy.v1"
            or certificate.get("status") != "certified"
            or certificate.get("problem_sha256") != problem.fingerprint
            or certificate.get("problem") != problem.specification()):
        return {"valid": False, "reason": "status_or_problem_mismatch"}
    strata = _strata(problem)
    trees = certificate.get("trees")
    if not isinstance(trees, list) or len(trees) != len(strata):
        return {"valid": False, "reason": "missing_strata"}
    evaluate = _RationalEvaluator(problem)
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
                    return {"valid": False, "reason": "leaf_not_verified_by_rational_backend", "stratum": label}
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
    return {"valid": True, "reason": "cover_recomputed_with_rational_Taylor_trigonometry",
            "nodes": nodes, "counts": counts, "mpmath_transcendentals_used": False}
