"""Finite-window Bloch scans with explicit reconstruction and evaluation scope.

Polynomial signatures are not source-solid equivalence classes. Failed or
missing evaluations have the reserved label -1, never a measured phase label.
The default scan uses a replayed cubical retraction without geometric smoothing.
No periodic face identification or analytic-to-voxel certificate is inferred.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from numbers import Integral
from typing import Callable, Sequence

import numpy as np
import sympy as sp

BlochVector = tuple[sp.Expr, sp.Expr, sp.Expr]
BlochVectorFactory = Callable[[float], Sequence[sp.Expr]]


def _vector(value: Sequence[sp.Expr]) -> BlochVector:
    if len(value) != 3:
        raise ValueError("a Bloch vector must contain exactly three components")
    return tuple(sp.sympify(component) for component in value)


def _axis(values, name):
    axis = np.asarray(values, dtype=float)
    if axis.ndim != 1 or not len(axis) or not np.isfinite(axis).all():
        raise ValueError(f"{name} must be non-empty, finite and one-dimensional")
    if len(np.unique(axis)) != len(axis):
        raise ValueError(f"{name} must contain distinct samples")
    if len(axis) > 1 and np.any(np.diff(axis) <= 0):
        raise ValueError(f"{name} must be strictly increasing")
    return axis.copy()


@dataclass(frozen=True)
class NodalBlochPath:
    start: BlochVectorFactory
    end: BlochVectorFactory
    start_name: str = "start"
    end_name: str = "end"

    def endpoints(self, gamma: float) -> tuple[BlochVector, BlochVector]:
        gamma = float(gamma)
        if not np.isfinite(gamma):
            raise ValueError("gamma must be finite")
        return _vector(self.start(gamma)), _vector(self.end(gamma))

    def at(self, gamma: float, lam: float) -> BlochVector:
        lam = float(lam)
        if not np.isfinite(lam) or not 0 <= lam <= 1:
            raise ValueError("lam must lie in [0, 1]")
        start, end = self.endpoints(gamma)
        return tuple(sp.expand((1 - lam) * a + lam * b) for a, b in zip(start, end))

    def at_components(self, gamma: float, weights: Sequence[float]) -> BlochVector:
        if len(weights) != 3:
            raise ValueError("weights must contain (lambda_x, lambda_y, lambda_z)")
        weights = tuple(float(value) for value in weights)
        if any(not np.isfinite(value) or not 0 <= value <= 1 for value in weights):
            raise ValueError("all component weights must be finite and lie in [0, 1]")
        start, end = self.endpoints(gamma)
        return tuple(sp.expand((1 - w) * a + w * b) for a, b, w in zip(start, end, weights))


@dataclass(frozen=True)
class NodalPhaseRecord:
    lam: float
    gamma: float
    yamada: sp.Expr | None
    phase_signature: str
    error: str | None = None
    evaluation_kind: str = "legacy-unscoped"
    normalization: str | None = None
    is_subcubic: bool | None = None
    projection: dict | None = None
    reconstruction: dict = field(default_factory=dict)

    @property
    def available(self) -> bool:
        return (self.error is None and self.yamada is not None
                and not self.phase_signature.startswith(("error:", "unavailable:")))


@dataclass
class NodalPhaseScanResult:
    lambdas: np.ndarray
    gammas: np.ndarray
    records: list[NodalPhaseRecord]

    def record_grid(self) -> np.ndarray:
        lambdas, gammas = _axis(self.lambdas, "lambdas"), _axis(self.gammas, "gammas")
        lookup = {}
        wanted = {(float(g), float(l)) for g in gammas for l in lambdas}
        for record in self.records:
            key = (float(record.gamma), float(record.lam))
            if key in lookup or key not in wanted:
                raise ValueError("duplicate or off-grid phase record")
            lookup[key] = record
        if set(lookup) != wanted:
            raise ValueError("phase records do not cover the requested grid")
        return np.asarray([[lookup[float(g), float(l)] for l in lambdas] for g in gammas], dtype=object)

    def phase_grid(self) -> tuple[np.ndarray, dict[int, str]]:
        """Operational polynomial labels; -1 is unavailable, not a phase."""
        grid = self.record_grid()
        signatures = sorted({r.phase_signature for r in grid.flat if r.available})
        ids = {signature: i for i, signature in enumerate(signatures)}
        labels = np.full(grid.shape, -1, dtype=int)
        for index, record in np.ndenumerate(grid):
            if record.available:
                labels[index] = ids[record.phase_signature]
        return labels, {i: signature for signature, i in ids.items()}

    def availability_grid(self) -> np.ndarray:
        return np.asarray([[r.available for r in row] for row in self.record_grid()], dtype=bool)

    def transition_intervals(self) -> list[dict]:
        """Report successful signature changes, never source-topology claims.

        An unavailable cell breaks an interval; it is neither bridged nor
        treated as a transition. Values with different evaluation scope are
        not compared to one another.
        """
        grid = self.record_grid()
        changes = []
        for row, gamma in enumerate(self.gammas):
            for column in range(1, len(self.lambdas)):
                left, right = grid[row, column - 1], grid[row, column]
                same_scope = (left.evaluation_kind == right.evaluation_kind
                              and left.normalization == right.normalization
                              and left.is_subcubic == right.is_subcubic)
                if (left.available and right.available and same_scope
                        and left.phase_signature != right.phase_signature):
                    changes.append({
                        "gamma": float(gamma),
                        "lambda_left": float(self.lambdas[column - 1]),
                        "lambda_right": float(self.lambdas[column]),
                        "phase_left": left.phase_signature,
                        "phase_right": right.phase_signature,
                        "evidence": "evaluated_spine_signature_change",
                        "source_topology_transition_certified": False,
                    })
        return changes


def _signature(polynomial: sp.Expr, audit: dict | None = None) -> str:
    polynomial = sp.sympify(polynomial)
    if polynomial.has(sp.nan, sp.zoo, sp.oo, -sp.oo):
        raise ValueError("nonfinite polynomial evaluation")
    canonical = sp.factor(sp.together(sp.expand(polynomial)))
    if audit is None:
        return "yamada:" + sp.srepr(canonical)
    return ":".join((audit["evaluation_kind"], audit["normalization"], sp.srepr(canonical)))


class NodalPhaseScan:
    def __init__(self, path: NodalBlochPath, *, lambdas: Sequence[float],
                 gammas: Sequence[float], dimension: int = 96, span=None,
                 normalize_yamada: bool = True, yamada_variable: sp.Symbol | None = None,
                 yamada_options: dict | None = None, continue_on_error: bool = True,
                 reconstruction: str = "cubical", cubical_options: dict | None = None) -> None:
        self.path = path
        self.lambdas = _axis(lambdas, "lambdas")
        self.gammas = _axis(gammas, "gammas")
        if np.any((self.lambdas < 0) | (self.lambdas > 1)):
            raise ValueError("all lambda samples must lie in [0, 1]")
        if isinstance(dimension, bool) or not isinstance(dimension, Integral) or dimension < 2:
            raise ValueError("dimension must be an integer >= 2")
        if reconstruction not in ("cubical", "guarded"):
            raise ValueError("reconstruction must be cubical or guarded")
        self.dimension, self.span = int(dimension), span
        self.normalize_yamada = bool(normalize_yamada)
        self.yamada_variable = yamada_variable if yamada_variable is not None else sp.Symbol("A")
        self.yamada_options = dict(yamada_options or {})
        self.continue_on_error = bool(continue_on_error)
        self.reconstruction = reconstruction
        self.cubical_options = dict(cubical_options or {})

    def run(self) -> NodalPhaseScanResult:
        from knotted_graph.applications.nodal.skeleton import NodalSkeleton
        from knotted_graph.applications.phase_maps import _compute_yamada_audited

        records = []
        for gamma in self.gammas:
            start, end = self.path.endpoints(float(gamma))
            for lam in self.lambdas:
                vector = tuple(sp.expand((1 - float(lam)) * a + float(lam) * b)
                               for a, b in zip(start, end))
                polynomial, error = None, None
                audit, reconstruction = {}, {}
                try:
                    kwargs = {"char": vector, "dimension": self.dimension}
                    if self.span is not None:
                        kwargs["span"] = self.span
                    model = NodalSkeleton(**kwargs)
                    graph = model.skeleton_graph(
                        smooth_epsilon=0, reconstruction=self.reconstruction,
                        cubical_options=self.cubical_options)
                    reconstruction = dict(graph.graph.get("reconstruction", {}))
                    options = {**self.yamada_options, "normalize": self.normalize_yamada}
                    polynomial, audit = _compute_yamada_audited(graph, self.yamada_variable, options)
                    signature = _signature(polynomial, audit)
                except Exception as exc:
                    if not self.continue_on_error:
                        raise
                    polynomial = None
                    error = f"{type(exc).__name__}: {exc}"
                    signature = "unavailable:" + error
                records.append(NodalPhaseRecord(
                    float(lam), float(gamma), polynomial, signature, error,
                    audit.get("evaluation_kind", "unavailable"), audit.get("normalization"),
                    audit.get("is_subcubic"), audit.get("projection"), reconstruction))
        return NodalPhaseScanResult(self.lambdas.copy(), self.gammas.copy(), records)


__all__ = ["NodalBlochPath", "NodalPhaseRecord", "NodalPhaseScan", "NodalPhaseScanResult"]
