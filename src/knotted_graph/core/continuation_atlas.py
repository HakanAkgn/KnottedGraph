"""Compose independently replayed analytic continuation certificates.

A connected verified cover gives POSITIVE source-solid equivalence evidence.
Disconnected cover components are not asserted to be different topologies.
This module never propagates a chosen spine's polynomial across a cover.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from hashlib import sha256
import json
import math

from .field_isotopy import FieldProblem, verify
from .field_isotopy_rational import verify_rational

__all__ = ["ContinuationAtlas"]


def _hash(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _source(problem):
    specification = problem.specification()
    specification["bounds_hex"] = specification["bounds_hex"][:3]
    return _hash(specification)


def _contains(bounds, point):
    return all(lo <= value <= hi for (lo, hi), value in zip(bounds, point))


def _intersection(left, right):
    low = tuple(max(a[0], b[0]) for a, b in zip(left, right))
    high = tuple(min(a[1], b[1]) for a, b in zip(left, right))
    return low if all(a <= b for a, b in zip(low, high)) else None


@dataclass(frozen=True)
class _Block:
    bounds: tuple
    problem_sha256: str
    certificate_sha256: str


class ContinuationAtlas:
    """Positive ambient-isotopy certificates within one fixed analytic family.

    The trusted FieldProblem is supplied by the caller, never parsed from an
    untrusted serialized expression. Both original and rational-trigonometric
    replay must succeed before a block can enter the equivalence cover.
    """
    def __init__(self):
        self._source_sha256 = None
        self._parameter_count = None
        self._blocks = []
        self._neighbors = []

    def add(self, problem: FieldProblem, certificate: dict, *, max_nodes=2000000) -> int:
        parameter_count = len(problem.variables) - 3
        if parameter_count < 1:
            raise ValueError("a continuation atlas requires at least one parameter")
        source = _source(problem)
        if self._source_sha256 is not None and source != self._source_sha256:
            raise ValueError("different analytic source, variables or spatial domain")
        primary = verify(problem, certificate, max_nodes=max_nodes)
        rational = verify_rational(problem, certificate, max_nodes=max_nodes)
        if not (primary["valid"] and rational["valid"]):
            raise ValueError("a block requires two successful full certificate replays")
        block = _Block(tuple(problem.bounds[3:]), problem.fingerprint, _hash(certificate))
        for index, existing in enumerate(self._blocks):
            if existing == block:
                return index
        index = len(self._blocks)
        neighbors = []
        for other, existing in enumerate(self._blocks):
            if _intersection(existing.bounds, block.bounds) is not None:
                neighbors.append(other)
                self._neighbors[other].append(index)
        self._blocks.append(block)
        self._neighbors.append(neighbors)
        self._source_sha256 = source
        self._parameter_count = parameter_count
        return index

    def _point(self, point):
        point = tuple(float(v) for v in point)
        if self._parameter_count is not None and len(point) != self._parameter_count:
            raise ValueError("parameter point has the wrong dimension")
        if not point or not all(math.isfinite(v) for v in point):
            raise ValueError("parameter point must be non-empty and finite")
        return point

    def compare(self, first, last):
        first, last = self._point(first), self._point(last)
        if len(first) != len(last):
            raise ValueError("parameter points have different dimensions")
        starts = [i for i, block in enumerate(self._blocks) if _contains(block.bounds, first)]
        ends = {i for i, block in enumerate(self._blocks) if _contains(block.bounds, last)}
        parents = {i: None for i in starts}
        queue = deque(starts)
        terminal = None
        while queue:
            index = queue.popleft()
            if index in ends:
                terminal = index
                break
            for neighbor in sorted(self._neighbors[index]):
                if neighbor not in parents:
                    parents[neighbor] = index
                    queue.append(neighbor)
        if terminal is None:
            return {"status": "unknown", "reason": "no_path_in_verified_cover",
                    "inequivalence_established": False}
        chain = []
        current = terminal
        while current is not None:
            chain.append(current)
            current = parents[current]
        chain.reverse()
        points = [first]
        for left, right in zip(chain, chain[1:]):
            points.append(_intersection(self._blocks[left].bounds, self._blocks[right].bounds))
        points.append(last)
        segments = []
        for offset, index in enumerate(chain):
            block = self._blocks[index]
            entry, exit_point = points[offset], points[offset + 1]
            if not (_contains(block.bounds, entry) and _contains(block.bounds, exit_point)):
                raise RuntimeError("invalid internal continuation path")
            segments.append({"block": index, "parameter_entry": entry, "parameter_exit": exit_point,
                             "problem_sha256": block.problem_sha256,
                             "certificate_sha256": block.certificate_sha256})
        return {"status": "equivalent", "reason": "path_in_doubly_replayed_closed_cover",
                "equivalence": "ambient_isotopy_of_specified_analytic_source_solids",
                "source_sha256": self._source_sha256, "segments": segments,
                "selected_graph_spines_validated": False}

    def cover_labels(self, points):
        """Labels of connected certified cover components; -1 is unresolved.

        These are not a count or a complete classification of topology classes.
        """
        components = {}
        for start in range(len(self._blocks)):
            if start in components:
                continue
            pending = [start]
            components[start] = start
            while pending:
                index = pending.pop()
                for neighbor in self._neighbors[index]:
                    if neighbor not in components:
                        components[neighbor] = start
                        pending.append(neighbor)
        labels = []
        for raw in points:
            point = self._point(raw)
            containing = [i for i, block in enumerate(self._blocks) if _contains(block.bounds, point)]
            labels.append(components[containing[0]] if containing else -1)
        return labels

    @property
    def block_count(self):
        return len(self._blocks)
