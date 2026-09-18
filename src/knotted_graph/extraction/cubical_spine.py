"""Witnessed deformation retractions of explicitly specified voxel solids.

This opt-in route does not use Lee thinning, junction clustering, endpoint
relocation or RDP. Every deletion is an elementary free-face collapse of a
closed cubical complex. A successful graph endpoint therefore has a checked
retraction from THIS voxel complex, not merely matching Betti numbers.

The certificate does not identify this complex with an analytic level set,
certify a manifold/regular neighborhood, or validate a legacy extracted spine.
Existing cavity handling and all polynomial evaluators are unaffected.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
import heapq
from itertools import product
import json
from numbers import Integral
import struct
from typing import Any

import networkx as nx
import numpy as np

Cell = tuple[int, int, int]
_OFFSETS = tuple(product((-1, 0, 1), repeat=3))
_SCHEMA = "knottedgraph.cubical_retraction.v1"


@dataclass
class CubicalRetraction:
    certificate: dict[str, Any]
    graph: nx.MultiGraph | None


def _positive(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def sample_grid_breaks(samples) -> np.ndarray:
    """Dual voxel boundaries clipped to the inclusive sampled finite interval.

    No padding, periodic identification, or extension beyond either endpoint
    is implicit. Interior boundaries are arithmetic midpoints of samples.
    """
    x = np.asarray(samples, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all() or np.any(np.diff(x) <= 0):
        raise ValueError("samples must be finite and strictly increasing")
    out = np.r_[x[0], x[:-1] / 2 + x[1:] / 2, x[-1]]
    if np.any(np.diff(out) <= 0):
        raise ValueError("midpoints are unresolved at the supplied precision")
    return out


def _input(mask, gridlines):
    image = np.asarray(mask)
    if image.dtype != np.bool_ or image.ndim != 3 or min(image.shape) < 1:
        raise ValueError("mask must be a nonempty-shape three-dimensional boolean array")
    if max(image.shape) >= 2**30:
        raise ValueError("grid exceeds the certificate coordinate range")
    if gridlines is None:
        axes = tuple(np.arange(n + 1, dtype=float) for n in image.shape)
    else:
        if len(gridlines) != 3:
            raise ValueError("provide exactly three arrays of grid boundaries")
        axes = tuple(np.asarray(a, dtype=float) for a in gridlines)
    for a, n in zip(axes, image.shape):
        if a.ndim != 1 or len(a) != n + 1 or not np.isfinite(a).all() or np.any(np.diff(a) <= 0):
            raise ValueError("each grid boundary array must have n+1 strictly increasing finite entries")
    spec = {
        "shape": list(image.shape),
        "occupied_voxels": int(image.sum()),
        "mask_sha256": sha256(image.astype(np.uint8).tobytes(order="C")).hexdigest(),
        "gridlines_hex": [[float(x).hex() for x in a] for a in axes],
        "realization": "union_of_closed_axis_aligned_voxel_cells",
        "periodic_identification": False,
    }
    return image, axes, spec


def _dim(cell: Cell) -> int:
    return sum(x & 1 for x in cell)


def _faces(cell: Cell):
    for axis in range(3):
        if cell[axis] & 1:
            for step in (-1, 1):
                face = list(cell)
                face[axis] += step
                yield tuple(face)


def _cofaces(cell: Cell):
    for axis in range(3):
        if not cell[axis] & 1:
            for step in (-1, 1):
                coface = list(cell)
                coface[axis] += step
                yield tuple(coface)


def _proper_subfaces(cell: Cell):
    choices = [(x - 1, x, x + 1) if x & 1 else (x,) for x in cell]
    for face in product(*choices):
        if face != cell:
            yield face


def _complex(image: np.ndarray, max_cells: int) -> set[Cell]:
    cells: set[Cell] = set()
    for index in np.argwhere(image):
        center = tuple(2 * int(x) + 1 for x in index)
        cells.update(tuple(center[a] + delta[a] for a in range(3)) for delta in _OFFSETS)
        if len(cells) > max_cells:
            raise ValueError("cubical-complex cell budget exceeded; no retraction was produced")
    return cells


def _free_parent(face: Cell, cells: set[Cell]) -> Cell | None:
    if face not in cells:
        return None
    parents = [c for c in _cofaces(face) if c in cells]
    if len(parents) != 1:
        return None
    parent = parents[0]
    if any(c in cells for c in _cofaces(parent)):
        return None
    return parent


def _digest(cells) -> str:
    h = sha256()
    for cell in sorted(cells):
        h.update(struct.pack("<iii", *cell))
    return h.hexdigest()


def _counts(cells) -> list[int]:
    counts = Counter(_dim(c) for c in cells)
    return [counts[d] for d in range(4)]


def _edge_key(a, b):
    return (a, b) if a < b else (b, a)


def _graph(cells: set[Cell], axes) -> nx.MultiGraph:
    """Suppress degree-two vertices without changing any embedded polyline."""
    if any(_dim(c) > 1 for c in cells):
        raise ValueError("terminal complex is not one-dimensional")
    raw = nx.Graph()
    raw.add_nodes_from(sorted(c for c in cells if _dim(c) == 0))
    for cell in sorted(c for c in cells if _dim(c) == 1):
        a, b = tuple(_faces(cell))
        if a not in raw or b not in raw:
            raise ValueError("terminal complex is not closed")
        raw.add_edge(a, b)
    anchors = {n for n, d in raw.degree() if d != 2}
    for component in nx.connected_components(raw):
        if not anchors.intersection(component):
            anchors.add(min(component))
    graph = nx.MultiGraph()

    def position(node):
        return tuple(float(axes[a][node[a] // 2]) for a in range(3))

    for node in sorted(anchors):
        graph.add_node(node, pos=position(node), cubical_vertex=node)
    used = set()
    for start in sorted(anchors):
        for nxt in sorted(raw[start]):
            if _edge_key(start, nxt) in used:
                continue
            points = [start, nxt]
            used.add(_edge_key(start, nxt))
            previous, current = start, nxt
            while current not in anchors:
                candidates = [n for n in raw[current] if n != previous]
                if len(candidates) != 1:
                    raise ValueError("invalid degree-two path")
                nxt = candidates[0]
                edge = _edge_key(current, nxt)
                if edge in used:
                    raise ValueError("path encountered a previously used edge")
                used.add(edge)
                points.append(nxt)
                previous, current = current, nxt
            graph.add_edge(start, current, pts=np.array([position(p) for p in points]),
                           cubical_vertices=points)
    if len(used) != raw.number_of_edges():
        raise ValueError("not every terminal edge was represented")
    return graph


def cubical_retract(mask, *, gridlines=None, max_cells: int = 3000000,
                    max_collapses: int = 3000000) -> CubicalRetraction:
    """Generate a replayable retraction, returning a graph only if reached.

    A higher-dimensional residual or exhausted budget is not a proof that no
    graph spine exists. No topology-changing fallback is ever applied. Zero-
    and one-dimensional free pairs may be removed, but no component is deleted.
    """
    max_cells = _positive(max_cells, "max_cells")
    max_collapses = _positive(max_collapses, "max_collapses")
    image, axes, spec = _input(mask, gridlines)
    cells = _complex(image, max_cells)
    initial_counts = _counts(cells)
    initial_hash = _digest(cells)
    queue = []
    pending = set()

    def enqueue(face):
        if face in pending:
            return
        parent = _free_parent(face, cells)
        if parent is not None:
            heapq.heappush(queue, (-_dim(parent), face))
            pending.add(face)

    for cell in sorted(cells):
        if _dim(cell) < 3:
            enqueue(cell)
    moves = []
    reason = "no_more_free_pairs"
    while queue:
        if len(moves) >= max_collapses:
            reason = "collapse_budget"
            break
        _, face = heapq.heappop(queue)
        pending.remove(face)
        parent = _free_parent(face, cells)
        if parent is None:
            continue
        cells.remove(face)
        cells.remove(parent)
        moves.append([list(face), list(parent)])
        # Faces which become maximal can make their own faces newly free.
        for affected in _proper_subfaces(parent):
            enqueue(affected)
    dimension = max((_dim(c) for c in cells), default=-1)
    kind = "empty" if dimension < 0 else "graph_retract" if dimension <= 1 else "partial_retract"
    certificate = {
        "schema": _SCHEMA, "source": spec, "kind": kind, "stop_reason": reason,
        "initial_cell_counts": initial_counts, "initial_complex_sha256": initial_hash,
        "collapses": moves, "terminal_cells": [list(c) for c in sorted(cells)],
        "terminal_cell_counts": _counts(cells), "terminal_complex_sha256": _digest(cells),
        "terminal_dimension": dimension, "analytic_source_correspondence_certified": False,
        "manifold_or_regular_neighborhood_certified": False,
        "legacy_graph_embedding_validated": False,
    }
    graph = _graph(cells, axes) if kind == "graph_retract" else None
    if graph is not None:
        components = nx.number_connected_components(graph)
        certificate["graph_summary"] = {
            "vertices": graph.number_of_nodes(), "edges": graph.number_of_edges(),
            "components": components,
            "cycle_rank": graph.number_of_edges() - graph.number_of_nodes() + components,
            "max_degree": max((d for _, d in graph.degree()), default=0),
        }
    return CubicalRetraction(certificate, graph)


def verify_cubical_retract(mask, certificate, *, gridlines=None,
                           max_cells: int = 3000000, max_collapses: int = 3000000) -> dict:
    """Rebuild from the caller's mask and replay every free-face condition.

    The input mask is never inferred from the certificate. Recorded summary
    statistics do not substitute for checking the complete sequence and endpoint.
    """
    max_cells = _positive(max_cells, "max_cells")
    max_collapses = _positive(max_collapses, "max_collapses")
    image, axes, spec = _input(mask, gridlines)
    if not isinstance(certificate, dict) or certificate.get("schema") != _SCHEMA:
        return {"valid": False, "reason": "schema"}
    if certificate.get("source") != spec:
        return {"valid": False, "reason": "source_mismatch"}
    cells = _complex(image, max_cells)
    if (certificate.get("initial_complex_sha256") != _digest(cells)
            or certificate.get("initial_cell_counts") != _counts(cells)):
        return {"valid": False, "reason": "initial_complex_mismatch"}
    moves = certificate.get("collapses")
    if not isinstance(moves, list) or len(moves) > max_collapses:
        return {"valid": False, "reason": "collapse_list_or_budget"}
    for index, pair in enumerate(moves):
        if (not isinstance(pair, list) or len(pair) != 2
                or any(not isinstance(c, list) or len(c) != 3
                       or any(type(x) is not int for x in c) for c in pair)):
            return {"valid": False, "reason": "malformed_collapse", "index": index}
        face, parent = (tuple(c) for c in pair)
        # Replay explicit incidence and maximality; no stored verdict is trusted.
        if face not in cells or parent not in cells or face not in set(_faces(parent)):
            return {"valid": False, "reason": "missing_or_nonincident_pair", "index": index}
        if [p for p in _cofaces(face) if p in cells] != [parent]:
            return {"valid": False, "reason": "face_not_free", "index": index}
        if any(p in cells for p in _cofaces(parent)):
            return {"valid": False, "reason": "parent_not_maximal", "index": index}
        cells.difference_update((face, parent))
    dimension = max((_dim(c) for c in cells), default=-1)
    kind = "empty" if dimension < 0 else "graph_retract" if dimension <= 1 else "partial_retract"
    if (certificate.get("terminal_cells") != [list(c) for c in sorted(cells)]
            or certificate.get("terminal_complex_sha256") != _digest(cells)
            or certificate.get("terminal_cell_counts") != _counts(cells)
            or certificate.get("terminal_dimension") != dimension
            or certificate.get("kind") != kind):
        return {"valid": False, "reason": "terminal_mismatch"}
    for claim in ("analytic_source_correspondence_certified",
                  "manifold_or_regular_neighborhood_certified", "legacy_graph_embedding_validated"):
        if certificate.get(claim) is not False:
            return {"valid": False, "reason": "unsupported_claim"}
    summary = None
    if kind == "graph_retract":
        graph = _graph(cells, axes)
        components = nx.number_connected_components(graph)
        summary = {"vertices": graph.number_of_nodes(), "edges": graph.number_of_edges(),
                   "components": components,
                   "cycle_rank": graph.number_of_edges() - graph.number_of_nodes() + components,
                   "max_degree": max((d for _, d in graph.degree()), default=0)}
        if certificate.get("graph_summary") != summary:
            return {"valid": False, "reason": "graph_summary_mismatch"}
    return {"valid": True, "reason": "every_elementary_collapse_replayed", "kind": kind,
            "collapses": len(moves), "terminal_dimension": dimension, "graph_summary": summary,
            "analytic_source_correspondence_certified": False}


def certificate_json(certificate: dict) -> str:
    """Canonical serialization suitable for gzip and SHA-256 manifests."""
    return json.dumps(certificate, sort_keys=True, separators=(",", ":"), allow_nan=False)
