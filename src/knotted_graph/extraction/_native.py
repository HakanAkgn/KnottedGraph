"""Native acceleration for exact 3-D skeleton extraction.

The compiled backend reproduces the existing sparse 26-neighbour adjacency,
multi-scale junction-zone tracing and embedded edge geometry.  Persistence
selection itself remains in Python/NetworkX so the scientific selection rule and
isomorphism checks are unchanged.  The Python implementation is retained as the
reference/fallback when the extension is unavailable.
"""

from __future__ import annotations

import networkx as nx
import numpy as np

from ._tracing import reduced_weighted, remove_leaves_for_diagnostic
from ._topology_optimized import _anomaly_from_reduced, _core_fingerprint

_NATIVE_IMPORT_ERROR: Exception | None = None
try:
    from . import _skeleton_native
except Exception as exc:  # pragma: no cover - compiler/platform fallback
    _skeleton_native = None
    _NATIVE_IMPORT_ERROR = exc


def native_skeleton_available() -> bool:
    return _skeleton_native is not None


def native_skeleton_import_error() -> Exception | None:
    return _NATIVE_IMPORT_ERROR


def sparse_adjacency_native(image: np.ndarray) -> tuple[np.ndarray, list[list[int]]]:
    if _skeleton_native is None:
        raise RuntimeError("native skeleton backend is unavailable")
    coords, adjacency = _skeleton_native.sparse_adjacency(
        np.asarray(image, dtype=bool, order="C")
    )
    return (
        np.asarray(coords, dtype=np.intp).reshape((-1, 3)),
        [list(map(int, row)) for row in adjacency],
    )


def _compact_graph(candidate: dict) -> nx.MultiGraph:
    graph = nx.MultiGraph()
    graph.add_nodes_from(int(node) for node in candidate["nodes"])
    for u, v, weight, point_count in candidate["edges"]:
        graph.add_edge(
            int(u),
            int(v),
            weight=float(weight),
            _point_count=int(point_count),
        )
    return graph


def _compact_summary(
    candidate: dict,
    *,
    max_degree: int | None,
    anomaly_ratio: float,
) -> tuple[nx.MultiGraph, bool, tuple, bool]:
    graph = _compact_graph(candidate)
    if max_degree is None:
        diagnostic_graph = nx.MultiGraph(graph)
    else:
        diagnostic_graph = remove_leaves_for_diagnostic(graph)

    reduced = reduced_weighted(diagnostic_graph)
    max_observed_degree = max(
        (degree for _, degree in reduced.degree()),
        default=0,
    )
    anomaly_count = _anomaly_from_reduced(reduced, anomaly_ratio)
    valence_ok = max_degree is None or max_observed_degree <= max_degree
    clean = (
        graph.number_of_nodes() > 0
        and valence_ok
        and anomaly_count == 0
        and bool(candidate["geometry_safe"])
    )

    if max_degree is None:
        one_hop_safe = True
    else:
        survivors = set(diagnostic_graph.nodes())
        one_hop_safe = not any(
            u != v
            and u in survivors
            and v in survivors
            and int(data.get("_point_count", 0)) < 5
            for u, v, data in graph.edges(data=True)
        )

    return reduced, clean, _core_fingerprint(reduced), one_hop_safe


def _select_hops(
    candidates: list[dict],
    *,
    max_degree: int | None,
    anomaly_ratio: float,
) -> int:
    summaries: list[tuple[int, tuple[nx.MultiGraph, bool, tuple, bool]]] = []
    clean_candidates: list[
        tuple[int, tuple[nx.MultiGraph, bool, tuple, bool]]
    ] = []

    for candidate in candidates:
        hops = int(candidate["hops"])
        summary = _compact_summary(
            candidate,
            max_degree=max_degree,
            anomaly_ratio=anomaly_ratio,
        )
        summaries.append((hops, summary))
        if summary[1] and summary[3]:
            clean_candidates.append((hops, summary))

    for (left_hops, left), (_right_hops, right) in zip(
        summaries[:-1],
        summaries[1:],
        strict=True,
    ):
        if (
            left[1]
            and left[3]
            and right[1]
            and right[3]
            and left[2] == right[2]
            and nx.is_isomorphic(left[0], right[0])
        ):
            return left_hops

    classes: list[
        list[tuple[int, tuple[nx.MultiGraph, bool, tuple, bool]]]
    ] = []
    for candidate in clean_candidates:
        _hops, summary = candidate
        for group in classes:
            representative = group[0][1]
            if (
                representative[2] == summary[2]
                and nx.is_isomorphic(representative[0], summary[0])
            ):
                group.append(candidate)
                break
        else:
            classes.append([candidate])

    eligible = [group for group in classes if len(group) >= 2]
    if eligible:
        best_group = min(
            eligible,
            key=lambda group: (-len(group), min(item[0] for item in group)),
        )
        return min(item[0] for item in best_group)

    return 0


def native_persistent_extract(
    image: np.ndarray,
    *,
    max_degree: int | None,
    max_hops: int,
    anomaly_ratio: float,
) -> nx.MultiGraph:
    if _skeleton_native is None:
        raise RuntimeError("native skeleton backend is unavailable")

    image = np.asarray(image, dtype=bool, order="C")
    candidates = list(_skeleton_native.plan_candidates(image, int(max_hops)))
    selected_hops = _select_hops(
        candidates,
        max_degree=max_degree,
        anomaly_ratio=anomaly_ratio,
    )
    materialized = _skeleton_native.materialize_candidate(
        image,
        int(selected_hops),
    )

    graph = nx.MultiGraph()
    for node, position in materialized["nodes"]:
        graph.add_node(
            int(node),
            pos=np.asarray(position, dtype=float),
        )
    for u, v, points, weight in materialized["edges"]:
        graph.add_edge(
            int(u),
            int(v),
            pts=np.asarray(points, dtype=float),
            weight=float(weight),
        )
    return graph
