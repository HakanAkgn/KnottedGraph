from __future__ import annotations

import logging
from typing import Any, Sequence

import networkx as nx
import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "EmbeddingValidationError",
    "as_point3",
    "as_polyline",
    "drop_consecutive_duplicates",
    "oriented_edge_polyline",
    "validate_embedding",
    "ensure_embedding",
    "is_embedding",
    "idx_to_coord",
    "get_all_edge_pts",
    "total_edge_pts",
    "smooth_edges",
    "remove_leaf_nodes",
    "simplify_edges",
    "contract_short_edges",
]


class EmbeddingValidationError(ValueError):
    """Raised when a graph cannot satisfy the embedded graph contract."""

    def __init__(self, issues: Sequence[str]):
        self.issues = list(issues)
        super().__init__("; ".join(self.issues))


_FLOAT_REL_TOL = 64.0 * np.finfo(float).eps
_ENDPOINT_REL_TOL = 1e-9


def _point_scale(*points: np.ndarray) -> float:
    """Return a translation-invariant local scale for geometric comparisons."""
    arrays = [np.asarray(point, dtype=float).reshape(-1, 3) for point in points]
    stacked = np.vstack(arrays)
    if len(stacked) <= 1:
        return 1.0
    center = stacked.mean(axis=0)
    scale = float(np.max(np.linalg.norm(stacked - center, axis=1)))
    return scale if np.isfinite(scale) and scale > 0.0 else 1.0


def _points_close(
    left: np.ndarray,
    right: np.ndarray,
    *,
    scale: float,
    rel_tol: float = _ENDPOINT_REL_TOL,
) -> bool:
    """Compare points using local extent rather than absolute coordinate size."""
    distance = float(
        np.linalg.norm(np.asarray(left, dtype=float) - np.asarray(right, dtype=float))
    )
    tolerance = max(_FLOAT_REL_TOL * scale, rel_tol * scale)
    return distance <= tolerance


def as_point3(value: Any, label: str) -> np.ndarray:
    """Return *value* as a finite 3D point."""

    point = np.asarray(value, dtype=float)
    if point.shape != (3,):
        raise ValueError(f"{label} must be a 3D point, got shape {point.shape}.")
    if not np.isfinite(point).all():
        raise ValueError(f"{label} contains NaN or infinite values.")
    return point.copy()


def as_polyline(value: Any, label: str) -> np.ndarray:
    """Return *value* as a finite polyline with shape ``(N, 3)``."""

    points = np.asarray(value, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"{label} must have shape (N, 3), got {points.shape}.")
    if points.shape[0] < 2:
        raise ValueError(f"{label} must contain at least two points.")
    if not np.isfinite(points).all():
        raise ValueError(f"{label} contains NaN or infinite values.")
    return points.copy()


def drop_consecutive_duplicates(
    points: np.ndarray,
    *,
    atol: float | None = None,
) -> np.ndarray:
    """Drop only numerically unresolved adjacent samples from a polyline."""

    points = np.asarray(points, dtype=float)
    if len(points) == 0:
        return points
    if atol is None:
        threshold = _FLOAT_REL_TOL * _point_scale(points)
    else:
        threshold = float(atol)
        if not np.isfinite(threshold) or threshold < 0.0:
            raise ValueError("atol must be finite and non-negative")

    keep = [0]
    for index in range(1, len(points)):
        if float(np.linalg.norm(points[index] - points[keep[-1]])) > threshold:
            keep.append(index)
    return points[np.asarray(keep, dtype=int)]


def oriented_edge_polyline(
    graph: nx.MultiGraph,
    u: Any,
    v: Any,
    key: Any,
    data: dict[str, Any],
) -> np.ndarray:
    """Return an edge polyline oriented from node *u* to node *v*."""

    start = as_point3(graph.nodes[u].get("pos"), f"node {u!r} 'pos'")
    end = as_point3(graph.nodes[v].get("pos"), f"node {v!r} 'pos'")

    if data.get("pts") is None:
        points = np.vstack([start, end])
    else:
        points = as_polyline(data["pts"], f"edge {(u, v, key)!r} 'pts'")

    forward = float(np.linalg.norm(points[0] - start) + np.linalg.norm(points[-1] - end))
    reverse = float(np.linalg.norm(points[0] - end) + np.linalg.norm(points[-1] - start))
    if reverse < forward:
        points = points[::-1].copy()

    points[0] = start
    points[-1] = end
    points = drop_consecutive_duplicates(points)
    if len(points) < 2:
        raise ValueError(f"edge {(u, v, key)!r} collapsed to fewer than two distinct points.")
    return points


def _point_segment_distance(
    point: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
) -> float:
    direction = end - start
    denom = float(np.dot(direction, direction))
    if denom <= np.finfo(float).tiny:
        return float(np.linalg.norm(point - start))
    t = float(np.clip(np.dot(point - start, direction) / denom, 0.0, 1.0))
    return float(np.linalg.norm(point - (start + t * direction)))


def _segment_segment_closest(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    d: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return distance and closest points for two closed 3-D segments."""
    u = b - a
    v = d - c
    w = a - c
    aa = float(np.dot(u, u))
    bb = float(np.dot(u, v))
    cc = float(np.dot(v, v))
    dd = float(np.dot(u, w))
    ee = float(np.dot(v, w))
    tiny = np.finfo(float).tiny

    if aa <= tiny and cc <= tiny:
        return float(np.linalg.norm(a - c)), a, c
    if aa <= tiny:
        t = float(np.clip(ee / cc, 0.0, 1.0))
        q = c + t * v
        return float(np.linalg.norm(a - q)), a, q
    if cc <= tiny:
        s = float(np.clip(-dd / aa, 0.0, 1.0))
        p = a + s * u
        return float(np.linalg.norm(p - c)), p, c

    denom = aa * cc - bb * bb
    if abs(denom) > _FLOAT_REL_TOL * max(aa * cc, tiny):
        s = float(np.clip((bb * ee - cc * dd) / denom, 0.0, 1.0))
    else:
        s = 0.0
    t = (bb * s + ee) / cc
    if t < 0.0:
        t = 0.0
        s = float(np.clip(-dd / aa, 0.0, 1.0))
    elif t > 1.0:
        t = 1.0
        s = float(np.clip((bb - dd) / aa, 0.0, 1.0))

    p = a + s * u
    q = c + t * v
    return float(np.linalg.norm(p - q)), p, q


def _geometric_embedding_issues(
    graph: nx.MultiGraph,
    positions: dict[Any, np.ndarray],
    polylines: dict[tuple[Any, Any, Any], np.ndarray],
) -> list[str]:
    """Detect 3-D contacts that are not represented by graph incidence."""
    if not polylines:
        return []

    scale = _point_scale(*positions.values(), *polylines.values())
    tolerance = max(_FLOAT_REL_TOL * scale, 1e-10 * scale)
    issues: list[str] = []
    segments: list[
        tuple[tuple[Any, Any, Any], int, np.ndarray, np.ndarray]
    ] = []

    for edge_ref, pts in polylines.items():
        for index in range(len(pts) - 1):
            segments.append((edge_ref, index, pts[index], pts[index + 1]))

    for left_index, (left_ref, left_seg, a, b) in enumerate(segments):
        lu, lv, _ = left_ref
        left_pts = polylines[left_ref]
        for right_ref, right_seg, c, d in segments[left_index + 1 :]:
            ru, rv, _ = right_ref

            if left_ref == right_ref:
                adjacent = abs(left_seg - right_seg) <= 1
                closed_adjacent = (
                    lu == lv
                    and {left_seg, right_seg} == {0, len(left_pts) - 2}
                )
                if adjacent or closed_adjacent:
                    continue

            if np.any(np.maximum(a, b) + tolerance < np.minimum(c, d)):
                continue
            if np.any(np.maximum(c, d) + tolerance < np.minimum(a, b)):
                continue

            distance, p, q = _segment_segment_closest(a, b, c, d)
            if distance > tolerance:
                continue

            common_nodes = set((lu, lv)).intersection((ru, rv))
            permitted = False
            for node in common_nodes:
                node_pos = positions.get(node)
                if node_pos is None:
                    continue
                if (
                    np.linalg.norm(p - node_pos) <= tolerance
                    and np.linalg.norm(q - node_pos) <= tolerance
                ):
                    left_other = (
                        b if np.linalg.norm(a - node_pos) <= tolerance else a
                    )
                    right_other = (
                        d if np.linalg.norm(c - node_pos) <= tolerance else c
                    )
                    if (
                        _point_segment_distance(left_other, c, d) > tolerance
                        and _point_segment_distance(right_other, a, b) > tolerance
                    ):
                        permitted = True
                        break

            if permitted:
                continue

            issues.append(
                f"edges {left_ref!r} and {right_ref!r} have an unmodeled 3D contact"
            )
            if len(issues) >= 20:
                return issues

    for node, point in positions.items():
        for edge_ref, pts in polylines.items():
            u, v, _ = edge_ref
            if node in (u, v):
                continue
            for index in range(len(pts) - 1):
                if (
                    _point_segment_distance(
                        point, pts[index], pts[index + 1]
                    )
                    <= tolerance
                ):
                    issues.append(
                        f"edge {edge_ref!r} passes through unincident vertex {node!r}"
                    )
                    break
            if len(issues) >= 20:
                return issues

    return issues


def validate_embedding(
    graph: nx.MultiGraph,
    *,
    check_geometry: bool = True,
) -> list[str]:
    """Return issues for an embedded MultiGraph(pos/pts)."""
    issues: list[str] = []
    if not isinstance(graph, nx.MultiGraph):
        return ["graph is not a networkx.MultiGraph"]
    if graph.is_directed():
        issues.append("graph must be undirected")
    if graph.number_of_nodes() == 0:
        issues.append("graph has no nodes")
        return issues

    valid_positions: dict[Any, np.ndarray] = {}
    valid_polylines: dict[tuple[Any, Any, Any], np.ndarray] = {}

    for node, data in graph.nodes(data=True):
        if "pos" not in data:
            issues.append(f"node {node!r} is missing 'pos'")
            continue
        try:
            valid_positions[node] = as_point3(
                data["pos"], f"node {node!r} pos"
            )
        except ValueError as exc:
            issues.append(str(exc))

    for u, v, key, data in graph.edges(keys=True, data=True):
        pts: np.ndarray | None = None
        if data.get("pts") is not None:
            try:
                pts = as_polyline(
                    data["pts"], f"edge {(u, v, key)!r} pts"
                )
            except ValueError as exc:
                issues.append(str(exc))
                continue

        if u not in valid_positions or v not in valid_positions:
            continue

        if pts is None:
            pts = np.vstack([valid_positions[u], valid_positions[v]])

        u_pos = valid_positions[u]
        v_pos = valid_positions[v]
        scale = _point_scale(pts, u_pos, v_pos)
        direct = (
            _points_close(pts[0], u_pos, scale=scale)
            and _points_close(pts[-1], v_pos, scale=scale)
        )
        reverse = (
            _points_close(pts[0], v_pos, scale=scale)
            and _points_close(pts[-1], u_pos, scale=scale)
        )
        if not (direct or reverse):
            issues.append(
                f"edge {(u, v, key)!r} endpoints do not match node positions"
            )
            continue

        if len(drop_consecutive_duplicates(pts)) < 2:
            issues.append(
                f"edge {(u, v, key)!r} collapsed to fewer than two distinct points"
            )
            continue

        valid_polylines[(u, v, key)] = pts

    if check_geometry and not issues:
        issues.extend(
            _geometric_embedding_issues(
                graph, valid_positions, valid_polylines
            )
        )

    return issues


def ensure_embedding(
    graph: nx.MultiGraph,
    *,
    copy: bool = True,
    normalize: bool = True,
    check_geometry: bool = False,
) -> nx.MultiGraph:
    """Validate and optionally normalize an embedded spatial graph."""
    issues = validate_embedding(
        graph, check_geometry=check_geometry
    )
    if issues:
        raise EmbeddingValidationError(issues)

    result = graph.copy() if copy else graph
    if not normalize:
        return result

    for _, data in result.nodes(data=True):
        data["pos"] = as_point3(data["pos"], "node 'pos'")

    for u, v, key, data in result.edges(keys=True, data=True):
        data["pts"] = oriented_edge_polyline(
            result, u, v, key, data
        )

    return result


def is_embedding(graph: nx.MultiGraph) -> bool:
    """Return whether graph satisfies the strict embedded-graph contract."""
    return not validate_embedding(graph, check_geometry=True)


def idx_to_coord(
    indices: ArrayLike,
    spacing: Sequence[float] = (1.0, 1.0, 1.0),
    origin: Sequence[float] = (0.0, 0.0, 0.0),
) -> NDArray:
    """Convert an array of 3D image indices to spatial coordinates."""

    array = np.asarray(indices)
    if array.shape[-1] != 3:
        raise ValueError("Input array must have shape (..., 3).")

    return array * spacing + origin


def get_all_edge_pts(G: nx.MultiGraph) -> NDArray:
    """Get all edge points from the graph as a single array."""

    graph = ensure_embedding(G, copy=False, normalize=False)
    edge_pts_list = [
        oriented_edge_polyline(graph, u, v, k, data)
        for u, v, k, data in graph.edges(keys=True, data=True)
    ]
    if not edge_pts_list:
        return np.empty((0, 3), dtype=float)
    return np.concatenate(edge_pts_list)


def total_edge_pts(G: nx.MultiGraph) -> int:
    """Count the total number of edge-polyline points."""

    return len(get_all_edge_pts(G))


def smooth_edges(
    G: nx.MultiGraph,
    epsilon: float = 0.0,
    copy: bool = True,
) -> nx.MultiGraph:
    """Simplify edge polylines without allowing topology-changing shortcuts."""
    epsilon = float(epsilon)
    if not np.isfinite(epsilon) or epsilon < 0.0:
        raise ValueError("epsilon must be finite and non-negative")

    H = ensure_embedding(G, copy=copy, normalize=True)
    if epsilon == 0.0 or H.number_of_edges() == 0:
        return H

    from knotted_graph.layout.repulsive.decimation import (
        DecimationOptions,
        decimate_curve_network,
    )

    vertices: list[np.ndarray] = []
    node_indices: dict[Any, int] = {}
    for node, data in H.nodes(data=True):
        node_indices[node] = len(vertices)
        vertices.append(np.asarray(data["pos"], dtype=float))

    edge_indices: dict[str, list[int]] = {}
    edge_refs: dict[str, tuple[Any, Any, Any]] = {}
    edge_order: list[str] = []

    for edge_number, (u, v, key, data) in enumerate(
        H.edges(keys=True, data=True)
    ):
        edge_id = f"edge_{edge_number}"
        points = oriented_edge_polyline(H, u, v, key, data)
        indices = [node_indices[u]]
        for point in points[1:-1]:
            vertices.append(np.asarray(point, dtype=float))
            indices.append(len(vertices) - 1)
        indices.append(node_indices[v])
        edge_indices[edge_id] = indices
        edge_refs[edge_id] = (u, v, key)
        edge_order.append(edge_id)

    vertex_array = np.asarray(vertices, dtype=float)
    scale = _point_scale(vertex_array)
    clearance = max(
        _FLOAT_REL_TOL * scale,
        np.finfo(float).tiny,
    )

    result = decimate_curve_network(
        vertex_array,
        edge_indices,
        tuple(edge_order),
        pinned_indices=set(node_indices.values()),
        options=DecimationOptions(
            max_passes=max(
                8,
                int(
                    np.ceil(
                        np.log2(max(2, len(vertex_array)))
                    )
                    + 4
                ),
            ),
            min_points_per_edge=2,
            clearance_fraction=0.0,
            min_clearance=clearance,
            max_deviation=epsilon,
            preserve_pinned_neighbors=True,
        ),
    )

    for edge_id in edge_order:
        u, v, key = edge_refs[edge_id]
        indices = np.asarray(
            result.edge_indices[edge_id], dtype=int
        )
        H.edges[u, v, key]["pts"] = oriented_edge_polyline(
            H,
            u,
            v,
            key,
            {"pts": result.vertices[indices]},
        )

    strict_issues = validate_embedding(
        H, check_geometry=True
    )
    if strict_issues:
        raise EmbeddingValidationError(strict_issues)

    return H


def remove_leaf_nodes(G: nx.MultiGraph) -> nx.MultiGraph:
    """Remove degree-1 leaves without deleting connected components.

    Components are labelled once before pruning. If all surviving vertices of
    one component are leaves in an iteration, one deterministic representative
    is retained. This preserves component count without per-component graph
    copies or repeated connectivity searches.
    """

    H = G.copy()
    component_of: dict[Any, int] = {}
    live_count: dict[int, int] = {}
    for component_id, nodes in enumerate(nx.connected_components(H)):
        component_nodes = list(nodes)
        live_count[component_id] = len(component_nodes)
        for node in component_nodes:
            component_of[node] = component_id

    while True:
        leaves = [node for node, degree in H.degree() if degree == 1]
        if not leaves:
            break

        leaves_by_component: dict[int, list[Any]] = {}
        for node in leaves:
            leaves_by_component.setdefault(component_of[node], []).append(node)

        to_remove: list[Any] = []
        for component_id, component_leaves in leaves_by_component.items():
            if len(component_leaves) == live_count[component_id]:
                keep = min(component_leaves, key=repr)
                component_leaves = [
                    node for node in component_leaves if node != keep
                ]
            to_remove.extend(component_leaves)
            live_count[component_id] -= len(component_leaves)

        if not to_remove:
            break
        H.remove_nodes_from(to_remove)

    return H


def contract_short_edges(
    G: nx.MultiGraph,
    min_length: float = 0.30,
    *,
    copy: bool = True,
) -> nx.MultiGraph:
    """Contract short non-loop edge occurrences without losing multigraph topology."""

    min_length = float(min_length)
    if not np.isfinite(min_length) or min_length < 0.0:
        raise ValueError("min_length must be finite and non-negative")

    H = ensure_embedding(G, copy=copy, normalize=True)

    def endpoint_distance(u: Any, v: Any) -> float:
        return float(np.linalg.norm(H.nodes[u]["pos"] - H.nodes[v]["pos"]))

    def safe_key(u: Any, v: Any, key: Any) -> Any:
        if not H.has_edge(u, v, key):
            return key
        suffix = 1
        candidate = ("contracted", key, suffix)
        while H.has_edge(u, v, candidate):
            suffix += 1
            candidate = ("contracted", key, suffix)
        return candidate

    def move_endpoint(
        pts: Any,
        old_endpoint: np.ndarray,
        new_endpoint: np.ndarray,
        other_endpoint: np.ndarray,
    ) -> np.ndarray:
        arr = np.asarray(pts, dtype=float).copy()
        if arr.ndim != 2 or arr.shape[1] != 3 or arr.shape[0] < 2:
            arr = np.vstack([old_endpoint, other_endpoint])
        if np.linalg.norm(arr[0] - old_endpoint) <= np.linalg.norm(arr[-1] - old_endpoint):
            arr[0] = new_endpoint
            arr[-1] = other_endpoint
        else:
            arr[-1] = new_endpoint
            arr[0] = other_endpoint
        return drop_consecutive_duplicates(arr)

    def move_loop(
        pts: Any,
        old_endpoint: np.ndarray,
        new_endpoint: np.ndarray,
    ) -> np.ndarray:
        arr = np.asarray(pts, dtype=float).copy()
        if arr.ndim != 2 or arr.shape[1] != 3 or arr.shape[0] < 2:
            arr = np.vstack([old_endpoint, old_endpoint])
        arr[0] = new_endpoint
        arr[-1] = new_endpoint
        arr = drop_consecutive_duplicates(arr)
        if len(arr) < 2:
            arr = np.vstack([new_endpoint, new_endpoint])
        return arr

    while True:
        candidates: list[tuple[float, str, str, Any, Any, Any]] = []
        for u, v, key in H.edges(keys=True):
            if u == v:
                continue
            length = endpoint_distance(u, v)
            if length < min_length:
                candidates.append((length, repr(u), repr(v), u, v, key))
        if not candidates:
            break

        _, _, _, u, v, contracted_key = min(candidates)
        degree_u, degree_v = H.degree[u], H.degree[v]
        if degree_u > degree_v:
            keep, kill = u, v
        elif degree_v > degree_u:
            keep, kill = v, u
        else:
            keep, kill = (u, v) if repr(u) <= repr(v) else (v, u)

        keep_pos = np.asarray(H.nodes[keep]["pos"], dtype=float)
        kill_pos = np.asarray(H.nodes[kill]["pos"], dtype=float)
        merged_pos = 0.5 * (keep_pos + kill_pos)
        incident = list(H.edges(kill, keys=True, data=True))
        keep_incident = list(H.edges(keep, keys=True, data=True))

        H.nodes[keep]["pos"] = merged_pos
        for a, b, key, data in keep_incident:
            if kill in (a, b):
                continue
            other = b if a == keep else a
            other_pos = np.asarray(H.nodes[other]["pos"], dtype=float)
            H.edges[a, b, key]["pts"] = move_endpoint(
                data.get("pts"), keep_pos, merged_pos, other_pos
            )

        contracted_removed = False
        for a, b, key, data in incident:
            other = b if a == kill else a
            if H.has_edge(a, b, key):
                H.remove_edge(a, b, key)

            if other == keep and key == contracted_key and not contracted_removed:
                contracted_removed = True
                continue

            edge_data = dict(data or {})
            if other == keep or other == kill:
                edge_data["pts"] = move_loop(
                    edge_data.get("pts"), kill_pos, merged_pos
                )
                new_key = safe_key(keep, keep, key)
                H.add_edge(keep, keep, key=new_key, **edge_data)
                continue

            other_pos = np.asarray(H.nodes[other]["pos"], dtype=float)
            edge_data["pts"] = move_endpoint(
                edge_data.get("pts"), kill_pos, merged_pos, other_pos
            )
            new_key = safe_key(keep, other, key)
            H.add_edge(keep, other, key=new_key, **edge_data)

        if kill in H:
            H.remove_node(kill)

    return ensure_embedding(H, copy=False, normalize=True)

def _append_edge_pts(path: list[np.ndarray], edge_pts: Any) -> None:
    if edge_pts is None or len(edge_pts) == 0:
        return

    pts = np.asarray(edge_pts, dtype=float)
    if np.array_equal(pts[-1], path[-1]):
        pts = pts[::-1]

    if np.array_equal(pts[0], path[-1]):
        path.extend(pts[1:])
        return

    raise RuntimeError(
        "Edge segment does not connect contiguously:\n"
        f"  current tail = {path[-1]}\n"
        f"  segment ends = ({pts[0]}, {pts[-1]})"
    )


def _edge_tag(u: Any, v: Any, key: Any) -> tuple[Any, Any, Any]:
    """Return a deterministic tag for an undirected multiedge."""
    return (u, v, key) if repr(u) <= repr(v) else (v, u, key)


def _has_cycles(G: nx.MultiGraph) -> bool:
    """Quick check for any cycle in *G*."""

    if G.number_of_edges() == 0:
        return False
    try:
        nx.find_cycle(G)
        return True
    except nx.NetworkXNoCycle:
        return False


def _collapse_component_with_junctions(
    G: nx.MultiGraph,
    comp: set[int],
    H: nx.MultiGraph,
) -> None:
    """Collapse chains inside a component that has junction nodes."""

    junctions = {node for node in comp if G.degree(node) > 2}
    for node in junctions:
        H.add_node(node, **G.nodes[node])

    seen_edges: set[tuple[int, int, int]] = set()

    for junction in junctions:
        for neighbor, edge_dict in G.adj[junction].items():
            for key, attrs in edge_dict.items():
                tag = _edge_tag(junction, neighbor, key)
                if tag in seen_edges:
                    continue
                seen_edges.add(tag)

                path_pts: list[np.ndarray] = [G.nodes[junction]["pos"]]
                _append_edge_pts(path_pts, attrs.get("pts", []))

                previous, current = junction, neighbor
                while current not in junctions and G.degree(current) == 2:
                    path_pts.append(G.nodes[current]["pos"])
                    next_candidates = [node for node in G.neighbors(current) if node != previous]
                    if not next_candidates:
                        break
                    nxt = next_candidates[0]

                    for key2, attrs2 in G[current][nxt].items():
                        tag2 = _edge_tag(current, nxt, key2)
                        if tag2 not in seen_edges:
                            seen_edges.add(tag2)
                            _append_edge_pts(path_pts, attrs2.get("pts", []))
                            break
                    previous, current = current, nxt

                path_pts.append(G.nodes[current]["pos"])
                if current not in H:
                    H.add_node(current, **G.nodes[current])

                H.add_edge(junction, current, pts=np.asarray(path_pts))


def _collapse_cycle_component(
    G: nx.MultiGraph,
    comp: set,
    H: nx.MultiGraph,
) -> None:
    """Collapse an Eulerian degree-two component to one embedded self-loop."""

    component = G.subgraph(comp).copy()
    rep = min(comp, key=repr)
    H.add_node(rep, **G.nodes[rep])

    if component.number_of_edges() == 0:
        return
    if not nx.is_eulerian(component):
        _copy_component(G, comp, H)
        return

    path_pts: list[np.ndarray] = [np.asarray(G.nodes[rep]["pos"], dtype=float)]
    source_edges: list[tuple[Any, Any, Any]] = []
    for u, v, key in nx.eulerian_circuit(component, source=rep, keys=True):
        attrs = G.edges[u, v, key]
        _append_edge_pts(path_pts, attrs.get("pts", []))
        source_edges.append((u, v, key))

    if not np.array_equal(path_pts[-1], path_pts[0]):
        path_pts.append(path_pts[0].copy())
    H.add_edge(
        rep,
        rep,
        pts=np.asarray(path_pts),
        source_edges=tuple(source_edges),
    )

def _copy_component(
    G: nx.MultiGraph,
    comp: set,
    H: nx.MultiGraph,
) -> None:
    """Copy one component, including all node and keyed-edge metadata."""

    for node in comp:
        H.add_node(node, **G.nodes[node])
    for u, v, key, data in G.subgraph(comp).edges(keys=True, data=True):
        H.add_edge(u, v, key=key, **data)


def simplify_edges(G: nx.MultiGraph) -> nx.MultiGraph:
    """Simplify degree-2 chains without discarding embedded connectivity.

    Components containing cycles or junctions are represented by embedded
    edges between their significant vertices. Acyclic components are returned
    normalized but otherwise unchanged so that paths, trees, and per-edge
    metadata cannot disappear implicitly. Use :func:`remove_leaf_nodes`
    explicitly when terminal branches should be removed.
    """

    G = ensure_embedding(G, copy=True, normalize=True)

    H = nx.MultiGraph()
    H.graph.update(G.graph)
    for comp in nx.connected_components(G):
        component = G.subgraph(comp)
        if not _has_cycles(component):
            logging.info(
                "Preserving one normalized acyclic component. Degree-2 chains "
                "are not collapsed because doing so could discard per-edge metadata."
            )
            _copy_component(G, comp, H)
        elif any(G.degree(node) > 2 for node in comp):
            _collapse_component_with_junctions(G, comp, H)
        else:
            _collapse_cycle_component(G, comp, H)

    return H
