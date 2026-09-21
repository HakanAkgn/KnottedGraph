from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed

import networkx as nx
import numpy as np
import pytest
import sympy as sp

from knotted_graph.applications._material_surface_base import (
    MaterialFermiSurface,
    _is_hermitian_symbolic_or_numeric,
)
from knotted_graph.applications.phase_maps import (
    _compute_yamada,
    _isolated_vertex_graph,
    volume_topology,
)
from knotted_graph.core import (
    BouquetGraph,
    contract_short_edges,
    ensure_embedding,
    simplify_edges,
    smooth_edges,
    validate_embedding,
)
from knotted_graph.core.graphs import weisfeiler_lehman_multigraph_hash
from knotted_graph.invariants.yamada import compute_graph_yamada_polynomial
from knotted_graph.projection import PDCode, compute_yamada_polynomial


A = sp.Symbol("A")


def _beta(graph: nx.MultiGraph) -> int:
    if graph.number_of_nodes() == 0:
        return 0
    return (
        graph.number_of_edges()
        - graph.number_of_nodes()
        + nx.number_connected_components(graph)
    )


def _crossing_graph(scale: float) -> nx.MultiGraph:
    s = float(scale)
    graph = nx.MultiGraph()
    points = {
        "a": np.array([-s, 0.0, 0.1 * s]),
        "b": np.array([s, 0.0, 0.1 * s]),
        "c": np.array([0.0, -s, -0.1 * s]),
        "d": np.array([0.0, s, -0.1 * s]),
    }
    for node, point in points.items():
        graph.add_node(node, pos=point)
    graph.add_edge("a", "b", pts=np.vstack([points["a"], points["b"]]))
    graph.add_edge("c", "d", pts=np.vstack([points["c"], points["d"]]))
    return graph


def _one_crossing_theta() -> nx.MultiGraph:
    graph = nx.MultiGraph()
    graph.add_node("u", pos=np.array([-2.0, 0.0, 0.0]))
    graph.add_node("v", pos=np.array([2.0, 0.0, 0.0]))
    curves = [
        np.array([[-2, 0, 0], [-1, -1, 0.5], [1, 1, 0.5], [2, 0, 0]], float),
        np.array([[-2, 0, 0], [-1, 1, -0.5], [1, -1, -0.5], [2, 0, 0]], float),
        np.array([[-2, 0, 0], [-1, 2, 0], [1, 2, 0], [2, 0, 0]], float),
    ]
    for points in curves:
        graph.add_edge("u", "v", pts=points)
    return graph


def _segment_distance(a, b, c, d):
    # Small independent closest-distance reference for the smoothing regression.
    u = b - a
    v = d - c
    w = a - c
    aa = float(np.dot(u, u))
    bb = float(np.dot(u, v))
    cc = float(np.dot(v, v))
    dd = float(np.dot(u, w))
    ee = float(np.dot(v, w))
    denom = aa * cc - bb * bb
    if denom > 1e-20:
        s = float(np.clip((bb * ee - cc * dd) / denom, 0.0, 1.0))
    else:
        s = 0.0
    t = float(np.clip((bb * s + ee) / cc, 0.0, 1.0))
    if t in (0.0, 1.0):
        s = float(np.clip((bb * t - dd) / aa, 0.0, 1.0))
    return float(np.linalg.norm((a + s * u) - (c + t * v)))


def test_empty_and_isolated_graph_semantics_remain_distinct():
    empty = nx.MultiGraph()
    isolated = nx.MultiGraph()
    isolated.add_node("v", pos=np.zeros(3))

    assert compute_graph_yamada_polynomial(empty, A) == 1
    assert compute_graph_yamada_polynomial(BouquetGraph(0), A) == -1
    assert BouquetGraph(0).number_of_nodes() == 1
    assert validate_embedding(isolated) == []

    processor = PDCode(isolated)
    assert processor.compute(rotation_angles=(0.0, 0.0, 0.0)) == "V[]"

    with pytest.raises(ValueError, match="graph has no nodes"):
        compute_yamada_polynomial(
            empty,
            A,
            rotation_angles=(0.0, 0.0, 0.0),
        )


def test_parallel_edge_contraction_preserves_cycle_rank_as_loop():
    graph = nx.MultiGraph()
    graph.add_node("u", pos=np.array([0.0, 0.0, 0.0]))
    graph.add_node("v", pos=np.array([0.1, 0.0, 0.0]))
    graph.add_edge(
        "u", "v", key="e0",
        pts=np.array([[0, 0, 0], [0.1, 0, 0]], float),
    )
    graph.add_edge(
        "u", "v", key="e1",
        pts=np.array([[0, 0, 0], [0.05, 0.05, 0], [0.1, 0, 0]], float),
    )

    reduced = contract_short_edges(graph, min_length=0.2)

    assert reduced.number_of_nodes() == 1
    assert reduced.number_of_edges() == 1
    assert _beta(reduced) == _beta(graph) == 1
    assert len(list(nx.selfloop_edges(reduced, keys=True))) == 1


def test_material_contraction_preserves_parallel_cycle_rank():
    obj = object.__new__(MaterialFermiSurface)
    graph = nx.MultiGraph()
    graph.add_node("u", pos=np.array([0.0, 0.0, 0.0]))
    graph.add_node("v", pos=np.array([0.1, 0.0, 0.0]))
    graph.add_edge(
        "u", "v", key="e0",
        pts=np.array([[0, 0, 0], [0.1, 0, 0]], float),
    )
    graph.add_edge(
        "u", "v", key="e1",
        pts=np.array([[0, 0, 0], [0.05, 0.1, 0], [0.1, 0, 0]], float),
    )
    before = _beta(graph)

    assert obj._contract_one_edge(graph, "u", "v", edge_key="e0")

    assert _beta(graph) == before == 1
    assert graph.number_of_nodes() == 1
    assert graph.number_of_edges() == 1


def test_parallel_cycle_simplification_keeps_both_geometric_branches():
    graph = nx.MultiGraph()
    graph.add_node("a", pos=np.array([0.0, 0.0, 0.0]))
    graph.add_node("b", pos=np.array([2.0, 0.0, 0.0]))
    graph.add_edge(
        "a", "b", key="upper",
        pts=np.array([[0, 0, 0], [1, 1, 0], [2, 0, 0]], float),
    )
    graph.add_edge(
        "a", "b", key="lower",
        pts=np.array([[0, 0, 0], [1, -1, 0], [2, 0, 0]], float),
    )

    simplified = simplify_edges(graph)

    assert simplified.number_of_nodes() == 1
    assert simplified.number_of_edges() == 1
    data = next(iter(simplified.edges(data=True)))[2]
    points = np.asarray(data["pts"], dtype=float)
    assert points[:, 1].max() >= 1.0
    assert points[:, 1].min() <= -1.0
    assert len(data["source_edges"]) == 2


def test_simplification_accepts_mixed_node_label_types():
    graph = nx.MultiGraph()
    positions = {
        0: np.array([0.0, 0.0, 0.0]),
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.5, 1.0, 0.0]),
    }
    for node, pos in positions.items():
        graph.add_node(node, pos=pos)
    for index, (u, v) in enumerate(((0, "a"), ("a", "b"), ("b", 0))):
        graph.add_edge(u, v, key=f"e{index}", pts=np.vstack([positions[u], positions[v]]))

    simplified = simplify_edges(graph)
    assert simplified.number_of_edges() == 1
    assert next(iter(simplified.nodes())) in positions


def test_topology_safe_smoothing_does_not_create_crossing():
    graph = nx.MultiGraph()
    positions = {
        "a": np.array([-2.0, 0.0, 0.0]),
        "b": np.array([2.0, 0.0, 0.0]),
        "c": np.array([0.0, 0.0, -1.0]),
        "d": np.array([0.0, 0.0, 1.0]),
    }
    for node, pos in positions.items():
        graph.add_node(node, pos=pos)
    graph.add_edge(
        "a", "b", key="bent",
        pts=np.array([positions["a"], [0.0, 2.0, 0.0], positions["b"]], float),
    )
    graph.add_edge(
        "c", "d", key="vertical",
        pts=np.vstack([positions["c"], positions["d"]]),
    )

    smoothed = smooth_edges(graph, epsilon=3.0)
    bent = np.asarray(smoothed.edges["a", "b", "bent"]["pts"], float)
    vertical = np.asarray(smoothed.edges["c", "d", "vertical"]["pts"], float)
    minimum = min(
        _segment_distance(bent[i], bent[i + 1], vertical[j], vertical[j + 1])
        for i in range(len(bent) - 1)
        for j in range(len(vertical) - 1)
    )
    assert minimum > 1e-6
    assert validate_embedding(smoothed, check_geometry=True) == []


@pytest.mark.parametrize("case", ["edge_contact", "overlap", "through_vertex"])
def test_strict_embedding_validation_rejects_unmodeled_3d_contacts(case):
    graph = nx.MultiGraph()
    if case == "edge_contact":
        positions = {
            "a": (-1, 0, 0), "b": (1, 0, 0),
            "c": (0, -1, 0), "d": (0, 1, 0),
        }
        edges = (("a", "b"), ("c", "d"))
    elif case == "overlap":
        positions = {
            "a": (0, 0, 0), "b": (2, 0, 0),
            "c": (1, 0, 0), "d": (3, 0, 0),
        }
        edges = (("a", "b"), ("c", "d"))
    else:
        positions = {
            "a": (-1, 0, 0), "b": (1, 0, 0),
            "v": (0, 0, 0), "w": (0, 1, 0),
        }
        edges = (("a", "b"), ("v", "w"))

    positions = {key: np.asarray(value, float) for key, value in positions.items()}
    for node, pos in positions.items():
        graph.add_node(node, pos=pos)
    for u, v in edges:
        graph.add_edge(u, v, pts=np.vstack([positions[u], positions[v]]))

    assert validate_embedding(graph, check_geometry=True)


def test_projection_crossing_count_is_uniform_scale_invariant():
    counts = []
    for scale in (1e-8, 1.0, 1e8):
        processor = PDCode(_crossing_graph(scale))
        counts.append(processor.count_crossings(rotation_angles=(0.0, 0.0, 0.0)))
    assert counts == [1, 1, 1]


def test_pdcode_concurrent_calls_are_deterministic_and_thread_safe():
    expected = PDCode(_one_crossing_theta()).compute(
        rotation_angles=(0.0, 0.0, 0.0)
    )
    outputs = []
    errors = []

    def one():
        return PDCode(_one_crossing_theta()).compute(
            rotation_angles=(0.0, 0.0, 0.0)
        )

    with ThreadPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(one) for _ in range(100)]
        for future in as_completed(futures):
            try:
                outputs.append(future.result())
            except Exception as exc:  # pragma: no cover - diagnostic assertion below
                errors.append(exc)

    assert errors == []
    assert set(outputs) == {expected}


def test_high_level_yamada_is_thread_safe():
    expected = sp.expand(
        compute_yamada_polynomial(
            _one_crossing_theta(),
            A,
            rotation_angles=(0.0, 0.0, 0.0),
            normalize=False,
        )
    )
    values = []
    errors = []

    def one():
        return sp.expand(
            compute_yamada_polynomial(
                _one_crossing_theta(),
                A,
                rotation_angles=(0.0, 0.0, 0.0),
                normalize=False,
            )
        )

    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(one) for _ in range(40)]
        for future in as_completed(futures):
            try:
                values.append(future.result())
            except Exception as exc:  # pragma: no cover
                errors.append(exc)

    assert errors == []
    assert all(sp.expand(value - expected) == 0 for value in values)


def test_phase_map_empty_and_fail_closed_semantics():
    empty_mask = np.zeros((7, 7, 7), dtype=bool)
    topology = volume_topology(empty_mask)
    graph = _isolated_vertex_graph(topology.boundary_components)
    assert topology.boundary_components == 0
    assert graph.number_of_nodes() == 0
    assert compute_graph_yamada_polynomial(graph, A) == 1

    triangle = nx.MultiGraph()
    positions = {
        "a": np.array([0.0, 0.0, 0.0]),
        "b": np.array([1.0, 0.0, 0.0]),
        "c": np.array([0.5, 1.0, 0.0]),
    }
    for node, pos in positions.items():
        triangle.add_node(node, pos=pos)
    for u, v in (("a", "b"), ("b", "c"), ("c", "a")):
        triangle.add_edge(u, v, pts=np.vstack([positions[u], positions[v]]))

    with pytest.raises(ValueError):
        _compute_yamada(triangle, A, {"rotation_order": "BAD"})


def test_material_hermiticity_check_does_not_certify_finite_sample_zero():
    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    a = sp.Float("0.5773502691896258")
    f = kx * (kx + a) * (kx - a) * (kx + sp.Float("0.26")) * (
        kx - sp.Float("0.66")
    )
    matrix = sp.Matrix([[sp.I * f, 0], [0, 0]])
    assert not _is_hermitian_symbolic_or_numeric(
        matrix,
        (kx, ky, kz),
        ((-1, 1), (-1, 1), (-1, 1)),
    )


def test_large_absolute_coordinates_do_not_hide_endpoint_mismatch():
    graph = nx.MultiGraph()
    graph.add_node("u", pos=np.array([1e9, 0.0, 0.0]))
    graph.add_node("v", pos=np.array([1e9 + 1e6, 0.0, 0.0]))
    graph.add_edge(
        "u", "v",
        pts=np.array([[1e9 + 1000.0, 0.0, 0.0], [1e9 + 1e6, 0.0, 0.0]]),
    )
    assert any(
        "endpoints do not match" in issue
        for issue in validate_embedding(graph, check_geometry=False)
    )


def test_tiny_embedding_remains_normalizable():
    graph = nx.MultiGraph()
    graph.add_node("u", pos=np.array([0.0, 0.0, 0.0]))
    graph.add_node("v", pos=np.array([1e-12, 0.0, 0.0]))
    graph.add_edge(
        "u", "v",
        pts=np.array([[0.0, 0.0, 0.0], [1e-12, 0.0, 0.0]]),
    )
    normalized = ensure_embedding(graph)
    assert normalized.number_of_edges() == 1


def test_simple_graph_wl_hash_is_supported():
    assert isinstance(weisfeiler_lehman_multigraph_hash(nx.path_graph(3)), str)
