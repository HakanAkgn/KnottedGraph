import networkx as nx
import numpy as np
import sympy as sp

from knotted_graph.applications.phase_maps import (
    boundary_filling_groups,
    make_yamada_phase_map,
    resolve_volume_mask,
    volume_topology,
)
from knotted_graph.inputs import KnotFunction


def _triangle_graph():
    graph = nx.MultiGraph()
    graph.add_node(0, pos=(0.0, 0.0, 0.0))
    graph.add_node(1, pos=(1.0, 0.0, 0.0))
    graph.add_node(2, pos=(0.0, 1.0, 0.0))
    graph.add_edge(0, 1, pts=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    graph.add_edge(1, 2, pts=[(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    graph.add_edge(2, 0, pts=[(0.0, 1.0, 0.0), (0.0, 0.0, 0.0)])
    return graph


def _patch_yamada(monkeypatch):
    import knotted_graph.invariants.yamada
    import knotted_graph.projection

    def fake_yamada(graph, variable, **kwargs):
        return variable + graph.number_of_edges()

    monkeypatch.setattr(
        knotted_graph.invariants.yamada,
        "compute_graph_yamada_polynomial",
        fake_yamada,
    )
    monkeypatch.setattr(knotted_graph.projection, "compute_yamada_polynomial", fake_yamada)


def test_unified_phase_map_accepts_nodal_bloch_factories(monkeypatch):
    from knotted_graph.applications.nodal.skeleton import NodalSkeleton

    _patch_yamada(monkeypatch)
    monkeypatch.setattr(NodalSkeleton, "skeleton_graph", lambda self, **kwargs: _triangle_graph())

    kx, ky, kz = sp.symbols("kx ky kz", real=True)

    def start(gamma):
        return (kx, ky, kz + sp.I * gamma)

    def end(gamma):
        return (kx + gamma, ky, kz + sp.I * gamma)

    result = make_yamada_phase_map(
        start,
        end,
        source_kind="nodal",
        lambdas=[0.0, 1.0],
        parameters=[0.1, 0.2],
        dimension=6,
        force_genus_zero_vertex=False,
    )

    assert result.source_kind == "nodal"
    assert result.parameter_name == "gamma"
    assert result.record_grid().shape == (2, 2)
    assert all(record.error is None for record in result.records)
    assert {record.edges for record in result.records} == {3}


def test_unified_phase_map_accepts_material_gap_and_energy_modes(monkeypatch):
    from knotted_graph.applications.material_surface import MaterialFermiSurface

    _patch_yamada(monkeypatch)
    monkeypatch.setattr(MaterialFermiSurface, "skeleton_graph", lambda self, **kwargs: _triangle_graph())

    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    h0 = sp.diag(kx, ky, kz)
    h1 = sp.diag(kx + sp.Rational(1, 10), ky, kz)

    gap_result = make_yamada_phase_map(
        h0,
        h1,
        source_kind="material",
        material_mode="gap",
        band_pair=(0, 1),
        lambdas=[0.0, 1.0],
        parameters=[0.02, 0.04],
        k_symbols=(kx, ky, kz),
        dimension=5,
        force_genus_zero_vertex=False,
    )
    energy_result = make_yamada_phase_map(
        h0,
        h1,
        source_kind="material",
        material_mode="energy",
        band_index=1,
        lambdas=[0.0, 1.0],
        parameters=[-0.1, 0.1],
        k_symbols=(kx, ky, kz),
        dimension=5,
        force_genus_zero_vertex=False,
    )

    assert gap_result.parameter_name == "gap_tol"
    assert energy_result.parameter_name == "energy"
    assert all(record.source_kind == "material" for record in gap_result.records)
    assert all(record.source_kind == "material" for record in energy_result.records)
    assert all(record.error is None for record in energy_result.records)

    alias_result = make_yamada_phase_map(
        h0,
        h1,
        source_kind="hamiltonian",
        material_mode="energy",
        band_index=1,
        lambdas=[0.0],
        parameters=[0.0],
        k_symbols=(kx, ky, kz),
        dimension=5,
        force_genus_zero_vertex=False,
    )
    assert alias_result.source_kind == "material"


def test_unified_phase_map_accepts_knot_functions(monkeypatch):
    _patch_yamada(monkeypatch)

    sample_calls = []

    def fake_sample(self, **kwargs):
        sample_calls.append(self.name)
        return object()

    def fake_graph(self, radius, *, sample=None, **kwargs):
        assert sample is not None
        return _triangle_graph()

    monkeypatch.setattr(KnotFunction, "sample", fake_sample)
    monkeypatch.setattr(KnotFunction, "to_spatial_graph", fake_graph)

    start = KnotFunction.from_function(lambda x, y, z: x + 1j * y, name="xy")
    end = KnotFunction.from_function(lambda x, y, z: x + 1j * z, name="xz")
    result = make_yamada_phase_map(
        start,
        end,
        source_kind="knot",
        lambdas=[0.0, 0.5, 1.0],
        parameters=[0.1, 0.2],
        dimension=6,
    )

    assert result.source_kind == "knot"
    assert result.parameter_name == "radius"
    assert result.record_grid().shape == (2, 3)
    assert len(sample_calls) == 3
    assert all(record.error is None for record in result.records)


def test_unified_phase_map_collapses_closed_genus_zero_masks_to_vertex(monkeypatch):
    from knotted_graph.applications.nodal.skeleton import NodalSkeleton

    _patch_yamada(monkeypatch)

    mask = __import__("numpy").zeros((7, 7, 7), dtype=bool)
    mask[2:5, 2:5, 2:5] = True
    monkeypatch.setattr(NodalSkeleton, "_interior_mask", property(lambda self: mask))

    def fail_skeleton_graph(self, **kwargs):
        raise AssertionError("closed genus-zero masks should collapse before graph extraction")

    monkeypatch.setattr(NodalSkeleton, "skeleton_graph", fail_skeleton_graph)

    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    result = make_yamada_phase_map(
        (kx, ky, kz),
        (kx + 1, ky, kz),
        source_kind="nodal",
        lambdas=[0.0],
        parameters=[0.0],
        dimension=7,
    )

    [record] = result.records
    assert record.error is None
    assert record.nodes == 1
    assert record.edges == 0
    assert record.yamada.free_symbols == {sp.Symbol("Y")}


def test_volume_topology_distinguishes_ball_shell_and_solid_torus():
    axis = np.arange(41, dtype=float) - 20.0
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    radius_squared = x * x + y * y + z * z

    ball = radius_squared <= 10.0**2
    shell = (radius_squared <= 10.0**2) & (radius_squared >= 5.0**2)
    torus = (np.sqrt(x * x + y * y) - 10.0) ** 2 + z * z <= 3.0**2

    ball_topology = volume_topology(ball)
    shell_topology = volume_topology(shell)
    torus_topology = volume_topology(torus)

    assert (
        ball_topology.connected_components,
        ball_topology.handle_rank,
        ball_topology.enclosed_voids,
        ball_topology.boundary_components,
    ) == (1, 0, 0, 1)
    assert (
        shell_topology.connected_components,
        shell_topology.handle_rank,
        shell_topology.enclosed_voids,
        shell_topology.boundary_components,
    ) == (1, 0, 1, 2)
    assert (
        torus_topology.connected_components,
        torus_topology.handle_rank,
        torus_topology.enclosed_voids,
        torus_topology.boundary_components,
    ) == (1, 1, 0, 1)

    [(outer_filling, inner_fillings)] = boundary_filling_groups(shell)
    assert volume_topology(outer_filling).boundary_components == 1
    assert len(inner_fillings) == 1
    assert volume_topology(inner_fillings[0]).boundary_components == 1


def test_resolve_volume_mask_removes_satellites_and_fills_tiny_voids():
    mask = np.zeros((25, 25, 25), dtype=bool)
    mask[4:20, 4:20, 4:20] = True
    mask[10, 10, 10] = False
    mask[22, 22, 22] = True

    resolved, diagnostics = resolve_volume_mask(
        mask,
        min_component_voxels=8,
        min_void_voxels=8,
        dominant_component_only=True,
    )

    assert resolved[10, 10, 10]
    assert not resolved[22, 22, 22]
    assert diagnostics == {
        "removed_component_voxels": 1,
        "filled_void_voxels": 1,
    }


def test_unified_phase_map_represents_a_shell_by_two_vertices(monkeypatch):
    from knotted_graph.applications.nodal.skeleton import NodalSkeleton

    _patch_yamada(monkeypatch)

    axis = np.arange(21, dtype=float) - 10.0
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    radius_squared = x * x + y * y + z * z
    shell = (radius_squared <= 7.0**2) & (radius_squared >= 3.0**2)
    monkeypatch.setattr(NodalSkeleton, "_interior_mask", property(lambda self: shell))

    def fail_skeleton_graph(self, **kwargs):
        raise AssertionError("spherical boundary components should collapse before extraction")

    monkeypatch.setattr(NodalSkeleton, "skeleton_graph", fail_skeleton_graph)

    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    result = make_yamada_phase_map(
        (kx, ky, kz),
        (kx + 1, ky, kz),
        source_kind="nodal",
        lambdas=[0.0],
        parameters=[0.0],
        dimension=21,
    )

    [record] = result.records
    assert record.error is None
    assert record.nodes == 2
    assert record.edges == 0


def test_unified_phase_map_collapses_no_core_skeletons_to_vertex(monkeypatch):
    from knotted_graph.applications.nodal.skeleton import NodalSkeleton
    from knotted_graph.core import EmbeddingValidationError

    _patch_yamada(monkeypatch)

    mask = __import__("numpy").zeros((7, 7, 7), dtype=bool)
    mask[0, 0, 0] = True
    monkeypatch.setattr(NodalSkeleton, "_interior_mask", property(lambda self: mask))

    def empty_skeleton_graph(self, **kwargs):
        raise EmbeddingValidationError(["graph has no edges"])

    monkeypatch.setattr(NodalSkeleton, "skeleton_graph", empty_skeleton_graph)

    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    result = make_yamada_phase_map(
        (kx, ky, kz),
        (kx + 1, ky, kz),
        source_kind="nodal",
        lambdas=[0.0],
        parameters=[0.0],
        dimension=7,
    )

    [record] = result.records
    assert record.error is None
    assert record.nodes == 1
    assert record.edges == 0
