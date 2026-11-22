import numpy as np
import pytest

from knotted_graph.periodic import analyze_periodic_mask, analyze_periodic_sublevel
from knotted_graph.periodic import _partition


def test_empty_and_full_torus():
    empty = analyze_periodic_mask(np.zeros((3, 4, 5), bool))
    full = analyze_periodic_mask(np.ones((3, 4, 5), bool))
    assert empty.occupied_betti == full.complementary_betti == (0, 0, 0, 0)
    assert full.occupied_betti == empty.complementary_betti == (1, 3, 3, 1)
    assert not empty.boundaries and not full.boundaries


def test_torus_minus_two_balls():
    field = np.full((10, 10, 10), -1.0)
    field[0, 0, 0] = field[5, 5, 5] = 1
    result = analyze_periodic_sublevel(field, 0)
    assert result.occupied_betti == (1, 3, 4, 0)
    assert result.complementary_betti == (2, 0, 0, 0)
    assert [(b.genus, b.winding_basis) for b in result.boundaries] == [(0, ()), (0, ())]
    assert result.has_same_invariants(analyze_periodic_sublevel(np.roll(field, (2, 4, 3), (0, 1, 2)), 0))


def test_two_nonseparating_boundary_surfaces_make_incidence_cycle():
    volume = np.zeros((7, 8, 9), bool)
    volume[2:5, :, :] = True
    result = analyze_periodic_mask(volume)
    assert result.occupied_betti == result.complementary_betti == (1, 2, 1, 0)
    assert len(result.regions) == len(result.boundaries) == 2
    for boundary in result.boundaries:
        assert boundary.genus == 1
        assert boundary.winding_basis == ((0, 1, 0), (0, 0, 1))
        assert boundary.occupied_region == result.boundaries[0].occupied_region
        assert boundary.complementary_region == result.boundaries[0].complementary_region


def test_equal_betti_with_different_marked_winding():
    straight = np.zeros((9, 9, 9), bool)
    straight[:, 4, 4] = True
    diagonal = np.zeros_like(straight)
    values = np.arange(9)
    diagonal[values, values, 4] = True
    straight_result = analyze_periodic_mask(straight)
    diagonal_result = analyze_periodic_mask(diagonal)
    assert straight_result.occupied_betti == diagonal_result.occupied_betti == (1, 1, 0, 0)
    assert straight_result.complementary_betti == diagonal_result.complementary_betti == (1, 3, 2, 0)
    assert straight_result.boundaries[0].genus == diagonal_result.boundaries[0].genus == 1
    assert straight_result.boundaries[0].winding_basis == ((1, 0, 0),)
    assert diagonal_result.boundaries[0].winding_basis == ((1, 1, 0),)
    assert not straight_result.has_same_invariants(diagonal_result)
    moved = np.roll(diagonal, (2, 4, 7), (0, 1, 2))
    assert diagonal_result.has_same_invariants(analyze_periodic_mask(moved))


def test_python_grid_and_surface_reference(monkeypatch):
    mask = np.random.default_rng(20261121).random((6, 7, 8)) < .4
    result = analyze_periodic_mask(mask)
    monkeypatch.setattr(_partition, "_periodic_native", None)
    assert result.has_same_invariants(analyze_periodic_mask(mask))


def test_relative_windings_are_integral_and_unsaturated():
    edges = [[0, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    shifts = [[1, 0, 0], [0, 0, 0], [1, 0, 0], [0, 0, 0], [1, 1, 0]]
    labels, wind = _partition.graph_partition(3, edges, shifts)
    assert labels.tolist() == [0, 0, 0]
    assert _partition.winding_by_component(wind, 1) == (((2, 0, 0), (1, 1, 0)),)


def test_threshold_must_avoid_vertex_energy():
    with pytest.raises(ValueError, match="equals a vertex"):
        analyze_periodic_sublevel(np.zeros((4, 4, 4)), 0)
