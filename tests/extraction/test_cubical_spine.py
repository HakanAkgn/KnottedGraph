"""Retraction witnesses, geometry, component preservation and adversarial replay."""
from copy import deepcopy
import json

import networkx as nx
import numpy as np
import pytest

from knotted_graph.extraction.cubical_spine import (
    cubical_retract, verify_cubical_retract, certificate_json, sample_grid_breaks,
)


def checked(mask, **kwargs):
    result = cubical_retract(mask, **kwargs)
    certificate = json.loads(certificate_json(result.certificate))
    options = {k: v for k, v in kwargs.items() if k == "gridlines"}
    assert verify_cubical_retract(mask, certificate, **options)["valid"]
    initial = certificate["initial_cell_counts"]
    final = certificate["terminal_cell_counts"]
    assert sum((-1)**d*n for d, n in enumerate(initial)) == sum((-1)**d*n for d, n in enumerate(final))
    return result


@pytest.mark.parametrize("shape", [(1, 1, 1), (2, 1, 1), (2, 2, 2), (3, 4, 2)])
def test_solid_blocks_have_checked_graph_retract(shape):
    result = checked(np.ones(shape, dtype=bool))
    assert result.certificate["kind"] == "graph_retract"
    assert result.certificate["graph_summary"]["components"] == 1
    assert result.certificate["graph_summary"]["cycle_rank"] == 0


def test_empty_is_not_an_invented_vertex():
    result = checked(np.zeros((2, 2, 2), dtype=bool))
    assert result.certificate["kind"] == "empty"
    assert result.graph is None


def test_tiny_separate_components_survive():
    mask = np.zeros((5, 3, 3), dtype=bool)
    mask[0, 0, 0] = True
    mask[3:5, 2, 2] = True
    result = checked(mask)
    assert nx.number_connected_components(result.graph) == 2
    assert result.certificate["graph_summary"]["cycle_rank"] == 0


def test_closed_cubes_record_diagonal_contact_without_manifold_claim():
    mask = np.zeros((2, 2, 2), dtype=bool)
    mask[0, 0, 0] = mask[1, 1, 1] = True
    result = checked(mask)
    assert nx.number_connected_components(result.graph) == 1
    assert not result.certificate["manifold_or_regular_neighborhood_certified"]


@pytest.mark.parametrize("holes", [1, 2])
def test_tunnel_count_and_full_polylines(holes):
    mask = np.ones((2 * holes + 3, 5, 2), dtype=bool)
    for i in range(holes):
        mask[2 + 2 * i, 2, :] = False
    result = checked(mask)
    assert result.graph is not None
    assert result.certificate["graph_summary"]["cycle_rank"] == holes
    for _, _, data in result.graph.edges(data=True):
        steps = np.diff(data["pts"], axis=0)
        assert np.all(np.count_nonzero(steps, axis=1) == 1)
        assert np.all(np.sum(np.abs(steps), axis=1) == 1)


def test_truncated_work_is_not_a_graph_certificate():
    mask = np.ones((3, 3, 3), dtype=bool)
    result = checked(mask, max_collapses=1)
    assert result.graph is None
    assert result.certificate["kind"] == "partial_retract"
    assert result.certificate["stop_reason"] == "collapse_budget"


def test_source_identity_includes_coordinates_and_mask():
    mask = np.ones((2, 2, 2), dtype=bool)
    result = checked(mask)
    altered = mask.copy()
    altered[0, 0, 0] = False
    assert not verify_cubical_retract(altered, result.certificate)["valid"]
    axes = (np.arange(3.) * 2, np.arange(3.), np.arange(3.))
    assert not verify_cubical_retract(mask, result.certificate, gridlines=axes)["valid"]


@pytest.mark.parametrize("change", ["truncate", "duplicate", "terminal", "summary", "claim", "bad_pair"])
def test_tampering_is_rejected(change):
    mask = np.ones((2, 2, 2), dtype=bool)
    cert = deepcopy(checked(mask).certificate)
    if change == "truncate":
        cert["collapses"].pop()
    elif change == "duplicate":
        cert["collapses"].append(cert["collapses"][0])
    elif change == "terminal":
        cert["terminal_cells"] = []
    elif change == "summary":
        cert["graph_summary"]["components"] = 7
    elif change == "claim":
        cert["analytic_source_correspondence_certified"] = True
    else:
        cert["collapses"][0][0][0] = 0.5
    assert not verify_cubical_retract(mask, cert)["valid"]


def test_irregular_physical_grid_is_preserved():
    mask = np.ones((3, 2, 2), dtype=bool)
    axes = (np.array([-2., -1., 0., 4.]), np.array([0., .1, 1.]), np.array([2., 3., 8.]))
    result = checked(mask, gridlines=axes)
    for _, data in result.graph.nodes(data=True):
        for i, x in enumerate(data["pos"]):
            assert x in axes[i]


def test_finite_domain_dual_cells_do_not_extend_or_glue_faces():
    breaks = sample_grid_breaks([-3., -1., 1., 3.])
    np.testing.assert_array_equal(breaks, [-3., -2., 0., 2., 3.])
    mask = np.zeros((4, 2, 2), dtype=bool)
    mask[0, :, :] = mask[-1, :, :] = True
    result = checked(mask, gridlines=(breaks, np.arange(3.), np.arange(3.)))
    assert result.certificate["source"]["periodic_identification"] is False
    assert result.certificate["graph_summary"]["components"] == 2


@pytest.mark.parametrize("bad", [np.ones((2, 2)), np.ones((2, 2, 2)), np.array([[[np.nan]]])])
def test_mask_type_is_not_silently_changed(bad):
    with pytest.raises(ValueError):
        cubical_retract(bad)


@pytest.mark.parametrize("budget", [0, -1, True, 1.5])
def test_invalid_budget(budget):
    with pytest.raises(ValueError):
        cubical_retract(np.ones((1, 1, 1), dtype=bool), max_cells=budget)


def test_cell_budget_is_explicit():
    with pytest.raises(ValueError, match="budget"):
        cubical_retract(np.ones((2, 2, 2), dtype=bool), max_cells=10)


@pytest.mark.parametrize("samples", [[1.], [0., 0.], [1., 0.], [0., np.inf]])
def test_invalid_sampling_grid(samples):
    with pytest.raises(ValueError):
        sample_grid_breaks(samples)
