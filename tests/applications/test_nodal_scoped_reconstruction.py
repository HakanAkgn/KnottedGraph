from copy import deepcopy
import inspect

import networkx as nx
import numpy as np
import pytest
import sympy as sp

from knotted_graph.applications.nodal.deformation import (
    NodalBlochPath, NodalPhaseRecord, NodalPhaseScan, NodalPhaseScanResult,
)
from knotted_graph.applications.nodal.skeleton import NodalSkeleton
from knotted_graph.extraction.cubical_spine import sample_grid_breaks, verify_cubical_retract


def model_for_mask(mask):
    model = NodalSkeleton.__new__(NodalSkeleton)
    model.__dict__["spectrum"] = np.where(mask, 1j, 1.).astype(complex)
    model.skeleton_graph_cache = None
    model.skeleton_graph_cache_args = None
    model._pv_data_args = None
    return model


def record(lam, value, *, error=None, kind="spatial-yamada"):
    return NodalPhaseRecord(lam, .3, None if value is None else sp.Integer(value),
                            f"v:{value}" if error is None else "error:" + error,
                            error, kind, "signed-minimum-degree-zero", True)


def test_errors_are_masked_and_never_become_phase_transitions():
    result = NodalPhaseScanResult(np.array([0., .5, 1.]), np.array([.3]),
                                 [record(0., 1), record(.5, None, error="projection"), record(1., 2)])
    labels, names = result.phase_grid()
    assert labels[0, 1] == -1 and len(names) == 2
    assert result.transition_intervals() == []
    assert result.availability_grid().tolist() == [[True, False, True]]


def test_missing_polynomial_without_error_is_still_unavailable():
    result = NodalPhaseScanResult(np.array([0.]), np.array([.3]), [record(0., None)])
    labels, names = result.phase_grid()
    assert labels.tolist() == [[-1]] and names == {}


def test_successful_signature_change_is_not_labeled_source_transition():
    result = NodalPhaseScanResult(np.array([0., 1.]), np.array([.3]), [record(0., 1), record(1., 2)])
    change, = result.transition_intervals()
    assert change["source_topology_transition_certified"] is False
    assert change["evidence"] == "evaluated_spine_signature_change"


def test_different_evaluation_scopes_are_not_compared():
    result = NodalPhaseScanResult(np.array([0., 1.]), np.array([.3]),
                                 [record(0., 1), record(1., 2, kind="diagram-yamada")])
    assert result.transition_intervals() == []


@pytest.mark.parametrize("records", [[record(0., 1), record(0., 2)], [], [record(.7, 1)]])
def test_duplicate_missing_and_off_grid_records_are_rejected(records):
    with pytest.raises(ValueError):
        NodalPhaseScanResult(np.array([0.]), np.array([.3]), records).record_grid()


@pytest.mark.parametrize("axis", [[float("nan")], [float("inf")], [0., 0.], [1., 0.], []])
def test_scan_rejects_invalid_parameter_arrays(axis):
    path = NodalBlochPath(lambda g: (0, 0, g), lambda g: (0, 0, g))
    with pytest.raises(ValueError):
        NodalPhaseScan(path, lambdas=axis, gammas=[.3])


@pytest.mark.parametrize("dimension", [True, 1, 2.5])
def test_scan_rejects_invalid_dimension(dimension):
    path = NodalBlochPath(lambda g: (0, 0, g), lambda g: (0, 0, g))
    with pytest.raises(ValueError):
        NodalPhaseScan(path, lambdas=[0.], gammas=[.3], dimension=dimension)


def test_scientific_default_does_not_smooth_geometry():
    assert inspect.signature(NodalSkeleton.skeleton_graph).parameters["smooth_epsilon"].default == 0
    assert inspect.signature(NodalPhaseScan).parameters["reconstruction"].default == "cubical"


def test_public_cubical_path_keeps_components_and_replay():
    mask = np.zeros((7, 7, 7), dtype=bool)
    mask[1:3, 1:3, 1:3] = True
    mask[5, 5, 5:7] = True
    model = model_for_mask(mask)
    graph = model.skeleton_graph(reconstruction="cubical")
    assert nx.number_connected_components(graph) == 2
    assert graph.graph["reconstruction"]["voxel_retraction_certified"]
    axes = tuple(sample_grid_breaks(np.arange(n, dtype=float)) for n in mask.shape)
    assert verify_cubical_retract(mask, model.spine_certificate, gridlines=axes)["valid"]
    assert graph.graph["reconstruction"]["analytic_source_correspondence_certified"] is False


def test_boundary_touching_source_is_not_silently_periodic():
    mask = np.ones((3, 3, 3), dtype=bool)
    graph = model_for_mask(mask).skeleton_graph(reconstruction="cubical")
    evidence = graph.graph["reconstruction"]
    assert evidence["source_touches_boundary"] is True
    assert evidence["periodic_identification"] is False
    assert graph.number_of_nodes() == 1


@pytest.mark.parametrize("kwargs", [{"smooth_epsilon": 1}, {"skeleton_image": np.ones((3, 3, 3))}])
def test_cubical_witness_cannot_be_attached_to_modified_geometry(kwargs):
    with pytest.raises(ValueError):
        model_for_mask(np.ones((3, 3, 3), dtype=bool)).skeleton_graph(reconstruction="cubical", **kwargs)


def test_graph_mutation_invalidates_cached_certificate():
    model = model_for_mask(np.ones((3, 3, 3), dtype=bool))
    first = model.skeleton_graph(reconstruction="cubical")
    original = deepcopy(first.nodes[next(iter(first))]["pos"])
    first.nodes[next(iter(first))]["pos"] = [999., 999., 999.]
    second = model.skeleton_graph(reconstruction="cubical")
    assert second is not first
    assert tuple(second.nodes[next(iter(second))]["pos"]) == tuple(original)


def test_clear_cache_drops_projection_and_witness():
    model = model_for_mask(np.ones((3, 3, 3), dtype=bool))
    model.skeleton_graph(reconstruction="cubical")
    model.__dict__["PDCode"] = "stale"
    model.clear_cache()
    assert model.spine_certificate is None and "PDCode" not in model.__dict__


def test_scan_evaluation_failure_remains_unavailable(monkeypatch):
    import knotted_graph.applications.phase_maps as maps
    def fail(*args, **kwargs):
        raise RuntimeError("deliberate projection failure")
    monkeypatch.setattr(NodalSkeleton, "skeleton_graph", lambda *a, **k: nx.MultiGraph([(0, 1)]))
    monkeypatch.setattr(maps, "_compute_yamada_audited", fail)
    x, y, z = sp.symbols("x y z", real=True)
    path = NodalBlochPath(lambda g: (x, sp.I*g, y+z), lambda g: (x, sp.I*g, y+z))
    result = NodalPhaseScan(path, lambdas=[0.], gammas=[.3], dimension=3).run()
    assert result.phase_grid()[0].tolist() == [[-1]]
    assert result.records[0].yamada is None
    assert "deliberate projection failure" in result.records[0].error
