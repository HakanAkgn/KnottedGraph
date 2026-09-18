import networkx as nx
import numpy as np

from knotted_graph.applications.nodal import skeleton as module
from knotted_graph.applications.nodal.skeleton import NodalSkeleton


def model_for(mask):
    model = NodalSkeleton.__new__(NodalSkeleton)
    model.__dict__["spectrum"] = np.where(mask, 1j, 1.).astype(complex)
    model.skeleton_graph_cache = None
    model.skeleton_graph_cache_args = None
    model._pv_data_args = None
    return model


def test_equal_bytes_with_different_shapes_do_not_share_cache(monkeypatch):
    calls = []
    def extract(image):
        calls.append(image.shape)
        graph = nx.MultiGraph()
        graph.add_node(0, pos=(0., 0., 0.))
        return graph
    monkeypatch.setattr(module, "skeleton_image_to_graph", extract)
    model = model_for(np.ones((2, 2, 3), dtype=bool))
    model.skeleton_graph(simplify=False, skeleton_image=np.ones((2, 2, 3), dtype=bool))
    model.skeleton_graph(simplify=False, skeleton_image=np.ones((2, 3, 2), dtype=bool))
    assert calls == [(2, 2, 3), (2, 3, 2)]


def test_modified_mask_cannot_reuse_old_thinning(monkeypatch):
    import knotted_graph.extraction as extraction
    calls = []
    def thin(mask):
        calls.append(mask.copy())
        return mask.copy()
    def extract(image, **kwargs):
        assert kwargs == {"expected_components": 1, "expected_cycle_rank": 0}
        graph = nx.MultiGraph()
        graph.add_node(0, pos=tuple(np.argwhere(image).mean(axis=0)))
        return graph
    monkeypatch.setattr(extraction, "skeletonize_volume", thin)
    monkeypatch.setattr(module, "skeleton_image_to_graph", extract)
    first = np.zeros((5, 5, 5), dtype=bool)
    first[1:3, 1, 1] = True
    last = np.zeros_like(first)
    last[2:4, 3, 3] = True
    model = model_for(first)
    model.__dict__["_skeleton_image"] = np.zeros_like(first)
    graph1 = model.skeleton_graph(simplify=False)
    model.__dict__["spectrum"] = np.where(last, 1j, 1.).astype(complex)
    graph2 = model.skeleton_graph(simplify=False)
    assert len(calls) == 2 and np.array_equal(calls[-1], last)
    assert graph1.nodes[0]["pos"] != graph2.nodes[0]["pos"]


def test_mutated_witness_is_recomputed():
    model = model_for(np.ones((3, 3, 3), dtype=bool))
    first = model.skeleton_graph(reconstruction="cubical")
    expected_hash = first.graph["reconstruction"]["certificate_sha256"]
    model.spine_certificate["collapses"].clear()
    last = model.skeleton_graph(reconstruction="cubical")
    assert last is not first
    assert last.graph["reconstruction"]["certificate_sha256"] == expected_hash
    assert model.spine_certificate["collapses"]


def test_false_source_claim_in_cached_metadata_is_not_reused():
    model = model_for(np.ones((3, 3, 3), dtype=bool))
    first = model.skeleton_graph(reconstruction="cubical")
    first.graph["reconstruction"]["analytic_source_correspondence_certified"] = True
    last = model.skeleton_graph(reconstruction="cubical")
    assert last is not first
    assert last.graph["reconstruction"]["analytic_source_correspondence_certified"] is False


def test_damaged_cached_geometry_is_recomputed():
    model = model_for(np.ones((3, 3, 3), dtype=bool))
    first = model.skeleton_graph(reconstruction="cubical")
    del first.nodes[next(iter(first))]["pos"]
    last = model.skeleton_graph(reconstruction="cubical")
    assert last is not first and "pos" in last.nodes[next(iter(last))]
