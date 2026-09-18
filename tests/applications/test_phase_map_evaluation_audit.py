import networkx as nx
import numpy as np
import pytest
import sympy as sp

from knotted_graph.applications import phase_maps
from knotted_graph.core import remove_leaf_nodes
from knotted_graph.extraction import skeleton_image_to_graph


def triangle():
    g = nx.MultiGraph()
    for node, pos in enumerate([(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)]):
        g.add_node(node, pos=pos)
    for u, v in [(0, 1), (1, 2), (2, 0)]:
        g.add_edge(u, v, pts=[g.nodes[u]["pos"], g.nodes[v]["pos"]])
    return g


def test_failed_spatial_evaluation_never_falls_back_to_abstract(monkeypatch):
    import knotted_graph.projection
    import knotted_graph.invariants.yamada
    def failed(*args, **kwargs):
        raise RuntimeError("nongeneric projection")
    def abstract(*args, **kwargs):
        pytest.fail("an abstract polynomial is not a spatial fallback")
    monkeypatch.setattr(knotted_graph.projection, "compute_yamada_polynomial", failed)
    monkeypatch.setattr(knotted_graph.invariants.yamada, "compute_graph_yamada_polynomial", abstract)
    with pytest.raises(RuntimeError, match="nongeneric"):
        phase_maps._compute_yamada(triangle(), sp.Symbol("Y"), {"normalize": True})


def test_failed_cell_keeps_graph_and_error(monkeypatch):
    g = triangle()
    monkeypatch.setattr(phase_maps, "_nodal_factories", lambda *a, **kw: (None, lambda *a: g))
    def failed(*args, **kwargs):
        raise RuntimeError("nongeneric projection")
    monkeypatch.setattr(phase_maps, "_compute_yamada_audited", failed)
    y = sp.Symbol("Y")
    result = phase_maps.make_yamada_phase_map((y,y,y), source_kind="nodal", lambdas=[0], parameters=[0])
    record = result.records[0]
    assert record.yamada is None and record.nodes == 3
    assert record.evaluation_kind == "failed"
    assert "nongeneric" in record.error
    assert record.phase_signature.startswith("error:")


def test_success_records_projection_and_normalization():
    y = sp.Symbol("Y")
    value, audit = phase_maps._compute_yamada_audited(triangle(), y, {"normalize": True})
    assert min(term.as_powers_dict().get(y, 0) for term in sp.expand(value).as_ordered_terms()) == 0
    assert audit["normalization"] == "signed-minimum-degree-zero"
    assert audit["projection"]["pd_code"]
    assert audit["projection"]["num_crossings"] == 0
    assert audit["is_subcubic"]


def test_higher_valence_is_a_diagram_signature():
    g = nx.MultiGraph()
    g.add_edges_from([(0,1)]*4)
    assert phase_maps._phase_signature(g, sp.Symbol("Y"), None).startswith("diagram-yamada:")


def test_isolated_vertices_have_exact_multiplicative_value():
    y = sp.Symbol("Y")
    for n in range(1, 5):
        g = nx.MultiGraph()
        g.add_nodes_from(range(n))
        value, audit = phase_maps._compute_yamada_audited(g, y, {"normalize": True})
        assert value == (-1)**n
        assert audit["evaluation_kind"] == "isolated-vertices"
        assert audit["projection"] is None


def test_skeleton_extraction_preserves_two_voxel_and_isolated_components():
    mask = np.zeros((20,20,20), dtype=bool)
    mask[2,2,2] = True
    mask[8:10,8,8] = True
    mask[14:18,15,15] = True
    graph = skeleton_image_to_graph(mask)
    assert nx.number_connected_components(graph) == 3
    assert sorted(len(c) for c in nx.connected_components(graph)) == [1,2,2]
    cleaned = remove_leaf_nodes(graph)
    assert len(cleaned) == 3
    assert cleaned.number_of_edges() == 0


def test_leaf_cleanup_retains_each_tree_component_and_cycle():
    g = triangle()
    g.add_edges_from([(10,11),(20,21),(21,22)])
    g.add_node(30)
    cleaned = remove_leaf_nodes(g)
    assert nx.number_connected_components(cleaned) == 4
    assert cleaned.number_of_edges() == 3


def test_isolated_component_label_uses_computed_polynomial():
    from knotted_graph.applications.phase_map_examples._tpms import short_signature_label

    assert "1" in short_signature_label("vertex:1", "vertex", "1")
    assert "-1" not in short_signature_label("vertex:1", "vertex", "1")
    assert "one-vertex" not in short_signature_label("vertex:1", "vertex", "1")


def test_raw_phase_map_does_not_claim_normalized_output(monkeypatch):
    graph = triangle()
    monkeypatch.setattr(phase_maps, "_nodal_factories", lambda *a, **kw: (None, lambda *a: graph))
    y = sp.Symbol("Y")
    result = phase_maps.make_yamada_phase_map(
        (y, y, y), source_kind="nodal", lambdas=[0], parameters=[0],
        normalize_yamada=False,
    )
    assert result.metadata["normalize_yamada"] is False
    assert result.records[0].normalization == "raw"
    assert "Only normalized subcubic" in result.metadata["classification_scope"]
