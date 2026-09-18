"""Adversarial tests for the experimental exact contraction primitive."""
import math
from fractions import Fraction

import networkx as nx

from knotted_graph.core.embedded_contraction import (
    certify_straight_edge_contraction, compare_by_certified_contractions,
    embedded_fingerprint, _point, _segment_intersection, _triangle_segment,
)


def graph_from_paths(paths):
    graph = nx.MultiGraph()
    for u,v,points in paths:
        graph.add_node(u,pos=points[0],label=str(u))
        graph.add_node(v,pos=points[-1],label=str(v))
        graph.add_edge(u,v,pts=points,color="kept")
    return graph


def clean_graph():
    return graph_from_paths([
        ("u","v",[(0,0,0),(1,0,0)]),
        ("u","a",[(0,0,0),(0,1,0)]),
        ("v","b",[(1,0,0),(2,0,1)]),
    ])


def test_exact_intersection_predicates_cover_skew_coplanar_and_contacts():
    p = _point
    assert _segment_intersection(p((0,0,0)),p((1,0,0)),p((.5,-1,0)),p((.5,1,0))) == (p((.5,0,0)),)
    assert not _segment_intersection(p((0,0,0)),p((1,0,0)),p((.5,-1,1)),p((.5,1,1)))
    assert len(_segment_intersection(p((0,0,0)),p((2,0,0)),p((1,0,0)),p((3,0,0)))) == 2
    triangle = tuple(map(p,[(0,0,0),(1,0,0),(0,1,0)]))
    assert _triangle_segment(triangle,p((.25,.25,-1)),p((.25,.25,1))) == (p((.25,.25,0)),)
    assert len(_triangle_segment(triangle,p((-.5,.25,0)),p((1,.25,0)))) == 2
    assert not _triangle_segment(triangle,p((.25,.25,1e-14)),p((.5,.25,1e-14)))


def test_clean_contraction_preserves_data_and_does_not_mutate_input():
    graph = clean_graph()
    before = embedded_fingerprint(graph)
    result = certify_straight_edge_contraction(graph,("u","v",0))
    assert result.status == "certified", result.reason
    assert embedded_fingerprint(graph) == before
    assert len(result.graph) == len(graph)-1
    assert result.graph.number_of_edges() == graph.number_of_edges()-1
    assert result.graph.nodes["v"]["contraction_node_provenance"][0]["node"] == "u"
    assert result.graph["v"]["a"][0]["color"] == "kept"
    assert result.witness["swept_triangles"]
    replay = certify_straight_edge_contraction(graph,result.witness["edge"])
    assert replay.witness["output_fingerprint"] == result.witness["output_fingerprint"]


def test_piercing_strand_blocks_move_even_when_endpoints_are_embedded():
    graph = clean_graph()
    graph.add_node("x",pos=(.25,.25,-1))
    graph.add_node("y",pos=(.25,.25,1))
    graph.add_edge("x","y",pts=[(.25,.25,-1),(.25,.25,1)])
    result = certify_straight_edge_contraction(graph,("u","v",0))
    assert result.status == "unknown"
    assert result.reason == "swept_triangle_meets_fixed_graph"


def test_near_miss_is_distinct_from_contact():
    graph = clean_graph()
    for z, expected in [(1e-14,"certified"),(0,"unknown")]:
        copy = graph.copy()
        copy.add_node("x",pos=(.2,.2,z))
        copy.add_node("y",pos=(.3,.2,z))
        copy.add_edge("x","y",pts=[(.2,.2,z),(.3,.2,z)])
        assert certify_straight_edge_contraction(copy,("u","v",0)).status == expected


def test_overlapping_fans_and_unintended_endpoint_contact_rejected():
    graph = clean_graph()
    graph.add_node("c",pos=(.2,.2,0))
    graph.add_edge("u","c",pts=[(0,0,0),(.2,.2,0)])
    assert certify_straight_edge_contraction(graph,("u","v",0)).status == "unknown"


def test_invalid_initial_crossing_is_not_certified():
    graph = clean_graph()
    graph.add_node("x",pos=(.5,-1,0))
    graph.add_node("y",pos=(.5,1,0))
    graph.add_edge("x","y",pts=[(.5,-1,0),(.5,1,0)])
    assert certify_straight_edge_contraction(graph,("u","v",0)).status == "invalid_input"


def test_budget_bent_edge_and_loop_return_unknown():
    graph = clean_graph()
    assert certify_straight_edge_contraction(graph,("u","v",0),max_predicates=0).status == "unknown"
    graph["u"]["v"][0]["pts"] = [(0,0,0),(.5,0,.2),(1,0,0)]
    assert certify_straight_edge_contraction(graph,("u","v",0)).reason == "selected_edge_not_single_straight_segment"
    loop = graph_from_paths([("x","x",[(0,0,0),(1,0,0),(0,1,0),(0,0,0)])])
    assert certify_straight_edge_contraction(loop,("x","x",0)).reason == "loop_contraction_forbidden"


def test_parallel_edges_are_preserved_when_one_becomes_a_loop():
    graph = graph_from_paths([
        ("u","v",[(0,0,0),(1,0,0)]),
        ("u","v",[(0,0,0),(0,1,0),(1,1,1),(1,0,0)]),
        ("u","v",[(0,0,0),(0,-1,0),(1,-1,-1),(1,0,0)]),
    ])
    result = certify_straight_edge_contraction(graph,("u","v",0))
    assert result.status == "certified", result.reason
    assert nx.number_of_selfloops(result.graph) == 2
    assert result.witness["max_degree_after"] == 4


def test_existing_loop_moves_both_halfedges_without_becoming_a_chord():
    graph = graph_from_paths([
        ("u","v",[(0,0,0),(1,0,0)]),
        ("u","u",[(0,0,0),(0,1,0),(0,1,1),(0,0,1),(0,0,0)]),
    ])
    result = certify_straight_edge_contraction(graph,("u","v",0))
    assert result.status == "certified", result.reason
    assert nx.number_of_selfloops(result.graph) == 1
    points = next(iter(result.graph.edges(data=True)))[2]['pts']
    assert len(points) == 5 and points[0] == points[-1] == (1.,0.,0.)


def test_rational_coordinates_are_not_silently_rounded_in_certified_output():
    graph = graph_from_paths([
        ("u","v",[(0,0,0),(1,0,0)]),
        ("u","a",[(0,0,0),(0,Fraction(1,3),0)]),
    ])
    result = certify_straight_edge_contraction(graph,("u","v",0))
    assert result.status == "unknown" and result.reason == "output_coordinate_rounding"


def test_true_witness_comparison_and_unknown_budget():
    left = clean_graph()
    right = certify_straight_edge_contraction(left,("u","v",0)).graph
    result = compare_by_certified_contractions(left,right,max_depth=1)
    assert result.status == "equivalent", result.reason
    assert result.witness["left"] and not result.witness["right"]
    assert compare_by_certified_contractions(left,right,max_attempts=0).status == "unknown"


def test_different_component_genera_are_a_handlebody_obstruction():
    tree = clean_graph()
    circle = graph_from_paths([("x","x",[(0,0,0),(1,0,0),(0,1,0),(0,0,0)])])
    assert compare_by_certified_contractions(tree,circle).status == "inequivalent"


def test_same_abstract_self_loop_unknot_and_trefoil_never_merged():
    count = 36
    circle = [(math.cos(2*math.pi*i/count),math.sin(2*math.pi*i/count),0) for i in range(count)]
    trefoil = []
    for i in range(count):
        t = 2*math.pi*i/count
        trefoil.append(((2+math.cos(3*t))*math.cos(2*t), (2+math.cos(3*t))*math.sin(2*t),math.sin(3*t)))
    circle.append(circle[0])
    trefoil.append(trefoil[0])
    left,right = graph_from_paths([("x","x",circle)]),graph_from_paths([("x","x",trefoil)])
    assert nx.is_isomorphic(left,right)
    result = compare_by_certified_contractions(left,right,max_depth=1)
    assert result.status == "unknown", result.reason


def test_same_abstract_two_loops_unlink_and_hopf_never_merged():
    count = 32
    first = [(math.cos(2*math.pi*i/count),math.sin(2*math.pi*i/count),0) for i in range(count)]
    linked = [(1+math.cos(2*math.pi*i/count),0,math.sin(2*math.pi*i/count)) for i in range(count)]
    unlinked = [(4+x,y,z) for x,y,z in linked]
    for points in (first,linked,unlinked):
        points.append(points[0])
    left = graph_from_paths([("x","x",first),("y","y",linked)])
    right = graph_from_paths([("x","x",first),("y","y",unlinked)])
    assert nx.is_isomorphic(left,right)
    assert compare_by_certified_contractions(left,right,max_depth=1).status == "unknown"
