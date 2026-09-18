import json
from pathlib import Path

import networkx as nx
import numpy as np
import pytest
from shapely import LineString
import sympy as sp

from knotted_graph.projection import compute_yamada_polynomial
from knotted_graph.projection.geom import Arc
from knotted_graph.projection.pd_code import PDCode


@pytest.mark.parametrize('start', [True, False])
def test_endpoint_tangent_ignores_roundoff_duplicate_sample(start):
    # Substring at an existing sample can insert a nearly identical endpoint.
    coords = [(0.0, 0.0, 2.0), (1e-17, 0.0, 2.0), (0.0, -1.0, 2.0)]
    arc = Arc((0, 1, 0), LineString(coords if start else coords[::-1]), 'x', 0, 'x', 1)
    assert PDCode._endpoint_angle(arc, start=start) == pytest.approx(-np.pi / 2)


def test_projected_point_arc_has_no_usable_endpoint_tangent():
    arc = Arc((0, 1, 0), LineString([(1., 2., 0.), (1., 2., 1.)]), 'x', 0, 'x', 1)
    with pytest.raises(ValueError, match='no resolved endpoint direction'):
        PDCode._endpoint_angle(arc, start=True)


def test_distinct_short_segment_keeps_its_local_direction():
    arc = Arc((0, 1, 0), LineString([(0., 0., 0.), (0., 1e-10, 0.), (1., 0., 0.)]), 'x', 0, 'x', 1)
    assert PDCode._endpoint_angle(arc, start=True) == pytest.approx(np.pi / 2)


@pytest.mark.parametrize('angles', [None, (0., 0., 0.), (-32.46117974981078, 23.556464309101237, 0.)])
def test_saved_handlebody_case_preserves_yamada_across_projections(angles):
    payload = json.loads((Path(__file__).parent / 'data' / 'handlebody_random_4920_recovered.json').read_text())
    graph = nx.MultiGraph()
    for node in payload['recovered_graph']['nodes']:
        graph.add_node(node['id'], pos=np.asarray(node['pos']))
    for edge in payload['recovered_graph']['edges']:
        graph.add_edge(edge['u'], edge['v'], key=edge['key'], pts=np.asarray(edge['pts']))
    result = compute_yamada_polynomial(graph, sp.Symbol('A'), rotation_angles=angles, num_rotation_samples=6, n_jobs=1)
    expected = sp.sympify(payload['expected_yamada'])
    assert sp.expand(result - expected) == 0


@pytest.mark.parametrize('order', ['xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx', 'xyx', 'xzx', 'yxy', 'yzy', 'zxz', 'zyz', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX', 'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ'])
def test_sampling_uses_actual_hemisphere_view_directions_for_each_convention(order):
    from knotted_graph.projection.rotations import generate_isotopy_angles, get_rotation_matrix
    n = 13
    directions = np.array([get_rotation_matrix(angles, order)[2] for angles in generate_isotopy_angles(n, order)])
    i = np.arange(n)
    phi = i * np.pi * (3 - np.sqrt(5))
    z = (i + .5) / n
    r = np.sqrt(1 - z**2)
    target = np.column_stack([r * np.cos(phi), r * np.sin(phi), z])
    np.testing.assert_allclose(directions, target, atol=1e-14)
    assert np.linalg.matrix_rank(directions) == 3


def test_rotation_convention_matches_standard_intrinsic_and_extrinsic_products():
    from knotted_graph.projection.rotations import get_rotation_matrix
    a, b, c = .2, .5, .7
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    np.testing.assert_allclose(get_rotation_matrix((a, b, c), 'XYZ', True), rx @ ry @ rz)
    np.testing.assert_allclose(get_rotation_matrix((a, b, c), 'xyz', True), rz @ ry @ rx)


@pytest.mark.parametrize('kind', ['theta', 'tetrahedron'])
def test_subcubic_graph_invariant_agrees_in_each_sampled_view(kind):
    from knotted_graph.projection import sample_projections
    graph = nx.MultiGraph()
    if kind == 'theta':
        first = np.array([-2., .17, .13])
        last = np.array([2., -.31, .24])
        graph.add_node(0, pos=first)
        graph.add_node(1, pos=last)
        for middle in ([0., 1., .4], [0., 0., -.2], [0., -1., .3]):
            graph.add_edge(0, 1, pts=np.array([first, middle, last]))
    else:
        points = np.array([[-1., -1., 0.], [1., -1., .5], [1., 1., -.4], [-1., 1., .7]])
        for i, point in enumerate(points):
            graph.add_node(i, pos=point)
        for i in range(4):
            for j in range(i + 1, 4):
                graph.add_edge(i, j, pts=points[[i, j]])
    A = sp.Symbol('A')
    expected = compute_yamada_polynomial(graph, A, rotation_angles=(0., 0., 0.), n_jobs=1)
    for projection in sample_projections(graph, num_rotation_samples=6):
        actual = projection.processor.compute_yamada(A, normalize=True, n_jobs=1)
        assert sp.expand(actual - expected) == 0
