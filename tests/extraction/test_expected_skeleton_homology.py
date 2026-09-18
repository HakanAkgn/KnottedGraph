import networkx as nx
import numpy as np
import pytest

from knotted_graph.extraction import skeleton_image_to_graph, topology_aware_skeleton_image_to_graph, skeletonize_volume


def cycle_rank(graph):
    return graph.number_of_edges()-len(graph)+nx.number_connected_components(graph)


def two_components():
    mask=np.zeros((12,12,12),dtype=bool)
    mask[2,2,2:6]=True
    mask[8,8,2:6]=True
    return mask


@pytest.mark.parametrize('entry',[skeleton_image_to_graph,topology_aware_skeleton_image_to_graph])
def test_expected_component_counts_forwarded_and_impossible_request_fails(entry):
    graph=entry(two_components(),expected_components=2,expected_cycle_rank=0)
    assert nx.number_connected_components(graph)==2
    assert cycle_rank(graph)==0
    with pytest.raises(ValueError,match='No persistent or zero-radius'):
        entry(two_components(),expected_components=1,expected_cycle_rank=0)
    with pytest.raises(ValueError,match='No persistent or zero-radius'):
        entry(two_components(),expected_components=2,expected_cycle_rank=99)


@pytest.mark.parametrize('option,value',[
    ('expected_cycle_rank',-1),('expected_components',-1),
    ('expected_cycle_rank',1.2),('expected_components',1.0),
])
def test_invalid_expected_count_is_not_silently_coerced(option,value):
    with pytest.raises(ValueError,match=option):
        skeleton_image_to_graph(two_components(),**{option:value})


def test_96_grid_gyroid_real_cycles_survive_persistence_selection():
    pytest.importorskip('skimage')
    from skimage.measure import euler_number,label
    axis=np.linspace(-2.25*np.pi,2.25*np.pi,96)
    x,y,z=np.meshgrid(axis,axis,axis,indexing='ij')
    mask=(np.sin(x)*np.cos(y)+np.sin(y)*np.cos(z)+np.sin(z)*np.cos(x)<=0)
    mask &= x*x+y*y+z*z <= (.72*2.25*np.pi)**2
    assert label(mask,connectivity=3,return_num=True)[1]==1
    assert euler_number(mask,connectivity=3)==-10
    skeleton=skeletonize_volume(mask)
    graph=skeleton_image_to_graph(skeleton,expected_components=1,expected_cycle_rank=11)
    assert nx.number_connected_components(graph)==1
    assert cycle_rank(graph)==11


def test_guard_does_not_change_unconstrained_results():
    mask=two_components()
    plain=skeleton_image_to_graph(mask)
    explicit=skeleton_image_to_graph(mask,expected_components=None,expected_cycle_rank=None)
    assert nx.is_isomorphic(plain,explicit)
    assert all(np.array_equal(plain.nodes[n]['pos'],explicit.nodes[n]['pos']) for n in plain)
