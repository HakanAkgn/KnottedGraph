"""Exact regressions for the independent finite skein realization."""
from __future__ import annotations
import hashlib,json
import networkx as nx
import pytest
import sympy as sp
from boundary_algebra import (STATES,STATE_INDEX,IDENTITY,Y,q,partitions,is_noncrossing,
    closed_polynomial,multiply_states,left_matrix,crossing_matrix,raw_realization)

@pytest.mark.parametrize('kind,expected',[
    ('empty',1),('vertex',-1),('circle',q-1),('bridge',0),('theta',-(q-1)*(q-2)),
    ('two_loops',-(q-1)**2),('two_circles',(q-1)**2)])
def test_terminal_values(kind,expected):
    g=nx.MultiGraph()
    if kind=='vertex':g.add_node(0)
    elif kind=='circle':g.add_edge(0,0)
    elif kind=='bridge':g.add_edge(0,1)
    elif kind=='theta':g.add_edges_from([(0,1)]*3)
    elif kind=='two_loops':g.add_edges_from([(0,0)]*2)
    elif kind=='two_circles':g.add_edges_from([(0,0),(1,1)])
    assert sp.expand(closed_polynomial(g)-expected)==0

def test_fifteen_boundary_states():
    assert len(set(STATES))==15
    assert all(is_noncrossing(p) for p in STATES)
    assert all(len(b)>1 for p in STATES for b in p)

@pytest.mark.parametrize('state',STATES)
def test_two_sided_identity(state):
    assert multiply_states(IDENTITY,state)=={state:1}
    assert multiply_states(state,IDENTITY)=={state:1}

@pytest.mark.parametrize('i',[0,1])
def test_reidemeister_two(i):
    forward=crossing_matrix(i).subs(q,Y+2+1/Y)
    backward=crossing_matrix(i,-1).subs(q,Y+2+1/Y)
    assert (forward*backward-sp.eye(15)).applyfunc(sp.expand)==sp.zeros(15)

def test_braid_relation_exact():
    a=crossing_matrix(0).subs(q,Y+2+1/Y)
    b=crossing_matrix(1).subs(q,Y+2+1/Y)
    assert (a*b*a-b*a*b).applyfunc(sp.expand)==sp.zeros(15)

@pytest.mark.parametrize('word,expected',[
    ('','5fa3ca52d0fde8491c06a69e86bf352535032da7aa01758684b642f92aff49a7'),
    ('A','dbb9aed5f1e90d9b815494a8ed21c24d5c5e5667bf347aea7d1366c9e0b7f5c0'),
    ('B','4abd598fb2bc75175d8b2297dc9928635b01962a9879aa6dc537358e46cc021a'),
    ('AA','fc737cbd094177862fec0cb40c504744d44cd63ed3e7744c0b56547dbdc804ae'),
    ('AB','2ea960a8ca36a23311e883a0e6e7ec4a558ff9789dabfbdcbebadbf8ad174330')])
def test_canonical_raw_calibration(word,expected):
    left,a,b,state=raw_realization()
    for ch in reversed(word): state=((a if ch=='A' else b)*state).applyfunc(sp.expand)
    expression=sp.expand((left*state)[0])
    coeff={int(t.as_coeff_exponent(Y)[1]):int(t.as_coeff_exponent(Y)[0])
           for t in sp.Add.make_args(expression)}
    digest=hashlib.sha256(json.dumps(coeff,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert digest==expected

def test_all_basis_products_associate():
    # 225 exact matrix identities test all 3375 basis triples over Z[q].
    matrices={state:left_matrix({state:1}) for state in STATES}
    for left in STATES:
        for right in STATES:
            expected=sp.zeros(15)
            for state,weight in multiply_states(left,right).items():expected+=weight*matrices[state]
            assert (matrices[left]*matrices[right]-expected).applyfunc(sp.expand)==sp.zeros(15)

from canonical_geometry import canonical_graph,validate_geometry

@pytest.mark.parametrize('word',['','A','B','AB','BA','AAB','ABA','AAABA','ABBAAB'])
@pytest.mark.parametrize('legacy',[True,False])
def test_spatial_constructor_contacts(word,legacy):
    result=validate_geometry(canonical_graph(word,legacy_closure=legacy))
    assert len(result['distinct_unintended_points'])==(6 if legacy else 0)
    assert result['valid_PL_embedding']==(not legacy)
