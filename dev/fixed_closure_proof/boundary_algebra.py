"""Exact planar boundary-state realization for a fixed three-strand Yamada closure.

This implementation is independent of the geometric projection and production
Yamada evaluators. All graph reductions use deletion/contraction and the raw
Yamada terminal polynomial. No trained transfer coefficients are used.
"""
from __future__ import annotations
from functools import lru_cache
from collections import defaultdict
from itertools import combinations
import sympy as sp
import networkx as nx

q, Y = sp.symbols('q Y')
# Cyclic order for top labels 0,1,2 and bottom labels 3,4,5.
CYCLIC = (0, 1, 2, 5, 4, 3)
IDENTITY = ((0, 3), (1, 4), (2, 5))
BASIS_WORDS = ['', 'A', 'B', 'AA', 'AB', 'BA', 'BB',
               'AAB', 'ABA', 'ABB', 'BAA', 'BAB', 'AABA', 'ABAA', 'ABAB']


def partitions(items):
    """Enumerate set partitions, with a canonical order and no duplication."""
    if not items:
        yield ()
        return
    first, *rest = items
    for p in partitions(rest):
        yield tuple(sorted(((first,),) + p))
        for i, block in enumerate(p):
            yield tuple(sorted(p[:i] + (tuple(sorted((first,) + block)),) + p[i+1:]))


def is_noncrossing(partition):
    positions = {label: i for i, label in enumerate(CYCLIC)}
    blocks = [set(positions[v] for v in block) for block in partition]
    for a, b in combinations(blocks, 2):
        for i, j, k, l in combinations(range(6), 4):
            if (i in a and k in a and j in b and l in b) or (i in b and k in b and j in a and l in a):
                return False
    return True


STATES = tuple(sorted(set(p for p in partitions(list(range(6)))
                          if all(len(b) >= 2 for b in p) and is_noncrossing(p))))
assert len(STATES) == 15
STATE_INDEX = {p: i for i, p in enumerate(STATES)}


def state_graph(partition, labels=None, prefix='s'):
    """Pairs become arcs; larger blocks become stars with one coupon."""
    if labels is None:
        labels = {i: i for i in range(6)}
    graph = nx.MultiGraph()
    graph.add_nodes_from(labels.values())
    for j, block in enumerate(partition):
        ports = [labels[i] for i in block]
        if len(ports) == 2:
            graph.add_edge(*ports)
        else:
            center = (prefix, j)
            graph.add_edges_from((center, port) for port in ports)
    return graph


def closed_polynomial(graph):
    """Raw planar Yamada polynomial in q=Y+2+Y^-1, by its subset sum.

    R(G)=sum_{S subset E} (-1)^(|V|+|S|) q^(|S|-|V|+k(S)).
    Isolated vertices, loops and parallel edges are retained.
    """
    nodes = list(graph)
    index = {v: i for i, v in enumerate(nodes)}
    edges = [(index[u], index[v]) for u, v in graph.edges()]
    n, m = len(nodes), len(edges)
    coeffs = defaultdict(int)
    for mask in range(1 << m):
        parent = list(range(n))
        components = n
        def root(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        count = 0
        for i, (u, v) in enumerate(edges):
            if mask >> i & 1:
                count += 1
                a, b = root(u), root(v)
                if a != b:
                    parent[a] = b
                    components -= 1
        power = count - n + components
        coeffs[power] += -1 if (n + count) % 2 else 1
    return sp.Add(*(value*q**power for power, value in coeffs.items() if value))


def _contract(graph, u, v, key):
    out = graph.copy()
    out.remove_edge(u, v, key)
    result = nx.MultiGraph()
    result.add_nodes_from(node for node in out if node != v)
    result.add_edges_from((u if a == v else a, u if b == v else b)
                         for a, b in out.edges())
    return result


def reduce_open(graph, boundary=frozenset(range(6))):
    """Resolve all internal edges. Return a linear combination of 15 stars."""
    graph = graph.copy()
    factor = sp.Integer(1)
    # Closed connected components factor out. One-port components vanish under
    # every closure, since their only port is attached through a cut edge.
    for comp in list(nx.connected_components(graph)):
        ports = boundary.intersection(comp)
        if not ports:
            factor *= closed_polynomial(graph.subgraph(comp))
            graph.remove_nodes_from(comp)
        elif len(ports) == 1:
            return {}
    # Remove loops; R(G with loop)=-(q-1) R(G without that loop).
    loops = list(nx.selfloop_edges(graph, keys=True))
    for u, v, key in loops:
        graph.remove_edge(u, v, key)
        factor *= 1-q
    # Suppress internal bivalent vertices, preserving loops and parallel arcs.
    changed = True
    while changed:
        changed = False
        for node in list(graph):
            if node in boundary or graph.degree(node) != 2:
                continue
            neighbors = [v if u == node else u for u, v in graph.edges(node)]
            if len(neighbors) == 2 and node not in neighbors:
                graph.remove_node(node)
                graph.add_edge(*neighbors)
                changed = True
                break
    # Suppression may have produced a loop; recurse once to use the loop rule.
    if nx.number_of_selfloops(graph):
        return {p: sp.expand(factor*c) for p, c in reduce_open(graph, boundary).items() if c != 0}
    for u, v, key in graph.edges(keys=True):
        if u not in boundary and v not in boundary:
            deleted = graph.copy()
            deleted.remove_edge(u, v, key)
            contracted = _contract(graph, u, v, key)
            out = defaultdict(lambda: sp.Integer(0))
            for term in (deleted, contracted):
                for state, coefficient in reduce_open(term, boundary).items():
                    out[state] += coefficient
            return {state: sp.expand(factor*coefficient) for state, coefficient in out.items()
                    if sp.expand(factor*coefficient) != 0}
    # Every remaining component is now a single boundary star or an arc.
    state = tuple(sorted(tuple(sorted(boundary.intersection(comp)))
                         for comp in nx.connected_components(graph)))
    if any(len(block) < 2 for block in state):
        raise AssertionError('Unreduced singleton or closed component')
    if state not in STATE_INDEX:
        raise AssertionError(f'Nonplanar/unexpected boundary state: {state}')
    return {state: sp.expand(factor)} if factor != 0 else {}


@lru_cache(maxsize=None)
def multiply_states(left, right):
    """Stack left above right, respecting the three ordered boundary ports."""
    # External labels remain 0..5. Interface nodes are not marked boundary.
    lmap = {i: i if i < 3 else ('interface', i-3) for i in range(6)}
    rmap = {i: ('interface', i) if i < 3 else i for i in range(6)}
    graph = state_graph(left, lmap, 'left')
    lower = state_graph(right, rmap, 'right')
    graph.add_nodes_from(lower.nodes())
    graph.add_edges_from(lower.edges())  # multiset union: never overwrite parallel interface arcs
    return reduce_open(graph)


def elementary(i):
    if i not in (0, 1):
        raise ValueError('Three-strand generator index must be 0 or 1')
    other = 2 if i == 0 else 0
    cupcap = tuple(sorted(((i, i+1), (i+3, i+4), (other, other+3))))
    vertex = tuple(sorted(((i, i+1, i+3, i+4), (other, other+3))))
    return cupcap, vertex


def left_matrix(terms):
    out = sp.zeros(15)
    for j, state in enumerate(STATES):
        for left, weight in terms.items():
            for result, coefficient in multiply_states(left, state).items():
                out[STATE_INDEX[result], j] += weight*coefficient
    return out.applyfunc(sp.expand)


def crossing_matrix(i, orientation=1):
    """One fixed crossing resolved into straight, cup-cap and vertex terms."""
    cupcap, vertex = elementary(i)
    return left_matrix({IDENTITY: Y**orientation, cupcap: Y**(-orientation), vertex: sp.Integer(1)})


def exterior_graph():
    """Exactly the nine non-braid edges in the archived eight-vertex closure."""
    graph = nx.MultiGraph()
    # 0 LT,1 LMB,2 LB,3 RT,4 RMB,5 RB,6 LMC,7 RMC.
    graph.add_nodes_from(range(8))
    graph.add_edges_from([(1,6),(4,7),(0,1),(3,4),(6,2),(7,5),
                          (3,0),(7,6),(5,2)])
    return graph


def closure_vector():
    out = []
    for p in STATES:
        base = exterior_graph()
        inserted = state_graph(p, prefix='inside')
        base.add_nodes_from(inserted.nodes())
        base.add_edges_from(inserted.edges())
        out.append(closed_polynomial(base))
    return sp.Matrix([out])


def raw_realization(orientation=1):
    replacement = Y+2+1/Y
    braid = [crossing_matrix(i, orientation).subs(q, replacement).applyfunc(sp.expand)
             for i in (0,1)]
    matrices = [(b*b).applyfunc(sp.expand) for b in braid]
    left = closure_vector().subs(q,replacement).applyfunc(sp.expand)
    right = sp.zeros(15,1)
    right[STATE_INDEX[IDENTITY]] = 1
    return left, matrices[0], matrices[1], right


if __name__ == '__main__':
    import json, time
    start=time.perf_counter()
    left, ta, tb, right = raw_realization()
    print('States:',len(STATES),'seconds:',time.perf_counter()-start,flush=True)
    print('empty:',sp.factor((left*right)[0]),flush=True)
    for name, matrix in [('A',ta),('B',tb)]:
        print(name,sp.factor((left*matrix*right)[0]),flush=True)
        cube=(matrix-Y**2*sp.eye(15))*(matrix-Y**-2*sp.eye(15))*(matrix-Y**-4*sp.eye(15))
        print('cubic',cube.applyfunc(sp.expand)==sp.zeros(15),flush=True)
    print('commutator_zero',(ta*tb-tb*ta).applyfunc(sp.expand)==sp.zeros(15),flush=True)
