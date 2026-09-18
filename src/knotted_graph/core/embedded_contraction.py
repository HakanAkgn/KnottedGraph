"""Conservative, exact PL edge-contraction witnesses (experimental).

Only a straight non-loop edge and disjoint swept fan triangles are supported.
Coordinates are interpreted as exact rationals, including the exact binary
values of supplied floats. No coordinate snapping or polyline simplification
is performed. This certifies the supplied PL graph, not reconstruction from a
volume. Equality of graph incidence or of Yamada polynomials is not used as a
certificate of embedding equivalence.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import combinations
import math
from numbers import Integral, Real
from typing import Any

import networkx as nx


Point = tuple[Fraction, Fraction, Fraction]


@dataclass
class ContractionResult:
    status: str
    reason: str
    graph: nx.MultiGraph | None = None
    witness: dict[str, Any] = field(default_factory=dict)
    predicate_count: int = 0


@dataclass
class ComparisonResult:
    status: str
    reason: str
    witness: dict[str, Any] = field(default_factory=dict)
    explored_states: int = 0


class _Limit(Exception):
    pass


@dataclass
class _Budget:
    limit: int
    used: int = 0

    def tick(self):
        self.used += 1
        if self.used > self.limit:
            raise _Limit


@dataclass(frozen=True)
class _Segment:
    edge: int
    index: int
    p: Point
    q: Point
    ta: Any
    tb: Any


def _point(value) -> Point:
    if len(value) != 3:
        raise ValueError("point must have three coordinates")
    values = []
    for item in value:
        if isinstance(item, Fraction):
            values.append(item)
        elif isinstance(item, Integral):
            values.append(Fraction(int(item)))
        else:
            if not isinstance(item, Real):
                raise ValueError("unsupported coordinate type; expected integers or real floating-point values")
            number = float(item)
            if not math.isfinite(number):
                raise ValueError("nonfinite coordinate")
            values.append(Fraction(number))
    return tuple(values)


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _scale(t, a):
    return tuple(t * x for x in a)


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _zero(a):
    return all(x == 0 for x in a)


def _in_box(p, points):
    return all(min(q[k] for q in points) <= p[k] <= max(q[k] for q in points) for k in range(3))


def _boxes_overlap(first, second):
    return all(max(min(p[k] for p in first), min(p[k] for p in second)) <=
               min(max(p[k] for p in first), max(p[k] for p in second)) for k in range(3))


def _on_segment(p, a, b):
    return _in_box(p, (a, b)) and _zero(_cross(_sub(p, a), _sub(b, a)))


def _segment_intersection(p, q, a, b):
    """Return zero, one or two endpoints of the exact intersection interval."""
    if not _boxes_overlap((p, q), (a, b)):
        return ()
    d, e, offset = _sub(q, p), _sub(b, a), _sub(a, p)
    n = _cross(d, e)
    if not _zero(n):
        if _dot(offset, n) != 0:
            return ()
        nn = _dot(n, n)
        t = _dot(_cross(offset, e), n) / nn
        s = _dot(_cross(offset, d), n) / nn
        if 0 <= t <= 1 and 0 <= s <= 1:
            return (_add(p, _scale(t, d)),)
        return ()
    if not _zero(_cross(offset, d)):
        return ()
    axis = next(i for i in range(3) if d[i] != 0)
    first, last = sorted(((a[axis]-p[axis])/d[axis], (b[axis]-p[axis])/d[axis]))
    lo, hi = max(Fraction(0), first), min(Fraction(1), last)
    if lo > hi:
        return ()
    start = _add(p, _scale(lo, d))
    return (start,) if lo == hi else (start, _add(p, _scale(hi, d)))


def _triangle_segment(triangle, p, q):
    """Intersection of a nondegenerate closed triangle with a closed segment."""
    if not _boxes_overlap(triangle, (p, q)):
        return ()
    a, b, c = triangle
    n = _cross(_sub(b, a), _sub(c, a))
    dp, dq = _dot(n, _sub(p, a)), _dot(n, _sub(q, a))
    direction = _sub(q, p)
    if dp != dq:
        t = dp / (dp - dq)
        if not 0 <= t <= 1:
            return ()
        x = _add(p, _scale(t, direction))
        if all(_dot(n, _cross(_sub(y, x0), _sub(x, x0))) >= 0
               for x0, y in ((a, b), (b, c), (c, a))):
            return (x,)
        return ()
    if dp != 0:
        return ()
    lo, hi = Fraction(0), Fraction(1)
    for x0, y in ((a, b), (b, c), (c, a)):
        edge = _sub(y, x0)
        value = _dot(n, _cross(edge, _sub(p, x0)))
        slope = _dot(n, _cross(edge, direction))
        if slope == 0:
            if value < 0:
                return ()
        elif slope > 0:
            lo = max(lo, -value / slope)
        else:
            hi = min(hi, -value / slope)
        if lo > hi:
            return ()
    start = _add(p, _scale(lo, direction))
    return (start,) if lo == hi else (start, _add(p, _scale(hi, direction)))


def _extract(graph, max_segments):
    if not isinstance(graph, nx.MultiGraph) or graph.is_directed():
        raise ValueError("expected an undirected networkx.MultiGraph")
    positions = {node: _point(data["pos"]) for node, data in graph.nodes(data=True)}
    if len(set(positions.values())) != len(positions):
        raise ValueError("distinct graph vertices have coincident positions")
    edges, segments = [], []
    for edge_id, (u, v, key, data) in enumerate(graph.edges(keys=True, data=True)):
        values = data.get("pts")
        if values is None:
            if u == v:
                raise ValueError("a self-loop needs an explicit closed polyline")
            values = (positions[u], positions[v])
        points = [_point(p) for p in values]
        if len(points) < 2:
            raise ValueError("an edge needs at least two samples")
        if (points[0], points[-1]) == (positions[v], positions[u]):
            points.reverse()
        if (points[0], points[-1]) != (positions[u], positions[v]):
            raise ValueError("polyline endpoints do not exactly equal vertex coordinates")
        edges.append((u, v, key, points, data))
        for index, (p, q) in enumerate(zip(points, points[1:])):
            if p == q:
                raise ValueError("zero-length segment")
            ta = ("node", u) if index == 0 else ("sample", edge_id, index)
            tb = ("node", v) if index == len(points)-2 else ("sample", edge_id, index+1)
            segments.append(_Segment(edge_id, index, p, q, ta, tb))
            if len(segments) > max_segments:
                raise _Limit
    return positions, edges, segments


def _validate(positions, segments, budget):
    for first, second in combinations(segments, 2):
        if not _boxes_overlap((first.p, first.q), (second.p, second.q)):
            continue
        budget.tick()
        intersection = _segment_intersection(first.p, first.q, second.p, second.q)
        allowed = {p for token, p in ((first.ta, first.p), (first.tb, first.q))
                   if token in (second.ta, second.tb)}
        if len(intersection) > 1 or any(p not in allowed for p in intersection):
            raise ValueError("segments intersect outside a prescribed common endpoint")
    for node, p in positions.items():
        for segment in segments:
            if ("node", node) in (segment.ta, segment.tb):
                continue
            if not _in_box(p, (segment.p, segment.q)):
                continue
            budget.tick()
            if _on_segment(p, segment.p, segment.q):
                raise ValueError("a graph vertex lies on a nonincident segment")


def certify_straight_edge_contraction(
    graph: nx.MultiGraph, edge: tuple[Any, Any, Any], *,
    max_segments: int = 2000, max_predicates: int = 200000,
) -> ContractionResult:
    """Contract edge[0] into edge[1] only when the exact empty-fan test passes.

    `certified` means regular-neighborhood equivalence of the supplied PL
    input/output. `unknown` means this move was not certified, not that the
    handlebodies differ. Invalid initial geometry receives `invalid_input`.
    No input graph or data arrays are mutated.
    """
    budget = _Budget(max_predicates)
    try:
        positions, edges, segments = _extract(graph, max_segments)
        _validate(positions, segments, budget)
    except _Limit:
        return ContractionResult("unknown", "input_validation_budget", predicate_count=budget.used)
    except (ValueError, KeyError, TypeError, OverflowError) as exc:
        return ContractionResult("invalid_input", str(exc), predicate_count=budget.used)
    u, v, selected_key = edge
    if u == v:
        return ContractionResult("unknown", "loop_contraction_forbidden", predicate_count=budget.used)
    selected = [i for i, (a,b,k,_,_) in enumerate(edges)
                if k == selected_key and ((a,b)==(u,v) or (a,b)==(v,u))]
    if len(selected) != 1:
        return ContractionResult("unknown", "selected_edge_missing", predicate_count=budget.used)
    selected_id = selected[0]
    if len(edges[selected_id][3]) != 2:
        return ContractionResult("unknown", "selected_edge_not_single_straight_segment", predicate_count=budget.used)
    start, end = positions[u], positions[v]
    moving = [s for s in segments if s.edge != selected_id and ("node",u) in (s.ta,s.tb)]
    moving_ids = {(s.edge, s.index) for s in moving}
    fixed = [s for s in segments if s.edge != selected_id and (s.edge,s.index) not in moving_ids]
    fans = []
    for s in moving:
        anchor, token = (s.q,s.tb) if s.ta == ("node",u) else (s.p,s.ta)
        triangle = (start,end,anchor)
        if _zero(_cross(_sub(end,start),_sub(anchor,start))):
            return ContractionResult("unknown", "degenerate_swept_triangle", predicate_count=budget.used)
        fans.append((triangle,token))
    try:
        for triangle, anchor_token in fans:
            for s in fixed:
                if not _boxes_overlap(triangle,(s.p,s.q)):
                    continue
                budget.tick()
                hit = _triangle_segment(triangle,s.p,s.q)
                allowed = {p for token,p in ((s.ta,s.p),(s.tb,s.q))
                           if token in (anchor_token,("node",v))}
                if len(hit)>1 or any(p not in allowed for p in hit):
                    return ContractionResult("unknown", "swept_triangle_meets_fixed_graph", predicate_count=budget.used)
            for node,p in positions.items():
                if node in (u,v) or anchor_token == ("node",node):
                    continue
                budget.tick()
                if _triangle_segment(triangle,p,p):
                    return ContractionResult("unknown", "swept_triangle_meets_other_vertex", predicate_count=budget.used)
        for (left,_),(right,_) in combinations(fans,2):
            for triangle,other in ((left,right),(right,left)):
                for a,b in ((other[0],other[1]),(other[1],other[2]),(other[2],other[0])):
                    budget.tick()
                    if any(not _on_segment(p,start,end) for p in _triangle_segment(triangle,a,b)):
                        return ContractionResult("unknown", "swept_fans_overlap_away_from_edge", predicate_count=budget.used)
        result = nx.MultiGraph()
        result.graph.update(deepcopy(graph.graph))
        for node,data in graph.nodes(data=True):
            if node != u:
                result.add_node(node,**deepcopy(data))
        result.nodes[v].setdefault("contraction_node_provenance",[])
        result.nodes[v]["contraction_node_provenance"].append({"node":u,"attributes":deepcopy(graph.nodes[u])})
        for index,(a,b,key,points,data) in enumerate(edges):
            if index == selected_id:
                continue
            updated = list(points)
            if a == u:
                updated[0] = end
            if b == u:
                updated[-1] = end
            attrs = deepcopy(data)
            # Every output coordinate already appeared in the input; float
            # conversion is exact for ordinary float input (no new rounding).
            if any(Fraction(float(x)) != x for p in updated for x in p):
                return ContractionResult("unknown", "output_coordinate_rounding", predicate_count=budget.used)
            attrs["pts"] = [tuple(float(x) for x in p) for p in updated]
            attrs.setdefault("contraction_edge_provenance",[])
            attrs["contraction_edge_provenance"].append((a,b,key))
            na,nb = (v if a==u else a),(v if b==u else b)
            new_key = key if not result.has_edge(na,nb,key) else ("contracted",index,key)
            while result.has_edge(na,nb,new_key):
                new_key = ("contracted",new_key)
            result.add_edge(na,nb,key=new_key,**attrs)
        # Recheck the actual emitted coordinates. The returned witness is about
        # this exact stored graph even if a caller supplies rational coordinates.
        output_positions,_,output_segments = _extract(result,max_segments)
        if output_positions != {node:p for node,p in positions.items() if node!=u}:
            return ContractionResult("unknown", "output_coordinate_rounding", predicate_count=budget.used)
        input_points = {pt for entry in edges for pt in entry[3]}
        for _,_,_,points,_ in _extract(result,max_segments)[1]:
            if any(p not in input_points for p in points):
                return ContractionResult("unknown", "output_coordinate_rounding", predicate_count=budget.used)
        _validate(output_positions,output_segments,budget)
    except _Limit:
        return ContractionResult("unknown", "move_predicate_budget", predicate_count=budget.used)
    except (ValueError,KeyError,TypeError) as exc:
        return ContractionResult("unknown", "output_validation:"+str(exc), predicate_count=budget.used)
    witness = {"method":"exact_rational_empty_fan_v1", "edge":edge,
               "motion":"u(t)=(1-t)u+t*v; move each incident first segment",
               "swept_triangles":[[[str(x) for x in p] for p in tri] for tri,_ in fans],
               "preserves":"ambient-isotopy class of a regular neighborhood",
               "max_degree_after":max((d for _,d in result.degree()),default=0),
               "input_coordinates":"exact binary float values; no snapping",
               "input_fingerprint":embedded_fingerprint(graph),
               "output_fingerprint":embedded_fingerprint(result)}
    return ContractionResult("certified","exact_empty_fan",result,witness,budget.used)


def embedded_fingerprint(graph):
    """An exact geometric signature, not a topological canonical form."""
    positions,edges,_ = _extract(graph,10**9)
    nodes = tuple(sorted(positions.values()))
    paths = tuple(sorted(min(tuple(points),tuple(reversed(points))) for _,_,_,points,_ in edges))
    return repr((nodes,paths))


def _component_genera(graph):
    return sorted(graph.subgraph(c).number_of_edges()-len(c)+1 for c in nx.connected_components(graph))


def compare_by_certified_contractions(
    left, right, *, max_depth=2, max_states=32, max_attempts=128,
    max_segments=2000, max_predicates=200000,
) -> ComparisonResult:
    """Bounded two-sided search; terminal match requires exact embedded geometry.

    No general isotopy solver is attempted. Exhaustion always returns unknown.
    Each primitive has its own predicate budget; total work is additionally
    bounded by max_attempts. Results retain a replayable sequence on each side.
    """
    for graph in (left,right):
        budget = _Budget(max_predicates)
        try:
            positions,_,segments = _extract(graph,max_segments)
            _validate(positions,segments,budget)
        except _Limit:
            return ComparisonResult("unknown","input_validation_budget")
        except (ValueError,KeyError,TypeError,OverflowError) as exc:
            return ComparisonResult("unknown","invalid_input:"+str(exc))
    if _component_genera(left) != _component_genera(right):
        return ComparisonResult("inequivalent","different_component_genus_multisets")
    levels = [{embedded_fingerprint(left):(left,[])},{embedded_fingerprint(right):(right,[])}]
    queues = [[(left,[])],[(right,[])]]
    attempts = states = 0
    reasons = {}
    for depth in range(max_depth+1):
        common = levels[0].keys() & levels[1].keys()
        if common:
            key = next(iter(common))
            return ComparisonResult("equivalent","certified_contractions_to_identical_PL_graph",
                                    {"left":levels[0][key][1],"right":levels[1][key][1],
                                     "terminal_fingerprint":key},states)
        if depth == max_depth:
            break
        next_queues = [[],[]]
        for side in (0,1):
            for graph,chain in queues[side]:
                states += 1
                if states > max_states:
                    return ComparisonResult("unknown","state_budget",{"rejections":reasons},states)
                for u,v,key,data in graph.edges(keys=True,data=True):
                    # Unsupported edges cannot yield this move's certificate;
                    # avoid repeating whole-graph validation for them.
                    if u == v:
                        reasons["loop_contraction_forbidden"] = reasons.get("loop_contraction_forbidden",0)+1
                        continue
                    if data.get("pts") is not None and len(data["pts"]) != 2:
                        reason = "selected_edge_not_single_straight_segment"
                        reasons[reason] = reasons.get(reason,0)+1
                        continue
                    for edge in ((u,v,key),(v,u,key)):
                        attempts += 1
                        if attempts > max_attempts:
                            return ComparisonResult("unknown","attempt_budget",{"rejections":reasons},states)
                        result = certify_straight_edge_contraction(graph,edge,max_segments=max_segments,max_predicates=max_predicates)
                        if result.status != "certified":
                            reasons[result.reason] = reasons.get(result.reason,0)+1
                            continue
                        fingerprint = embedded_fingerprint(result.graph)
                        if fingerprint not in levels[side]:
                            item = (result.graph,chain+[result.witness])
                            levels[side][fingerprint] = item
                            next_queues[side].append(item)
                            if fingerprint in levels[1-side]:
                                return ComparisonResult("equivalent","certified_contractions_to_identical_PL_graph",
                                    {"left":levels[0][fingerprint][1],"right":levels[1][fingerprint][1],
                                     "terminal_fingerprint":fingerprint},states)
        queues = next_queues
    return ComparisonResult("unknown","no_witness_in_bounded_move_set",{"rejections":reasons,"attempts":attempts},states)
