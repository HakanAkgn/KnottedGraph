"""Exact local PL-manifold test for a union of closed Cartesian voxels.

A lattice vertex link consists of the occupied octants' triangles in the
boundary of an octahedron. A nonempty link must be a PL 2-sphere (all eight
triangles), or a connected PL disk. For the latter we check all link-vertex
links are paths/cycles, one boundary circle, and Euler characteristic one.
This is finite incidence arithmetic, not a sampled smoothness test. Passing
is not an analytic-source correspondence or a regular-neighborhood certificate.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations, product

import numpy as np


def connected(edges, vertices):
    vertices = set(vertices)
    if not vertices:
        return False
    adjacency = {v: set() for v in vertices}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen = {next(iter(vertices))}
    stack = list(seen)
    while stack:
        for neighbor in adjacency[stack.pop()]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return seen == vertices


def pattern_kind(pattern):
    if pattern == 0:
        return 'empty'
    triangles = [tuple(2*a + signs[a] for a in range(3))
                 for i, signs in enumerate(product((0, 1), repeat=3)) if pattern & (1 << i)]
    vertices = {v for t in triangles for v in t}
    counts = Counter(tuple(sorted(e)) for t in triangles for e in combinations(t, 2))
    if not connected(counts, vertices):
        return 'nonmanifold'
    for vertex in vertices:
        link_edges = [tuple(v for v in t if v != vertex) for t in triangles if vertex in t]
        degree = Counter(v for edge in link_edges for v in edge)
        if (not connected(link_edges, degree) or any(n not in (1, 2) for n in degree.values())
                or sum(n == 1 for n in degree.values()) not in (0, 2)):
            return 'nonmanifold'
    chi = len(vertices) - len(counts) + len(triangles)
    boundary = [e for e, n in counts.items() if n == 1]
    if pattern == 255:
        if chi != 2 or boundary:
            raise AssertionError('octahedral sphere incidence error')
        return 'sphere'
    degrees = Counter(v for e in boundary for v in e)
    if chi != 1 or not connected(boundary, degrees) or any(n != 2 for n in degrees.values()):
        return 'nonmanifold'
    return 'disk'


@lru_cache(maxsize=1)
def pattern_table():
    return tuple(pattern_kind(i) for i in range(256))


def audit_mask(mask):
    mask = np.asarray(mask)
    if mask.ndim != 3 or mask.dtype != np.bool_:
        raise ValueError('a three-dimensional Boolean closed-voxel mask is required')
    padded = np.pad(mask, 1)
    shape = tuple(n + 1 for n in mask.shape)
    patterns = np.zeros(shape, dtype=np.uint8)
    for bit, offset in enumerate(product((0, 1), repeat=3)):
        region = tuple(slice(s, s+n) for s, n in zip(offset, shape))
        patterns |= padded[region].astype(np.uint8) << bit
    histogram = np.bincount(patterns.ravel(), minlength=256)
    bad = np.array([v == 'nonmanifold' for v in pattern_table()], dtype=bool)
    bad_mask = bad[patterns]
    count = int(bad_mask.sum())
    indices = np.argwhere(bad_mask)[:32]
    return {'schema': 'knottedgraph.closed_voxel_vertex_links.v1',
            'is_pl_3_manifold_with_boundary': bool(mask.any() and count == 0),
            'empty': not bool(mask.any()), 'nonmanifold_lattice_vertices': count,
            'pattern_histogram': {str(i): int(n) for i, n in enumerate(histogram) if n},
            'first_nonmanifold_vertices': [
                {'index': p.tolist(), 'octant_pattern': int(patterns[tuple(p)])} for p in indices],
            'analytic_source_correspondence_proved': False}
