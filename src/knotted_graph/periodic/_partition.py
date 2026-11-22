"""Compute connected components and integral winding from lifted edges."""

from itertools import product

import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

try:
    from . import _periodic_native
except ImportError:
    _periodic_native = None


OFFSETS = np.array(tuple(product((0, 1), repeat=3)))


def winding_basis(vectors):
    """Basis for the subgroup of the marked ambient Z^3, without saturation."""
    if not vectors:
        return ()
    columns = hermite_normal_form(Matrix(vectors).T).T
    return tuple(tuple(map(int, columns.row(i))) for i in range(columns.rows))


def winding_by_component(cycles, count):
    vectors = [[] for _ in range(count)]
    for label, x, y, z in cycles:
        vectors[label].append((x, y, z))
    return tuple(winding_basis(items) for items in vectors)


def _graph_partition_python(count, edges, shifts):
    adjacency = [[] for _ in range(count)]
    for (a, b), shift in zip(edges.tolist(), shifts.tolist(), strict=True):
        x, y, z = shift
        adjacency[a].append((b, x, y, z))
        adjacency[b].append((a, -x, -y, -z))
    labels = np.full(count, -1, dtype=np.int64)
    lift = [None] * count
    cycles = []
    label = 0
    for seed in range(count):
        if lift[seed] is not None:
            continue
        lift[seed] = (0, 0, 0)
        seen_cycles = set()
        stack = [seed]
        while stack:
            vertex = stack.pop()
            labels[vertex] = label
            x, y, z = lift[vertex]
            for other, dx, dy, dz in adjacency[vertex]:
                candidate = (x + dx, y + dy, z + dz)
                if lift[other] is None:
                    lift[other] = candidate
                    stack.append(other)
                elif candidate != lift[other]:
                    winding = tuple(a - b for a, b in zip(candidate, lift[other], strict=True))
                    seen_cycles.add(winding)
        cycles.extend((label, *row) for row in sorted(seen_cycles))
        label += 1
    return labels, cycles


def graph_partition(count, edges, shifts):
    edges = np.ascontiguousarray(edges, dtype=np.int64)
    shifts = np.ascontiguousarray(shifts, dtype=np.int64)
    if _periodic_native is None:
        return _graph_partition_python(count, edges, shifts)
    return _periodic_native.graph_partition(count, edges, shifts)


def _grid_partition_python(occupied):
    shape = np.array(occupied.shape)
    edges = []
    shifts = []
    for offset in OFFSETS[1:]:
        start = np.argwhere(occupied == np.roll(occupied, -offset, axis=(0, 1, 2)))
        endpoint = start + offset
        left = np.ravel_multi_index(start.T, occupied.shape)
        right = np.ravel_multi_index((endpoint % shape).T, occupied.shape)
        edges.append(np.column_stack((left, right)))
        shifts.append(endpoint // shape)
    labels, cycles = _graph_partition_python(
        occupied.size, np.concatenate(edges), np.concatenate(shifts)
    )
    return labels.reshape(occupied.shape), cycles


def grid_partition(occupied):
    if _periodic_native is None:
        return _grid_partition_python(occupied)
    return _periodic_native.grid_partition(occupied)
