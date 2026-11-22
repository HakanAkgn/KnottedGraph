"""Spine-independent topology of a piecewise-linear band sublevel on T^3."""

from dataclasses import dataclass
from itertools import combinations, permutations

import networkx as nx
import numpy as np
from numpy.typing import ArrayLike

from ._partition import OFFSETS, graph_partition, grid_partition, winding_basis, winding_by_component


@dataclass(frozen=True)
class Region:
    """One connected occupied or complementary region of the periodic grid."""

    occupied: bool
    winding_basis: tuple[tuple[int, int, int], ...]
    euler_characteristic: int


@dataclass(frozen=True)
class Boundary:
    """One boundary surface, with incident region indices and marked winding."""

    occupied_region: int
    complementary_region: int
    genus: int
    winding_basis: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class PeriodicVolume:
    """Computed invariants for the regular PL sublevel and its complement.

    Winding vectors refer to the three marked period generators, in array-axis
    order. Region indices label a calculation and can change with the cut.
    Use has_same_invariants to compare the decorated incidence graphs.
    """

    regions: tuple[Region, ...]
    boundaries: tuple[Boundary, ...]
    occupied_betti: tuple[int, int, int, int]
    complementary_betti: tuple[int, int, int, int]

    def incidence_graph(self) -> nx.Graph:
        """Return a bipartite graph with regions and connected boundary surfaces."""
        graph = nx.Graph()
        for index, region in enumerate(self.regions):
            graph.add_node(("region", index), decoration=region)
        for index, boundary in enumerate(self.boundaries):
            node = ("boundary", index)
            graph.add_node(node, decoration=(boundary.genus, boundary.winding_basis))
            graph.add_edge(node, ("region", boundary.occupied_region))
            graph.add_edge(node, ("region", boundary.complementary_region))
        return graph

    def has_same_invariants(self, other: "PeriodicVolume") -> bool:
        """Compare the computed data modulo relabeling of regions and boundaries."""
        return (
            self.occupied_betti == other.occupied_betti
            and self.complementary_betti == other.complementary_betti
            and nx.is_isomorphic(
                self.incidence_graph(), other.incidence_graph(),
                node_match=nx.algorithms.isomorphism.categorical_node_match("decoration", None),
            )
        )


TETRAHEDRA = tuple((0, a, a | b, 7) for a, b, _ in permutations((4, 2, 1)))
TETRAHEDRON_EDGES = tuple(combinations(range(4), 2))


def _boundary_polygons(occupied, labels):
    """Yield each cut tetrahedron's incident regions and lifted crossing vertices."""
    shape = np.array(occupied.shape)
    shifted = [np.roll(occupied, -offset, axis=(0, 1, 2)) for offset in OFFSETS]
    polygon_start = 0
    for tetra in TETRAHEDRA:
        states = np.stack([shifted[index].ravel() for index in tetra])
        mixed = np.flatnonzero(states.any(axis=0) & ~states.all(axis=0))
        if not mixed.size:
            continue
        cells = np.column_stack(np.unravel_index(mixed, occupied.shape))
        signs = states[:, mixed]
        corners = OFFSETS[list(tetra)]
        inside = cells + corners[np.argmax(signs, axis=0)]
        outside = cells + corners[np.argmax(~signs, axis=0)]
        incidence = np.column_stack((labels[tuple((inside % shape).T)],
                                     labels[tuple((outside % shape).T)]))
        vertices = []
        for a, b in TETRAHEDRON_EDGES:
            crossed = np.flatnonzero(signs[a] != signs[b])
            start = cells[crossed] + corners[a]
            base = np.ravel_multi_index((start % shape).T, occupied.shape)
            # Globally identify each crossing by its positive lattice edge.
            key = 8 * base + (tetra[a] ^ tetra[b])
            vertices.append((polygon_start + crossed, key, start // shape))
        yield incidence, vertices
        polygon_start += mixed.size


def _components_of_boundary(occupied, region_labels):
    polygons = list(_boundary_polygons(occupied, region_labels))
    if not polygons:
        return ()
    incidence = np.concatenate([pair for pair, _ in polygons])
    groups = [group for _, groups in polygons for group in groups]
    polygon = np.concatenate([ids for ids, _, _ in groups])
    key = np.concatenate([key for _, key, _ in groups])
    lift = np.concatenate([lift for _, _, lift in groups])
    keys, vertex = np.unique(key, return_inverse=True)
    first = np.full(len(incidence), len(polygon), dtype=np.int64)
    np.minimum.at(first, polygon, np.arange(len(polygon)))
    reference = first[polygon]
    connected = reference != np.arange(len(polygon))
    endpoints = np.column_stack((vertex[reference[connected]], vertex[connected]))
    labels, cycles = graph_partition(len(keys), endpoints, lift[connected] - lift[reference[connected]])

    count = int(labels.max()) + 1
    bases = winding_by_component(cycles, count)
    face_labels = labels[vertex[first]]
    vertices = np.bincount(labels, minlength=count)
    faces = np.bincount(face_labels, minlength=count)
    edge_ends = np.bincount(face_labels[polygon], minlength=count)
    chi = vertices - edge_ends // 2 + faces
    if np.any((edge_ends % 2) | (chi % 2) | (chi > 2)):
        raise RuntimeError("periodic section failed the closed orientable surface check")
    representatives = np.full(count, len(incidence), dtype=np.int64)
    np.minimum.at(representatives, face_labels, np.arange(len(incidence)))
    adjacent = incidence[representatives]
    if np.any(incidence != adjacent[face_labels]):
        raise RuntimeError("inconsistent region labels on one connected section")
    return tuple(Boundary(int(a), int(b), int(1 - chi[i] // 2), bases[i])
                 for i, (a, b) in enumerate(adjacent))


def _betti(regions, occupied):
    own = [region for region in regions if region.occupied == occupied]
    complement = [region for region in regions if region.occupied != occupied]
    if not own:
        return (0, 0, 0, 0)
    if not complement:
        return (1, 3, 3, 1)
    # Via excision and Poincare-Lefschetz duality, b2 is the sum of
    # (b0(complement)-1) and the nullity of restriction H^1(T^3) -> H^1(complement).
    rank = len(winding_basis([v for region in complement for v in region.winding_basis]))
    b2 = len(complement) - 1 + 3 - rank
    b1 = len(own) + b2 - sum(region.euler_characteristic for region in own)
    return (len(own), b1, b2, 0)


def analyze_periodic_mask(occupied: ArrayLike) -> PeriodicVolume:
    """Analyze a periodic vertex mask, with a fixed Freudenthal triangulation.

    Boolean values are signs on distinct grid vertices. Interpolation uses
    values -1 and +1, respectively, inside and outside. Axes cover exactly one
    period without repeated endpoint samples. Each axis needs at least three
    samples. The result refers to this PL model.
    """
    mask = np.asarray(occupied, dtype=bool, order="C")
    if mask.ndim != 3 or min(mask.shape) < 3:
        raise ValueError("the periodic mask must have three axes of length at least three")

    labels, cycles = grid_partition(mask)
    count = int(labels.max()) + 1
    bases = winding_by_component(cycles, count)
    boundaries = _components_of_boundary(mask, labels)
    euler = [0] * count
    for boundary in boundaries:
        euler[boundary.occupied_region] += 1 - boundary.genus
        euler[boundary.complementary_region] += 1 - boundary.genus
    inside_region = np.zeros(count, dtype=bool)
    inside_region[labels[mask]] = True
    regions = tuple(Region(bool(flag), basis, chi)
                    for flag, basis, chi in zip(inside_region, bases, euler, strict=True))
    return PeriodicVolume(regions, boundaries, _betti(regions, True), _betti(regions, False))


def analyze_periodic_sublevel(energy: ArrayLike, level: float) -> PeriodicVolume:
    """Analyze the regular sublevel of a periodic scalar band's PL interpolant.

    This computes the topology at the supplied threshold for the chosen
    triangulation. Convergence to a smooth band should be checked separately.
    A threshold equal to a vertex value can be singular and is rejected.
    """
    values = np.asarray(energy)
    if (np.iscomplexobj(values) or np.iscomplexobj(level) or not np.isfinite(level)
        or not np.all(np.isfinite(values))):
        raise ValueError("energies and level must be finite real numbers")
    if np.any(values == level):
        raise ValueError("the level equals a vertex value")
    return analyze_periodic_mask(values < level)
