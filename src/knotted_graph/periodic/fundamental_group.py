"""Fundamental groups of periodic PL regions via the sign subcomplex."""

from dataclasses import dataclass
from heapq import heappop, heappush
from itertools import permutations, product
from operator import index

import networkx as nx
import numpy as np
from numpy.typing import ArrayLike

from ._partition import OFFSETS


@dataclass(frozen=True)
class Presentation:
    """Finite group presentation with signed generator indices starting at one."""

    generator_count: int
    relators: tuple[tuple[int, ...], ...]


def _reduced_word(sequence):
    letters = []
    for letter in sequence:
        if not letter:
            continue
        if letters and letters[-1] == -letter:
            letters.pop()
        else:
            letters.append(letter)
    first, last = 0, len(letters)
    while last - first > 1 and letters[first] == -letters[last - 1]:
        first += 1
        last -= 1
    return tuple(letters[first:last])


def _inverse(word):
    return tuple(-letter for letter in reversed(word))


def _simplify(generators, initial):
    """Eliminate a generator whenever one relator solves for it directly."""
    words = {}
    occurrences = {generator: set() for generator in generators}
    pending = []
    for label, raw in enumerate(initial):
        word = _reduced_word(raw)
        if word:
            words[label] = word
            heappush(pending, (len(word), label))
            for g in set(map(abs, word)):
                occurrences[g].add(label)

    live = set(generators)
    while pending:
        length, label = heappop(pending)
        word = words.get(label)
        if word is None or len(word) != length:
            continue
        support = list(map(abs, word))
        choices = [g for g in dict.fromkeys(support) if support.count(g) == 1]
        if not choices:
            continue
        g = min(choices, key=lambda value: (len(occurrences[value]), value))
        pivot = support.index(g)
        tail = word[pivot + 1:] + word[:pivot]
        replacement = _inverse(tail) if word[pivot] > 0 else tail
        inverse = _inverse(replacement)

        for touched in sorted(occurrences[g]):
            old = words.pop(touched)
            for value in set(map(abs, old)):
                occurrences[value].remove(touched)
            if touched == label:
                continue
            new = _reduced_word(part for letter in old
                                for part in (replacement if letter == g else inverse if letter == -g else (letter,)))
            if new:
                words[touched] = new
                heappush(pending, (len(new), touched))
                for value in set(map(abs, new)):
                    occurrences[value].add(touched)
        live.remove(g)
    renumber = {g: i for i, g in enumerate(sorted(live), start=1)}
    relators = {tuple(renumber[abs(letter)] * (1 if letter > 0 else -1) for letter in word)
                for word in words.values()}
    return Presentation(len(live), tuple(sorted(relators, key=lambda word: (len(word), word))))


def periodic_region_presentations(region: ArrayLike) -> tuple[Presentation, ...]:
    """Return pi_1 presentations of all components of a periodic PL region.

    Input semantics match analyze_periodic_mask. Pass the negated mask for the
    complement. Every triangle of the induced sign subcomplex contributes a
    relation after collapsing a spanning tree in each component.
    """
    mask = np.asarray(region, dtype=bool, order="C")
    if mask.ndim != 3 or min(mask.shape) < 3:
        raise ValueError("the periodic mask must have three axes of length at least three")

    inside = mask.ravel()
    vertices = np.flatnonzero(inside)
    tree = nx.utils.UnionFind(vertices.tolist())
    vertex_ids = np.arange(mask.size).reshape(mask.shape)
    shifted = [np.roll(vertex_ids, -offset, axis=(0, 1, 2)).ravel() for offset in OFFSETS]
    edge_ids = np.zeros((8, mask.size), dtype=np.int64)
    generator_starts = {}
    for offset, translated in enumerate(shifted[1:], start=1):
        for u in vertices[inside[translated[vertices]]].tolist():
            v = int(translated[u])
            if tree[u] != tree[v]:
                tree.union(u, v)
            else:
                generator = len(generator_starts) + 1
                edge_ids[offset, u] = generator
                generator_starts[generator] = u

    roots = dict.fromkeys(tree[int(u)] for u in vertices)
    relations = {root: [] for root in roots}
    generators = {root: [] for root in roots}
    for generator, u in generator_starts.items():
        generators[tree[u]].append(generator)
    for far in range(1, 8):
        for mid in range(1, far):
            if mid & far != mid:
                continue
            a, b = shifted[mid], shifted[far]
            lower = vertices[inside[a[vertices]] & inside[b[vertices]]]
            words = np.column_stack((edge_ids[mid, lower], edge_ids[far ^ mid, a[lower]],
                                     -edge_ids[far, lower]))
            nonzero = np.any(words, axis=1)
            for u, word in zip(lower[nonzero].tolist(), words[nonzero].tolist(), strict=True):
                relations[tree[u]].append(word)
    return tuple(_simplify(generators[r], relations[r]) for r in roots)


def s3_table() -> np.ndarray:
    """Multiplication table of the symmetric group on three letters."""
    elements = tuple(permutations(range(3)))
    ids = {element: i for i, element in enumerate(elements)}
    return np.array([[ids[tuple(left[j] for j in right)] for right in elements] for left in elements])


def count_representations(presentation: Presentation, multiplication: ArrayLike, *,
                          max_assignments: int = 2_000_000) -> int:
    """Count exact homomorphisms into a small finite group with identity index zero.

    A direct bounded search is used on the reduced non-free generators. An
    oversized search raises instead of returning an estimate.
    """
    table = np.asarray(multiplication)
    if table.ndim != 2:
        raise ValueError("expected a two-dimensional multiplication table")
    order = len(table)
    ids = np.arange(order)
    if (order == 0 or table.shape != (order, order) or not np.issubdtype(table.dtype, np.integer)
        or (table < 0).any() or (table >= order).any()
        or not np.array_equal(table[0, :], ids) or not np.array_equal(table[:, 0], ids)
        or not np.all((table == 0).sum(axis=1) == 1)):
        raise ValueError("expected a finite group multiplication table with identity 0")
    if not np.all(table[table[:, :, None], ids] == table[ids[:, None, None], table]):
        raise ValueError("finite group multiplication must be associative")

    count = index(presentation.generator_count)
    support = sorted({abs(g) for word in presentation.relators for g in word})
    if count < 0 or any(g < 1 or g > count for g in support):
        raise ValueError("relators must use signed generator indices")
    assignment_count = order ** len(support)
    if assignment_count > max_assignments:
        raise ValueError(f"exact representation search needs {assignment_count} assignments")
    positions = {g: i for i, g in enumerate(support)}
    words = [[(positions[abs(g)], g < 0) for g in word] for word in presentation.relators]
    inverse = (table == 0).argmax(axis=1)

    def satisfies(values, word):
        value = 0
        for position, negated in word:
            element = values[position]
            value = table[value, inverse[element] if negated else element]
        return value == 0

    return (sum(all(satisfies(values, word) for word in words)
                for values in product(range(order), repeat=len(support)))
            * order ** (count - len(support)))
