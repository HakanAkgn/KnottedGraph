"""Lift piecewise-linear level contours on a periodic plane section."""

from collections import defaultdict
from dataclasses import dataclass
from itertools import pairwise

import numpy as np
from numpy.typing import ArrayLike, NDArray


# All lifted triangles use positive horizontal, vertical and diagonal edges.
_STEPS = ((1, 0), (0, 1), (1, 1))
_TRIANGLES = (
    (((0, 0), (1, 0), (1, 1)), ((0, 0, 0), (1, 0, 1), (0, 0, 2))),
    (((0, 0), (1, 1), (0, 1)), ((0, 0, 2), (0, 1, 0), (0, 0, 1))),
)


@dataclass(frozen=True)
class SectionContour:
    """One lifted contour, in fractional coordinates of the section basis.

    The final point is the first point plus winding. For nonzero winding
    the trajectory is open in the repeated zone and its enclosed area is
    undefined. Winding orientation depends on the traversal.
    """

    winding: tuple[int, int]
    points: NDArray[np.float64]

    @property
    def area(self) -> float | None:
        """Area in units of the section parallelogram, for closed lifts."""
        if self.winding != (0, 0):
            return None
        x, y = (self.points - self.points[0]).T
        return float(abs(x[:-1] @ y[1:] - x[1:] @ y[:-1]) / 2)


def periodic_plane_contours(energy: ArrayLike, level: float) -> tuple[SectionContour, ...]:
    """Trace the regular PL level of a periodic plane sampled on a 2-D grid.

    Each axis omits its repeated endpoint and must have at least three
    samples. Coordinates are fractional, so use the area of the chosen
    reciprocal section parallelogram to convert closed-orbit areas.
    """
    values = np.asarray(energy)
    if values.ndim != 2 or min(values.shape) < 3:
        raise ValueError("the section must have two axes of length at least three")
    if (np.iscomplexobj(values) or np.iscomplexobj(level) or not np.isfinite(level)
        or not np.all(np.isfinite(values))):
        raise ValueError("energies and level must be finite real numbers")
    if np.any(values == level):
        raise ValueError("the level equals a vertex value")

    ni, nj = values.shape
    grid = np.array(values.shape)
    sign = values < level
    neighbors = defaultdict(list)
    for corners, edges in _TRIANGLES:
        states = [np.roll(sign, (-di, -dj), axis=(0, 1)) for di, dj in corners]
        changes = [first != second for first, second in pairwise(states + states[:1])]
        for i, j in np.argwhere(changes[0] | changes[1] | changes[2]):
            crossed = [(int(i + di), int(j + dj), k)
                       for change, (di, dj, k) in zip(changes, edges, strict=True)
                       if change[i, j]]
            (ai, aj, ak), (bi, bj, bk) = crossed
            a = (ai % ni, aj % nj, ak)
            b = (bi % ni, bj % nj, bk)
            deck = (bi // ni - ai // ni, bj // nj - aj // nj)
            neighbors[a].append((b, deck))
            neighbors[b].append((a, (-deck[0], -deck[1])))
    if any(len(edges) != 2 for edges in neighbors.values()):
        raise RuntimeError("unexpected vertex degree in periodic plane contour")

    def vertex_position(crossing):
        x, y, edge = crossing
        dx, dy = _STEPS[edge]
        left = float(values[x, y])
        right = float(values[(x + dx) % ni, (y + dy) % nj])
        a, b = abs(left - float(level)), abs(right - float(level))
        if not np.isfinite(a + b):
            scale = max(abs(left), abs(right), abs(float(level)))
            a, b = abs(left / scale - level / scale), abs(right / scale - level / scale)
        fraction = a / (a + b)
        return np.array([x + fraction * dx, y + fraction * dy]) / grid

    contours = []
    seen = set()
    for start in neighbors:
        if start in seen:
            continue
        initial = vertex_position(start)
        points = [initial]
        tx = ty = 0
        previous, current = None, start
        while True:
            seen.add(current)
            following, (dx, dy) = next(edge for edge in neighbors[current] if edge[0] != previous)
            tx, ty = tx + dx, ty + dy
            if following == start:
                break
            points.append(vertex_position(following) + (tx, ty))
            previous, current = current, following
        points.append(initial + (tx, ty))
        coords = np.asarray(points)
        coords.setflags(write=False)
        contours.append(SectionContour((tx, ty), coords))
    return tuple(contours)
