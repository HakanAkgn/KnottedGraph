"""Frozen pre-revision TPMS display filter; retained only to reproduce sensitivity.

This filter is not used by the revised main figure. It can overwrite errors and
depends on the arbitrary numerical ordering of labels. Do not classify with it.
"""
from collections import Counter, deque
import numpy as np


def connected_components_for_label(labels, phase_id):
    h, w = labels.shape
    seen = np.zeros(labels.shape, dtype=bool)
    components = []
    for i in range(h):
        for j in range(w):
            if seen[i, j] or int(labels[i, j]) != int(phase_id):
                continue
            queue = deque([(i, j)])
            seen[i, j] = True
            component = []
            while queue:
                ci, cj = queue.popleft()
                component.append((ci, cj))
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ni, nj = ci + di, cj + dj
                    if 0 <= ni < h and 0 <= nj < w and not seen[ni, nj] and int(labels[ni, nj]) == int(phase_id):
                        seen[ni, nj] = True
                        queue.append((ni, nj))
            components.append(component)
    return components


def stable_labels(labels, *, min_cells):
    if min_cells <= 1:
        return labels.copy(), 0
    stable = labels.copy()
    changed_cells = 0
    for _ in range(3):
        changed = False
        for phase_id in sorted(set(int(v) for v in stable.ravel())):
            for component in connected_components_for_label(stable, phase_id):
                if len(component) >= min_cells:
                    continue
                neighbor_counts = Counter()
                for i, j in component:
                    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        ni, nj = i + di, j + dj
                        if 0 <= ni < stable.shape[0] and 0 <= nj < stable.shape[1]:
                            neighbor = int(stable[ni, nj])
                            if neighbor != phase_id:
                                neighbor_counts[neighbor] += 1
                if not neighbor_counts:
                    continue
                replacement = neighbor_counts.most_common(1)[0][0]
                for i, j in component:
                    stable[i, j] = replacement
                    changed_cells += 1
                changed = True
        if not changed:
            break
    return stable, changed_cells
