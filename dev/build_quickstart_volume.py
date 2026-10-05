"""Build the sampled tube volume used by the volume-to-Yamada example."""

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def trefoil_theta_polylines():
    """Join two trefoil arcs and an elevated connecting arc at two vertices."""
    t = np.linspace(0, 2 * np.pi, 721)
    knot = 0.42 * np.column_stack([
        (2 + np.cos(3 * t)) * np.cos(2 * t),
        (2 + np.cos(3 * t)) * np.sin(2 * t),
        np.sin(3 * t),
    ])
    split = (len(knot) - 1) // 2
    left, right = knot[0], knot[split]
    c1 = left + [0, -0.45, 2.10]
    c2 = right + [0, 0.45, 2.10]
    u = np.linspace(0, 1, 241)[:, None]
    connector = ((1 - u)**3 * left + 3 * (1 - u)**2 * u * c1
                 + 3 * (1 - u) * u**2 * c2 + u**3 * right)
    return knot[:split + 1], knot[split:][::-1], connector


def build_volume(dimension=96, radius=0.2564841814517642):
    """Sample the radius-neighborhood of the three edge polylines."""
    edges = trefoil_theta_polylines()
    extent = np.linalg.norm(np.concatenate(edges), axis=1).max() + 2.25 * radius
    origin = np.full(3, -extent)
    spacing = 2 * extent / (dimension - 1)
    volume = np.zeros((dimension,) * 3, dtype=bool)
    for edge in edges:
        for a, b in zip(edge[:-1], edge[1:]):
            lower = np.maximum(np.floor((np.minimum(a, b) - radius - origin) / spacing).astype(int), 0)
            upper = np.minimum(np.ceil((np.maximum(a, b) + radius - origin) / spacing).astype(int), dimension - 1)
            x = origin[0] + spacing * np.arange(lower[0], upper[0] + 1)[:, None, None]
            y = origin[1] + spacing * np.arange(lower[1], upper[1] + 1)[None, :, None]
            z = origin[2] + spacing * np.arange(lower[2], upper[2] + 1)[None, None, :]
            delta = b - a
            denominator = float(delta @ delta)
            if denominator <= 1e-20:
                distance_squared = (x - a[0])**2 + (y - a[1])**2 + (z - a[2])**2
            else:
                t = np.clip(((x - a[0]) * delta[0] + (y - a[1]) * delta[1]
                             + (z - a[2]) * delta[2]) / denominator, 0, 1)
                distance_squared = ((x - a[0] - t * delta[0])**2
                                    + (y - a[1] - t * delta[1])**2
                                    + (z - a[2] - t * delta[2])**2)
            slices = tuple(slice(int(lo), int(hi) + 1) for lo, hi in zip(lower, upper))
            volume[slices] |= distance_squared <= radius**2
    return volume, origin, float(spacing), float(radius)


def main():
    volume, origin, spacing, radius = build_volume()
    destination = ROOT / "examples" / "data" / "quickstart-handlebody.npz"
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, volume=volume, origin=origin, spacing=spacing,
                        radius=radius, case=np.asarray("control_trefoil_theta"))
    print(destination)


if __name__ == "__main__":
    main()
