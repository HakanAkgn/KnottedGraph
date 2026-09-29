"""Build a 3D trefoil, project it and compute its normalized Yamada polynomial.

Run with ``uv run python examples/quickstart.py`` from the source checkout.
Only the base installation is needed.
"""

import numpy as np
import sympy as sp

from knotted_graph.inputs import from_coordinate_chain
from knotted_graph.projection import compute_yamada_polynomial


def build_trefoil():
    """Return a trefoil sampled as a closed curve in three dimensions."""
    t = np.linspace(0, 2 * np.pi, 120, endpoint=False)
    points = np.column_stack([
        (2 + np.cos(3 * t)) * np.cos(2 * t),
        (2 + np.cos(3 * t)) * np.sin(2 * t),
        np.sin(3 * t),
    ])
    # Close this known periodic curve with its final last-to-first segment.
    # The input adapter stores it as a graph with one node and one loop edge.
    return from_coordinate_chain(points, closure="direct").graph


def compute_quickstart():
    """Return the embedded graph and its polynomial/projection result."""
    graph = build_trefoil()
    Y = sp.Symbol("Y")
    result = compute_yamada_polynomial(
        graph,
        Y,
        # Euler angles are in degrees; this fixed view has three crossings.
        rotation_angles=(10, 20, 30),
        normalize=True,
        # One worker is sufficient for this small introductory example.
        n_jobs=1,
        return_result=True,
    )
    return graph, result


def main() -> None:
    """Print the result shown in the README and online Quick Start."""
    _, result = compute_quickstart()
    print("Projection crossings:", result.projection.num_crossings)
    print("Normalized Yamada:", result.polynomial)


if __name__ == "__main__":
    main()
