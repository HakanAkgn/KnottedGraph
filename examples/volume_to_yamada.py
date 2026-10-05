"""Extract a spatial graph from a tube volume and evaluate its planar diagram.

Run with ``uv run --extra knot-fields python examples/volume_to_yamada.py``.
"""

from pathlib import Path

import numpy as np
import sympy as sp

from knotted_graph.core import smooth_edges
from knotted_graph.extraction import skeleton_image_to_graph, skeletonize_volume
from knotted_graph.projection import compute_yamada_polynomial


VOLUME_PATH = Path(__file__).resolve().parent / "data" / "quickstart-handlebody.npz"


def load_volume():
    """Load the sampled trefoil-theta tube volume."""
    with np.load(VOLUME_PATH, allow_pickle=False) as sample:
        return sample["volume"]


def extract_volume_graph(volume):
    """Skeletonize the volume and simplify its edge polylines in voxel units."""
    skeleton = skeletonize_volume(volume)
    graph = skeleton_image_to_graph(skeleton, max_junction_degree=3)
    return smooth_edges(graph, epsilon=1.5)


def compute_volume_quickstart():
    """Return the input volume, extracted graph and computed projection result."""
    volume = load_volume()
    graph = extract_volume_graph(volume)
    Y = sp.Symbol("Y")
    result = compute_yamada_polynomial(
        graph, Y, rotation_angles=(10, 20, 30),
        normalize=True, n_jobs=1, return_result=True,
    )
    return volume, graph, result


def main():
    """Print the diagram and polynomial shown in the README and Quick Start."""
    _, graph, result = compute_volume_quickstart()
    print(f"Graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    print("Projection crossings:", result.projection.num_crossings)
    print("PD code:", result.projection.pd_code)
    print("Normalized Yamada:", result.polynomial)


if __name__ == "__main__":
    main()
