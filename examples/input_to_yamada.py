"""Load two CSV files, validate a spatial graph, project it and evaluate Yamada.

Uses the base installation and temporary files. Run from a source checkout with
``uv run python examples/input_to_yamada.py``. The same script can also be run
outside the checkout against an installed wheel.
"""

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import sympy as sp

from knotted_graph.core.embedding import ensure_embedding
from knotted_graph.inputs import from_spatial_graph_csv
from knotted_graph.projection import compute_yamada_polynomial


def main() -> None:
    # Three distinct planar arcs joining the same endpoints form Theta_3.
    # Keeping the arcs as separate CSV rows preserves its parallel edges.
    with TemporaryDirectory(prefix="knotted_graph_input_") as temporary:
        directory = Path(temporary)
        nodes = directory / "nodes.csv"
        edges = directory / "edges.csv"
        with nodes.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["node_id", "x", "y", "z"])
            writer.writerows([["u", -2, 0, 0], ["v", 2, 0, 0]])
        with edges.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["edge_id", "source", "target", "points_json"])
            for index, height in enumerate([1, 0, -1]):
                points = [[-2, 0, 0], [-1, height, 0], [1, height, 0], [2, 0, 0]]
                writer.writerow([f"arc{index}", "u", "v", json.dumps(points)])

        loaded = from_spatial_graph_csv(nodes, edges, graph_id="csv-theta")
        print("Input kind:", loaded.graph.graph["input_kind"])
        print("Input issues:", loaded.issues)
        if loaded.issues:
            raise ValueError(f"Review input issues before proceeding: {loaded.issues}")
        graph = ensure_embedding(loaded.graph)
        print("Graph nodes / edges:", graph.number_of_nodes(), graph.number_of_edges())

        Y = sp.Symbol("Y")
        result = compute_yamada_polynomial(
            graph,
            Y,
            rotation_angles=(0.0, 0.0, 0.0),
            normalize=False,
            n_jobs=1,
            method="recursive",
            return_result=True,
        )
        expected = -(Y**2) - Y - 2 - Y**-1 - Y**-2
        if sp.simplify(result.polynomial - expected) != 0:
            raise RuntimeError(
                "The CSV example no longer gives the expected Theta_3 value."
            )
        if result.projection.num_crossings != 0:
            raise RuntimeError("The chosen planar projection should be crossing-free.")
        print("Selected projection crossings:", result.projection.num_crossings)
        print("PD code:", result.projection.pd_code)
        print("Upsilon(Theta_3; Y) =", result.polynomial)


if __name__ == "__main__":
    main()
