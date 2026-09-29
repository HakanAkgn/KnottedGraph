# Quick Start

Follow a small trefoil from **3D coordinates to a planar diagram and its
Yamada polynomial**. This example uses the 0.2.0 development API and only the
base installation; follow {doc}`installation` to set up the source checkout.

```{figure} assets/site_figures/quickstart-trefoil.png
:alt: A three-dimensional trefoil and its planar projection with three over-under crossings

The same sampled curve in 3D (left) and in the projection used for the
calculation (right). Gaps in the diagram show which strand passes underneath.
```

## 1. Create a curve in three dimensions

Each row of `points` contains an x, y and z coordinate. The formulas below
sample a trefoil, a closed knotted curve. `closure="direct"` adds the short
segment from the final sample back to the first.

```python
import numpy as np
import sympy as sp

from knotted_graph.inputs import from_coordinate_chain
from knotted_graph.projection import compute_yamada_polynomial

t = np.linspace(0, 2 * np.pi, 120, endpoint=False)
points = np.column_stack([
    (2 + np.cos(3 * t)) * np.cos(2 * t),
    (2 + np.cos(3 * t)) * np.sin(2 * t),
    np.sin(3 * t),
])
graph = from_coordinate_chain(points, closure="direct").graph
```

The adapter returns an embedded graph with one node and one loop edge. That
edge stores the sampled 3D curve in `pts`, so its geometry is available to the
projection step.

## 2. Project the curve and compute Yamada

A projection shows where strands cross and uses their depth to determine
which passes over the other. The selected view below has three crossings.

```python
Y = sp.Symbol("Y")
result = compute_yamada_polynomial(
    graph, Y, rotation_angles=(10, 20, 30),
    normalize=True, n_jobs=1, return_result=True,
)
print("Projection crossings:", result.projection.num_crossings)
print("Normalized Yamada:", result.polynomial)
```

The rotation angles are in degrees. Using a fixed view makes this example
repeatable; for your own data, you can omit `rotation_angles` to let the
library choose a view. `return_result=True` keeps the diagram with the
polynomial, and `n_jobs=1` is sufficient for this small calculation.

## 3. Read the result

Expected output:

```text
Projection crossings: 3
Normalized Yamada: -Y**11 + Y**9 + Y**8 + Y**7 - Y**4 - Y**3 - Y**2 - Y - 1
```

`result.projection` contains the viewing angles, crossing count and PD code
(the encoded planar diagram). `result.polynomial` is the exact symbolic
expression in `Y`. Here `normalize=True` shifts the lowest exponent to zero.
The {doc}`user_guide/projection_yamada` guide explains these conventions and
how to inspect other projections of the same embedding.

Run the complete, maintained example from the repository root with:

```bash
uv run python examples/quickstart.py
```

## Explore the graph in 3D

For a rotatable view, install the visualization extra:

```bash
uv sync --extra viz
```

After constructing `graph` above, run:

```python
from knotted_graph.visualization import plot_3D_graph_plotly

figure = plot_3D_graph_plotly(graph)
figure.show()
```

The {doc}`api/visualization` page also covers exporting plots and viewing
results in notebooks.

## Use your own data

Replace `points` with your own ordered 3D coordinates, or choose a file reader
from {doc}`user_guide/input_adapters`. That guide explains curve closure and
how to prepare molecular backbones, graph tables and surface meshes.

For a graph with several nodes and edges, run the companion CSV example:

```bash
uv run python examples/input_to_yamada.py
```

It loads a planar theta graph from two small temporary CSV files, then prints
2 nodes, 3 edges, zero projection crossings and its Yamada polynomial. Replace
the temporary CSV creation with your own node and edge file paths to use that
route. These scripts are included in the source checkout.

Choose your next example from {doc}`feature_status` or the
{doc}`applications/index`. For installation or display help, see
{doc}`troubleshooting`.
