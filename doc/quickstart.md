# Quick Start

Follow a **3D volume → embedded spatial graph → planar diagram → PD code and
Yamada polynomial**. The sample is a tubular solid built around a trefoil
with an extra connecting strand: extracting its skeleton produces a graph
with two branching vertices and three edges.

```{figure} assets/site_figures/quickstart-volume.png
:alt: The surface of a tubular volume, its extracted spatial graph, and its planar diagram with numbered arcs
:target: _images/quickstart-volume.png

Left: the surface encloses the input volume. Middle: skeletonization extracts
the spatial graph; red dots mark its two vertices. Right: projection resolves
four crossings, with gaps marking the underpassing strands. The arc numbers
correspond to the PD code below.
```

Computed from the displayed volume:

```text
Graph: 2 nodes, 3 edges
Projection crossings: 4
PD code: V[7,5,0];V[4,6,10];X[2,9,1,8];X[9,0,10,1];X[8,3,7,2];X[6,4,5,3]
Normalized Yamada: Y**12 - Y**8 - Y**6 - Y**4 - Y**3 - Y**2 - Y - 1
```

## 1. Load the volume

This example uses the 0.2.0 development API. Follow {doc}`installation` to
set up the source checkout, then install the volume-processing dependency:

```bash
uv sync --extra knot-fields
```

The checkout includes
[`examples/data/quickstart-handlebody.npz`](https://github.com/HakanAkgn/KnottedGraph/blob/main/examples/data/quickstart-handlebody.npz).
You can also {download}`download the sample volume <../examples/data/quickstart-handlebody.npz>`.
Run the code from the repository root so the relative file path resolves.

```python
import numpy as np
import sympy as sp

with np.load("examples/data/quickstart-handlebody.npz", allow_pickle=False) as data:
    volume = data["volume"]
```

`volume` is a 96³ boolean array: occupied voxels describe the solid interior
bounded by the surface in the first panel. The
{download}`sample data guide <../examples/data/README.md>` gives its construction
and physical coordinate mapping.

## 2. Extract the spatial graph

Skeletonization reduces the solid to its central strands. Graph extraction
then identifies their junctions and connects them with ordered 3D polylines.

```python
from knotted_graph.core import smooth_edges
from knotted_graph.extraction import skeleton_image_to_graph, skeletonize_volume

skeleton = skeletonize_volume(volume)
graph = skeleton_image_to_graph(skeleton, max_junction_degree=3)
graph = smooth_edges(graph, epsilon=1.5)
```

`max_junction_degree=3` supplies the expected maximum junction degree for this sample.
`epsilon=1.5` simplifies the edge polylines with a tolerance of 1.5 voxels.
The extracted node positions and edge coordinates are expressed in voxel
coordinates; they form the spatial graph shown in the middle panel.

## 3. Project and compute Yamada

Projection records both the crossings and their depth order. The PD code
encodes this diagram, which is then used to evaluate the Yamada polynomial.

```python
from knotted_graph.projection import compute_yamada_polynomial

Y = sp.Symbol("Y")
result = compute_yamada_polynomial(
    graph, Y, rotation_angles=(10, 20, 30),
    normalize=True, n_jobs=1, return_result=True,
)
print(f"Graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
print("Projection crossings:", result.projection.num_crossings)
print("PD code:", result.projection.pd_code)
print("Normalized Yamada:", result.polynomial)
```

The rotation angles are in degrees. A fixed view makes this example repeatable;
you can omit `rotation_angles` to let the library choose a view.
`return_result=True` keeps the diagram alongside the polynomial, and one worker
is sufficient for this small calculation.

Run the complete example from the repository root with:

```bash
uv run --extra knot-fields python examples/volume_to_yamada.py
```

## Read the result

In the PD code, `V[...]` lists the arcs incident to a graph vertex, while
`X[...]` records the strands at a crossing in the diagram's prescribed order.
For example, `V[7,5,0]` identifies the junction of arcs 7, 5 and 0 in the
right panel. The two `V` entries describe the graph's vertices; the four `X`
entries describe its projection crossings.

`result.projection` contains the viewing angles, crossing count and PD code.
`result.polynomial` is the exact symbolic expression in `Y`.
Here `normalize=True` shifts the lowest exponent to zero. The
{doc}`user_guide/projection_yamada` guide explains these conventions and how to
inspect other projections of the same embedding.

## Explore the graph in 3D

For a rotatable view, install the visualization extra alongside volume processing:

```bash
uv sync --extra knot-fields --extra viz
```

After constructing `graph` above, run:

```python
from knotted_graph.visualization import plot_3D_graph_plotly

figure = plot_3D_graph_plotly(graph)
figure.show()
```

The {doc}`api/visualization` page covers exporting plots and viewing results
in notebooks.

## Use your own data

Replace `volume` with your own 3D boolean array. The
{doc}`user_guide/workflow_overview` guide explains volume extraction and graph
cleanup. For ordered coordinates, molecular backbones, graph tables or surface
meshes, choose a route from {doc}`user_guide/input_adapters`.

The companion CSV example loads a graph from node and edge tables:

```bash
uv run python examples/input_to_yamada.py
```

It prints the graph size, projection crossing count and Yamada polynomial.
Replace its temporary CSV creation with your own node and edge file paths to
use that route.

## Coordinate curves: a trefoil

When the 3D curve is already known, construct its graph directly from ordered
coordinates. This companion example uses the base installation.

```{figure} assets/site_figures/quickstart-trefoil.png
:alt: A three-dimensional trefoil and its planar projection with three over-under crossings

A sampled trefoil in 3D and its three-crossing projection.
```

<span id="create-a-curve-in-three-dimensions"></span>

### Build the coordinate curve

Each row contains an x, y and z coordinate. `closure="direct"` connects the
final sample back to the first.

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
curve_graph = from_coordinate_chain(points, closure="direct").graph
```

<span id="project-the-curve-and-compute-yamada"></span>

### Evaluate the curve's diagram

```python
Y = sp.Symbol("Y")
curve_result = compute_yamada_polynomial(
    curve_graph, Y, rotation_angles=(10, 20, 30),
    normalize=True, n_jobs=1, return_result=True,
)
print("Projection crossings:", curve_result.projection.num_crossings)
print("Normalized Yamada:", curve_result.polynomial)
```

```text
Projection crossings: 3
Normalized Yamada: -Y**11 + Y**9 + Y**8 + Y**7 - Y**4 - Y**3 - Y**2 - Y - 1
```

Run this coordinate example with `uv run python examples/quickstart.py`.
Choose your next workflow from {doc}`feature_status` or
{doc}`applications/index`. For installation or display help, see
{doc}`troubleshooting`.
