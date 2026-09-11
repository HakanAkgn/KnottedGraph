# Quick Start

This page follows the shortest reproducible path from an abstract graph to an
embedded spatial graph and its Yamada polynomial. It uses only the base
installation; see {doc}`installation` if the current 0.2.0 development API is
not installed yet.

## 1. Compute a crossing-free graph directly

Start with the theta graph consisting of three parallel edges. It is small
enough to inspect immediately but, unlike a single open edge, it has no bridge
and therefore has a nonzero Yamada polynomial.

```python
import sympy as sp

from knotted_graph.core import ThetaGraph
from knotted_graph.invariants.yamada import compute_graph_yamada_polynomial

Y = sp.Symbol("Y")
theta = ThetaGraph(3)
abstract_polynomial = sp.expand(
    compute_graph_yamada_polynomial(theta, Y)
)
print(f"Upsilon(Theta_3; Y) = {abstract_polynomial}")
```

Expected output:

```text
Upsilon(Theta_3; Y) = -Y**2 - Y - 2 - 1/Y - 1/Y**2
```
Here `Y` is the polynomial variable and
\(\Upsilon(\Theta_3;Y)\) denotes the Yamada polynomial. The direct graph entry
point uses the fastest available exact backend because this abstract graph has
no crossing data to resolve.

## 2. Give the same graph a spatial embedding

An embedded graph stores a 3D position on each node and a sampled 3D polyline
in the `pts` attribute of each edge. The following construction bends the
three parallel edges into distinct planar arcs:

```python
import networkx as nx
import numpy as np

embedded_theta = nx.MultiGraph()
embedded_theta.add_node("u", pos=np.array([-2.0, 0.0, 0.0]))
embedded_theta.add_node("v", pos=np.array([2.0, 0.0, 0.0]))

curves = [
    np.array([
        [-2.0, 0.0, 0.0],
        [-1.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
    ]),
    np.array([
        [-2.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
    ]),
    np.array([
        [-2.0, 0.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [2.0, 0.0, 0.0],
    ]),
]

for points in curves:
    embedded_theta.add_edge("u", "v", pts=points)
```

## 3. Project and evaluate the embedding

Use an explicit rotation instead of randomly sampling projection directions.
This makes the tutorial deterministic. `normalize=False` preserves the same
Laurent-polynomial convention used by the direct crossing-free calculation,
and `n_jobs=1` keeps this small example single-worker and predictable.

```python
from knotted_graph.projection import compute_yamada_polynomial

result = compute_yamada_polynomial(
    embedded_theta,
    Y,
    rotation_angles=(0.0, 0.0, 0.0),
    normalize=False,
    n_jobs=1,
    method="recursive",
    return_result=True,
)

print(f"Upsilon(Theta_3; Y) = {result.polynomial}")
print(f"selected projection crossings = {result.projection.num_crossings}")
```

Expected output:

```text
Upsilon(Theta_3; Y) = -Y**2 - Y - 2 - 1/Y - 1/Y**2
selected projection crossings = 0
```

The abstract and embedded calculations agree because the selected projection
is crossing-free. A nonzero value confirms that this smoke test did not
accidentally reduce to the bridge identity.

## Run the maintained example

The complete, tested version of the code above is stored in
`examples/quickstart.py`:

```bash
uv run python examples/quickstart.py
```

## Start from files

The companion `examples/input_to_yamada.py` completes the same route from
paired node/edge CSV files. It creates a small temporary input, calls the
existing `from_spatial_graph_csv()` adapter, checks `.issues`, validates the
`pos`/`pts` embedding with `ensure_embedding()`, and prints the selected PD
code and the polynomial:

```bash
uv run python examples/input_to_yamada.py
```

Expect 2 nodes, 3 edges, zero projection crossings and the same nonzero
polynomial shown above. This route uses the base installation and no external
data or network access. The script is a source-repository example; copy it
separately if you installed only the wheel. Replace its temporary CSV creation
with your two file paths when using your own graph, and inspect input issues
before continuing. See {doc}`user_guide/input_adapters` for the CSV schema.

After this succeeds, choose the next workflow from the
{doc}`feature_status` matrix or the {doc}`user_guide/workflow_overview`. If it
fails, continue with the
{doc}`troubleshooting` guide.
