# Projection, PD Codes, And Yamada Polynomials

The workflow connects three objects:

1. an embedded spatial graph in three dimensions;
2. a selected regular planar diagram with over/under information; and
3. the Yamada polynomial computed from the graph/diagram data.

<div class="kg-link-row">
  <a href="../quickstart.html">Run the Quick Start</a>
  <a href="workflow_overview.html">Review the full workflow</a>
  <a href="../api/projection.html">Projection API</a>
  <a href="../api/yamada.html">Yamada API</a>
</div>

<div class="kg-wide-figure">
  <img src="../site_figures/pdcode_to_yamada.png" alt="Planar diagram and PD-code data leading to a Yamada-polynomial calculation">
</div>

## Choose the correct entry point

### Abstract graph with no crossing data

Use the direct evaluator when the intended object is an undirected abstract
Graph/MultiGraph:

```python
import sympy as sp
from knotted_graph.core import ThetaGraph
from knotted_graph.invariants.yamada import compute_graph_yamada_polynomial

Y = sp.Symbol("Y")
polynomial = compute_graph_yamada_polynomial(ThetaGraph(3), Y)
```

This route evaluates the abstract graph directly from its connectivity.

### Embedded spatial graph

Use the projection helper when node `pos` and edge `pts` geometry matter:

```python
from knotted_graph.projection import compute_yamada_polynomial

result = compute_yamada_polynomial(
    graph,
    Y,
    num_rotation_samples=16,
    n_jobs=1,
    return_result=True,
)

print(result.polynomial)
print(result.projection.rotation_angles)
print(result.projection.num_crossings)
```

`return_result=True` retains the selected projection with the polynomial,
including its rotation angles, crossing count and PD code.

## What projection selection does

Projection selection determines the diagram used in the computation. Candidate
rotations are evaluated for regularity and diagram complexity. A candidate may
fail because of overlapping projected segments, a vertex/edge degeneracy, or
an ambiguous crossing event.

For a deterministic example, supply explicit angles:

```python
result = compute_yamada_polynomial(
    graph,
    Y,
    rotation_angles=(0.0, 0.0, 0.0),
    n_jobs=1,
    return_result=True,
)
```

For explicit angles, check that the view is regular. To search for a view,
let `select_projection` sample orientations and retain its diagnostic result.

## How to read PD-code data

The planar diagram records:

- projected graph arcs;
- true graph vertices;
- transverse crossing locations;
- which arc passes over and under; and
- cyclic/local incidence information needed by the evaluator.

Crossing records describe intersections in the selected view; graph vertices
retain their identity from the spatial graph. Changing the viewing direction
may change the number and locations of crossings while leaving the embedded
graph and invariant unchanged.

If a result differs from what you expect, inspect the selected projection and
PD records first. They help locate unresolved crossings or a degenerate view
before comparing the polynomial expressions.

## Variable and normalization conventions

The user-facing examples write the polynomial as
\(\Upsilon(G;Y)\). Some backend, benchmark, or literature-comparison code uses
`A` internally. The variable name is a symbolic convention; use the same
normalization convention when comparing expressions.

Record:

- the symbolic variable;
- whether normalization was enabled;
- the selected evaluation method;
- the exact expression before/after expansion or factorization; and
- the package/backend version.

The Quick Start uses `normalize=True` for the embedded trefoil, shifting its
lowest exponent to zero. The companion CSV theta-graph example uses
`normalize=False` to compare its Laurent expression with the crossing-free
formula.

## Interpreting zero

A graph containing a bridge has zero Yamada polynomial under the implemented
convention. To interpret an exact zero,
inspect graph connectivity, bridges, leaf cleanup, and whether the intended
object was open or closed.

Also distinguish:

- a successfully evaluated exact zero;
- a failed projection with no polynomial;
- a skipped or filtered graph; and
- a missing backend or timeout.

Record these outcomes separately in tables and figures.

## Cost and worker policy

For \(c\) diagram crossings, a direct state expansion has approximately
\(3^c\) states before reductions. Geometry sampling and projection selection add
their own costs.

- Start with `n_jobs=1` and inspect the crossing count.
- Increase workers within your allocated resources.
- Record the compiled native-backend status.
- Retain failed-view diagnostics with the selected projection.

## Formula-discovery and publication notebooks

The formula-discovery notebook constructs graph families and tests proposed
identities on held-out cases. Its setup records the reference source commit.
For exact publication regeneration, enable strict mode:

```bash
export KNOTTEDGRAPH_STRICT_PUBLICATION_REGENERATION=1
```

Strict mode additionally requires the audited source revision, a clean library
tree, the expected editable checkout, and the optimized factorized native
backend. Exploratory runs report differences from the reference environment.

Start with {doc}`../quickstart`, then use the
[Advanced and Reproduction notebook](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/03_advanced_and_reproduction.ipynb).
Use formula discovery for the held-out reconstruction and verification
workflow; its setup lists the native backend and cache requirements.
