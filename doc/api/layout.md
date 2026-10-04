# Layout

The repulsive-layout interface accepts an embedded graph and returns a
`GraphLayoutResult`. The direct graph call imports with the base package. The
`repulsion` extra adds Biopython/Plotly helpers used by protein examples and
HTML rendering:

```bash
uv sync --extra repulsion
```

Install the external C++ Repulsor solver and its native libraries using the
{doc}`../user_guide/repulsive_layout` setup. Layout updates the embedded geometry;
validate `result.graph`, then continue with projection and invariant evaluation.

## Public Python interface

```{eval-rst}
.. automodule:: knotted_graph.layout.repulsive
   :members: GraphLayoutResult, DriverConfig, SolverOptions, relax_spatial_graph
```
