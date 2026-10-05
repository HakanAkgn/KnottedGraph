# Architecture

The repository is organized around generic mathematical components first, with
applications layered on top.

## Submodules

- `knotted_graph.core`: embedded graph validation, graph constructors, and graph
  simplification helpers.
- `knotted_graph.inputs`: adapters that normalize external data into core
  graph or mesh objects.
- `knotted_graph.extraction`: production sparse skeleton-image conversion.
- `knotted_graph.projection`: rotations, projection sampling, and PD-code
  objects.
- `knotted_graph.invariants.yamada`: Yamada polynomial evaluation for resolved
  graphs and diagrams.
- `knotted_graph.layout.repulsive`: optional curve-network relaxation.
- `knotted_graph.visualization`: generic graph plotting helpers.
- `knotted_graph.applications.nodal`: non-Hermitian nodal-skeleton workflow.

## Import Boundary

`import knotted_graph` must not import optional application stacks such as
PyVista, scikit-image, poly2graph, minorminer, or Plotly. Application workflows
are imported explicitly from their application namespace.

## Skeleton Extraction Policy

`knotted_graph.extraction.skeleton_image_to_graph` is the canonical
skeleton-to-graph entry point and uses the sparse extractor.
Empty image margins are cropped before foreground indexing, 26-neighbour
adjacency is generated deterministically, and returned coordinates remain in
the original global voxel frame. Use Git history and the isolated
`02_application_regression_checks.ipynb` worktree comparison to inspect
historical extraction behavior.

By default, optimized extraction preserves the zero-radius historical topology.
A caller that knows a genuine valence bound may pass `max_junction_degree` to
enable persistence-based junction repair. Without a supplied valence bound,
the extractor preserves generic higher-valence junctions.

## Yamada Evaluation Policy

`Yamada.compute()` has one diagram-level production algorithm: the exact
factorized-connectivity dynamic program implemented by
`factorized_frontier.py` and `_yamada_factorized_frontier`. This same
algorithm handles every diagram-level production evaluation.

The compiled factorized extension is required for diagram-level production
execution. A stale or incomplete installation fails with an explicit rebuild
instruction. The compact
resolved-graph kernel remains available for crossing-free graphs and
`compute_yamada_from_states()` compatibility. Exhaustive/raw evaluators are
retained only as exact validation or arithmetic-overflow safety references.

## Yamada Projection Policy

`knotted_graph.projection.compute_yamada_polynomial` is the canonical embedded
spatial-graph entry point. If `rotation_angles` is supplied, it computes exactly
that projection. Otherwise it samples `num_rotation_samples=10` rotations,
chooses the valid projection with the fewest crossings, and emits a
`RuntimeWarning` when the selected diagram has at least 10 crossings.

## Inspection-Oriented API Direction

Keep intermediate objects available for inspection and export at the
stages used in the paper figures: imported input, surface or mesh,
skeleton image, raw spatial graph, simplified spatial graph, sampled
projections, selected planar diagram, PD code, and invariant output.

Graph-to-Yamada workflows expose this pattern through
`sample_projections`, `select_projection`, and
`compute_yamada_polynomial(..., return_result=True)`. Future high-level
pipeline helpers should preserve these lower-level APIs while collecting the
same intermediate objects into a convenient result object for tutorials,
figures, and reproducible analyses.
