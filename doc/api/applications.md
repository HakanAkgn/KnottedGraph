# Applications

Application APIs assemble domain-specific models around the graph, projection
and invariant core. Use the {doc}`../feature_status` reference to choose the
input object, dependencies and public call for your workflow.

If you already have an embedded `networkx.MultiGraph`, continue directly with
the core graph and projection tools.

## Mathematical graph families

```{eval-rst}
.. automodule:: knotted_graph.applications.mathematical
   :members: graph_family_names, build_graph_case
```

## Analytic knot-field deformations

`KnotFunction` and `KnotFunctionPath` are constructed from the base
{mod}`knotted_graph.inputs` API. Running a sampled level-set deformation needs
the `knot-fields` extra because each cell extracts a 3-D level set and spatial
graph. The returned result records each cell's phase signature and any errors.

```{eval-rst}
.. automodule:: knotted_graph.applications.knot_deformation
   :members: KnotDeformationRecord, KnotDeformationScan, KnotDeformationScanResult
```

## Unified Yamada phase maps

`make_yamada_phase_map` provides one finite-grid interface for analytic knot
fields, nodal Bloch-vector paths, and in-memory material Hamiltonians. Install
`knot-fields` for knot sources or `nodal` for Hamiltonian/material sources.
Supply the in-memory source object expected by the chosen application.

The default polynomial variable is `Y`. With `force_genus_zero_vertex=True`,
each closed spherical boundary component becomes an isolated vertex: a ball
has one, while a shell has two. The volume helpers use 26-connectivity for
occupied voxels and dual 6-connectivity for the complement; their results refer
to the sampled volume. Running them needs scikit-image, included in `nodal` and
`knot-fields`. `resolve_volume_mask` returns the resolved mask and removal/fill
diagnostics; `boundary_filling_groups` returns outer and nested fillings.

```{eval-rst}
.. automodule:: knotted_graph.applications.phase_maps
   :members: YamadaPhaseRecord, YamadaPhaseMapResult, VolumeTopology, MaterialBandEnergySurface, align_material_hamiltonians, pad_material_hamiltonian, make_yamada_phase_map, enclosed_void_masks, resolve_volume_mask, volume_topology, boundary_filling_groups
```

## Saved material examples

These helpers load saved CSV/JSON phase-map records for inspection and raw
plotting with the base installation.
New scans need the optional scientific dependencies and explicit compute
resources described in {doc}`../applications/material_phase_maps`.
Use the public functions listed below. Underscore-prefixed research engines
support these workflows internally.

```{eval-rst}
.. automodule:: knotted_graph.applications.phase_map_examples
   :members: PhaseMapData, load_phase_map, read_phase_map_records, plot_phase_map
```

## Nodal and material workflows

Executing `NodalSkeleton` or `MaterialFermiSurface` requires the `nodal` extra;
the listed material symbolic-Hamiltonian constructors themselves use the base
SymPy stack. These workflows operate on in-memory SymPy Hamiltonians or Bloch
vectors. Prepare file-based models in that representation before using them.

```{eval-rst}
.. automodule:: knotted_graph.applications.nodal.models
   :members:

.. automodule:: knotted_graph.applications.nodal.skeleton
   :members: NodalSkeleton

.. automodule:: knotted_graph.applications.materials
   :members: H_Ti3Al_sympy, H_D6_sympy, H_YH3_sympy

.. py:class:: MaterialFermiSurface
   :module: knotted_graph.applications.materials

   Analyze an in-memory Hermitian multiband Hamiltonian and expose its sampled
   surface and embedded skeleton.  This class is loaded lazily and requires the
   ``nodal`` optional dependencies; see the feature-status matrix before using
   it in a workflow.
```
