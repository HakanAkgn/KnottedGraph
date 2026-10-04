# Material Fermi Surface Fingerprints

<div class="kg-hero">
  <p class="kg-lead">The physics applications notebook explains how material Fermi surfaces are converted into knotted spatial graphs and then summarized by Yamada-polynomial fingerprints. This page uses the feature image as the visual summary of the material-to-graph-to-invariant pipeline.</p>
  <div class="kg-link-row">
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/01_physics_applications.ipynb">Open 01_physics_applications.ipynb</a>
  </div>
</div>

<div class="kg-wide-figure">
  <img src="../site_figures/feature_image1.png" alt="Material Fermi-surface fingerprints and Yamada polynomials">
</div>

## What this route accepts

Start from an in-memory symbolic Bloch Hamiltonian $H(k_x,k_y,k_z)$. Install the `nodal`
extra before constructing `MaterialFermiSurface`; install `viz` as well for the
interactive Plotly views used by the notebook.

The supplied examples include symbolic constructors for $Ti_3Al$, a
$D_6$/$TiB_2$ model, and $YH_3$. The notebook also shows a six-band
$Co_2MnGa$ construction and a template for a user-defined SymPy matrix. See the
{doc}`../feature_status` reference for the calls and dependencies used by each route.

## Data flow

```text
symbolic Hermitian Hamiltonian
    -> sample selected band gap on an explicit k-space grid
    -> threshold and inspect the complete surface
    -> skeletonize the occupied region
    -> construct and validate an embedded MultiGraph
    -> optionally clean the graph
    -> select a projection and compute an invariant
```

The continuous isosurface shows the physical geometry. Its extracted
spatial graph supplies the node and edge data used for projection and
invariant evaluation.

<a id="decisions-you-must-make-explicitly"></a>

## Choose and record the model settings

For a reproducible fingerprint, record and test:

- the momentum-space `span` and grid `dimension`;
- the selected `band_pair` and `gap_tol`;
- whether the extracted surface touches a domain boundary;
- the skeletonization and graph-cleanup settings;
- graph node/edge/degree diagnostics before and after cleanup; and
- a resolution or threshold perturbation check.

Compare graph diagnostics before and after leaf removal, short-edge
contraction or smoothing to check that the processed graph represents the
intended scientific object.

## How to read the notebook

Sections 2--7 build and diagnose non-Hermitian nodal examples, including the
surface, medial axis, graph, physical fields, and Yamada comparisons. Section 8
then applies the analogous surface-to-graph workflow to the material models.
Many examples use $300^3$ or $400^3$ grids and are publication-scale workloads;
read the Markdown and inspect cached figures before running those cells on a
proper compute resource.

Continue with {doc}`../user_guide/workflow_overview` for the reusable graph
stages and {doc}`../user_guide/projection_yamada` before interpreting a final
polynomial.
