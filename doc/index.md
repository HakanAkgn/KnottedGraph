# KnottedGraph

<div class="kg-hero">
  <p class="kg-lead"><strong>KnottedGraph</strong> studies three-dimensional geometry through embedded spatial graphs. Extract a graph from a sampled volume, resolve its planar projection, and compute its PD code and Yamada polynomial. Explore the complete example below, application tutorials, or the paper's figures and saved data.</p>
  <div class="kg-link-row">
    <a href="paper_results.html">Paper Figures &amp; Data</a>
    <a href="benchmarks.html">Benchmarks</a>
    <a href="sanity_checks.html">Sanity Checks</a>
    <a href="citing.html">Cite Our Work</a>
    <a href="installation.html">Install</a>
    <a href="quickstart.html">Quick Start</a>
    <a href="feature_status.html">Choose A Workflow</a>
    <a href="user_guide/index.html">User Guide</a>
    <a href="applications/index.html">Applications</a>
    <a href="troubleshooting.html">Troubleshooting</a>
  </div>
</div>

## Start with a 3D example

The {doc}`quickstart` starts with a tube volume, extracts a graph with two
vertices and three edges, and computes its planar diagram, PD code and exact
Yamada polynomial. The first panel shows the volume's surface; the second
shows the extracted graph; the third labels the arcs used in the PD code.

```{figure} assets/site_figures/quickstart-volume.png
:alt: Surface of a branched tube volume, its extracted spatial graph and its four-crossing planar diagram
:target: _images/quickstart-volume.png

{doc}`Open the step-by-step Quick Start <quickstart>`.
```

```text
PD code: V[7,5,0];V[4,6,10];X[2,9,1,8];X[9,0,10,1];X[8,3,7,2];X[6,4,5,3]
Normalized Yamada: Y**12 - Y**8 - Y**6 - Y**4 - Y**3 - Y**2 - Y - 1
```

## Read the paper, inspect the evidence

The {doc}`paper_results` guide places the published figures beside their source
data and notebooks. The {doc}`benchmarks` page explains the 500-crossing scaling
measurements and 4,400-case construction benchmark, with direct CSV downloads.
Figure 4 and Supplementary Figure 10 have downloadable coefficient
tables, exact matrices and a complete data/source archive beside their figures.
Use {doc}`sanity_checks` for the small correctness checks and their expected
outcomes. Please {doc}`cite the software and application papers <citing>` when
using this work.

**{download}`Download Figure 4 / Supplementary Figure 10 data and sources <assets/data/figure4-suppfig10-source-data.zip>`** ·
[Browse figure panels and individual files](paper_results.md#figure-4-discovering-and-testing-family-laws)

## New here?

1. Follow {doc}`installation` for the 0.2 development API.
2. Run {doc}`quickstart` to extract, project and analyse a small 3D volume.
3. Choose your real starting object from {doc}`feature_status`.
4. Follow {doc}`user_guide/workflow_overview` to connect the steps for your data.

| Starting point | Go to |
| --- | --- |
| Coordinate, biomolecular, polymer, spatial-CSV, or surface data | {doc}`user_guide/input_adapters` |
| Existing embedded graph | {doc}`user_guide/workflow_overview` |
| Projection, PD code, or invariant question | {doc}`user_guide/projection_yamada` |
| Analytic knot or braid | {doc}`applications/analytic_knot_fields` |
| Nodal/material/phase-map workflow | {doc}`applications/index` |
| Import or native-backend problem | {doc}`troubleshooting` |

Each notebook includes a setup and runtime card. Use it to prepare the
dependencies and saved data for the example you want to explore.

## How the workflows connect

<div class="kg-wide-figure">
  <img src="site_figures/architecture.png" alt="KnottedGraph package architecture">
</div>

The modules share an embedded graph representation. Follow
{doc}`user_guide/workflow_overview` to trace the input, extraction, projection,
evaluation and visualization stages.

```{toctree}
:hidden:
:maxdepth: 2

Paper & data <paper_results>
Benchmarks <benchmarks>
Sanity checks <sanity_checks>
Cite our work <citing>
installation
quickstart
feature_status
user_guide/index
applications/index
troubleshooting
api/index
developer/index
```
