# User Guide

<div class="kg-hero">
  <p class="kg-lead">Start with a small tested route, follow input handling, graph inspection, projection and invariant evaluation, then choose an application or publication workflow for your task.</p>
  <div class="kg-link-row">
    <a href="../feature_status.html">Choose A Supported Route</a>
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/00_user_guide.ipynb">Open the notebook map</a>
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/01_getting_started.ipynb">Getting Started</a>
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/02_core_workflows.ipynb">Core Workflows</a>
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/03_advanced_and_reproduction.ipynb">Advanced And Reproduction</a>
  </div>
</div>

The {doc}`../feature_status` page distinguishes public adapters, optional
features, application APIs, and external native backends before you choose a
notebook. New users can start with the {doc}`../quickstart`, then choose the
guide that matches their data.

## Recommended order

| Step | Read or run | What you will learn |
| --- | --- | --- |
| 1 | {doc}`../installation` and {doc}`../quickstart` | active version/environment and one expected exact result |
| 2 | {doc}`input_adapters` | source representation, selections and preparation steps |
| 3 | {doc}`workflow_overview` | graph contract, extraction/cleanup settings and provenance |
| 4 | {doc}`projection_yamada` | projection regularity, PD-code meaning, zero versus failure, and cost |
| 5 | {doc}`../applications/index` | one domain-specific application or advanced reproduction route |

The first three notebooks are progressive tutorials. Application notebooks
apply those concepts to scientific models; benchmark notebooks contain
correctness and performance checks with their saved results.

<div class="kg-wide-figure">
  <img src="../site_figures/knot_to_spatial_graph.png" alt="Knot and surface workflow leading to a spatial graph">
</div>

```{toctree}
:hidden:
:maxdepth: 2

workflow_overview
input_adapters
projection_yamada
repulsive_layout
../applications/index
```
