# Workflow Overview

Start at the stage that matches your data. Curves, graph tables, surfaces and
scientific models meet at a shared embedded spatial graph, which can then be
projected, visualized and analysed.

<div class="kg-link-row">
  <a href="../feature_status.html">Choose by starting object</a>
  <a href="input_adapters.html">Load external data</a>
  <a href="projection_yamada.html">Projection and Yamada</a>
  <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/02_core_workflows.ipynb">Open Core Workflows notebook</a>
</div>

<div class="kg-wide-figure">
  <img src="../site_figures/skeletonization_steps.png" alt="Surface and volume skeletonization stages leading to an embedded spatial graph">
</div>

## The reusable path

```text
external data / analytic field / in-memory model
    -> adapter or application-specific extraction
    -> embedded MultiGraph(pos, pts)
    -> validate and simplify without changing intended topology
    -> choose a regular planar projection
    -> inspect crossings and PD-code records
    -> evaluate Upsilon(G; Y)
    -> retain settings, diagnostics, and provenance
```

An abstract graph can go directly to the crossing-free graph evaluator. For
a loaded `PolyData` surface, choose an extraction method to obtain a graph.
Field and Hamiltonian workflows may start with 3-D sampling and skeletonization.

<a id="stage-1-identify-what-you-actually-have"></a>

## Stage 1: choose your starting object

Start by distinguishing these cases:

- **Ordered coordinates** describe one sampled curve. For a branching
  network, supply node and edge connections.
- **Node and edge data** can represent a spatial MultiGraph directly.
- **Surface or volume data** require an extraction decision; multiple skeletons
  can be plausible for the same physical geometry.
- **Analytic knot fields and Hamiltonians** are application objects that must be
  sampled at an explicit domain, resolution, and level/energy.
- **Abstract graphs** contain topology but no embedding or crossing data.

Use {doc}`input_adapters` for direct public input support and
{doc}`../feature_status` for application routes.

## Stage 2: inspect the graph contract

```python
from knotted_graph.core import ensure_embedding

graph = ensure_embedding(graph)
print("nodes:", graph.number_of_nodes())
print("edges:", graph.number_of_edges())
print("degrees:", sorted(dict(graph.degree()).values()))
```

Before simplifying anything, inspect:

- connected components;
- vertex degrees and parallel edges;
- leaf/bridge structure;
- edge sampling density;
- coordinate units and scale; and
- whether boundaries or periodic faces were touched during extraction.

Smoothing, short-edge contraction and leaf removal act on different parts of
the graph. Record the operation and compare graph diagnostics before and after
it to check the intended topology.

## Stage 3: select a regular projection

```python
from knotted_graph.projection import select_projection

projection = select_projection(
    graph,
    num_rotation_samples=16,
)
print("angles:", projection.rotation_angles)
print("crossings:", projection.num_crossings)
```

A valid projection must avoid degeneracies such as a projected vertex lying on
an unrelated edge or multiple events collapsing to the same point. The selector
samples candidate views and returns the chosen rotation and diagram
diagnostics for inspection.

Projected crossings are diagram events, separate from graph vertices. Read the
crossing and arc records before interpreting a PD code.

## Stage 4: compute the invariant with provenance

```python
import sympy as sp
from knotted_graph.projection import compute_yamada_polynomial

Y = sp.Symbol("Y")
result = compute_yamada_polynomial(
    graph,
    Y,
    n_jobs=1,
    return_result=True,
)

print("Upsilon(G; Y) =", result.polynomial)
print("projection crossings =", result.projection.num_crossings)
```

For a reproducible result, retain the selected projection, rotation policy, normalization convention,
worker count, backend status, and any failed candidate views.

State evaluation grows approximately as \(3^c\) with projected crossing count
\(c\). Use one worker by default, inspect the chosen view, and request cluster
resources explicitly for expensive cases.

## Extraction and skeletonization are broader than invariant calculation

Skeletonization can be used to obtain compact centerline graphs for geometry
inspection, comparison, routing, or visualization even when no Yamada
polynomial is needed. The extracted representation depends on the field,
sampling and extraction method, so record these choices when comparing graphs.

For sampled surfaces/volumes, report at least:

- sampling domain and grid resolution;
- threshold, level, energy, or gap parameter;
- boundary contact and periodicity policy;
- skeletonization/extraction method;
- cleanup parameters;
- graph diagnostics before projection; and
- resolution or perturbation checks used to establish stability.

## Tutorial versus reproduction notebooks

The notebooks have different purposes:

| Level | Use | Expectation |
| --- | --- | --- |
| `01_getting_started` | first successful graph/projection/invariant run | modest runtime, copyable cells, expected outputs |
| `02_core_workflows` | inspect each processing stage | reduced-resolution tutorial defaults |
| `03_advanced_and_reproduction` | robustness and diagnostic policy | more computation and domain knowledge |
| application notebooks | domain-specific scientific workflows | optional extras; may include paper-sized modes |
| benchmark notebooks | correctness/performance evidence | native backends, caches, and compute resources may be required |

Start with the Quick Start to explore a small 3D example, then follow the
notebook level that matches your task. Reproduction notebooks list the data
and dependencies used for the corresponding paper calculation.

## Reproducibility checklist

For a result intended for comparison or publication, save:

1. package version and Git commit;
2. input identifier, file checksum, or analytic constructor parameters;
3. coordinate units, closure, chain/model, and domain choices;
4. extraction/skeletonization and cleanup settings;
5. graph fingerprint and node/edge/degree diagnostics;
6. projection angles, crossing count, and failed-view diagnostics;
7. Yamada method, normalization, variable, worker count, and backend status;
8. warnings/issues returned by adapters or applications; and
9. runtime environment and relevant optional/native dependency versions.

Continue with {doc}`projection_yamada` for diagram/invariant details or open the
[Core Workflows notebook](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/02_core_workflows.ipynb)
for a staged executable example.
