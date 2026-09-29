# Application And Reproduction Workflows

<div class="kg-hero">
  <p class="kg-lead">Explore examples built around mathematical graphs, molecular structures, knot fields and material models. Choose a route that matches your data; each guide lists the setup, expected outputs and links to its notebook.</p>
  <div class="kg-link-row">
    <a href="../feature_status.html">Check Feature Status</a>
    <a href="analytic_knot_fields.html">Analytic Knot Fields</a>
    <a href="yamada_formula_discovery.html">Formula Discovery</a>
    <a href="hamiltonian_yamada_phase_maps.html">Hamiltonian Phase Maps</a>
  </div>
</div>

## Choose the right level

| Route | Intended reader | Start here when... | Execution class |
| --- | --- | --- | --- |
| {doc}`mathematical_investigations` | first-time application user | you have a named abstract graph family | guided; base API |
| {doc}`analytic_knot_fields` | application user | you have a knot/link name, torus type, or braid | guided, then compute-intensive extraction |
| {doc}`material_fingerprints` | domain user | you already have an in-memory nodal/material model | optional `nodal`; grid work |
| {doc}`protein_derived_spatial_graphs` | input user | you want to load a PDB/mmCIF backbone and prepare a spatial graph | guided backbone input; choose domain connections for interaction networks |
| {doc}`yamada_formula_discovery` | researcher reproducing a result | you need the exact dataset/held-out symbolic checks | advanced publication reproduction |
| {doc}`hamiltonian_yamada_phase_maps` | domain researcher | you need a two-parameter Hamiltonian topology scan | advanced, cached, compute-intensive |
| {doc}`material_phase_maps` | reader of the new material plots | you want to inspect saved records, then try a coarse scan | base for saved records; optional and compute-intensive for new scans |

For a first run, try the {doc}`../quickstart`. To inspect the published results,
start with {doc}`../paper_results`, then open the corresponding application
guide for its saved files, setup and calculation steps.

<div class="kg-wide-figure">
  <img src="../site_figures/feature_image1.png" alt="Material Fermi-surface fingerprints and Yamada polynomials">
</div>

```{toctree}
:hidden:
:maxdepth: 1

material_fingerprints
protein_derived_spatial_graphs
mathematical_investigations
analytic_knot_fields
yamada_formula_discovery
hamiltonian_yamada_phase_maps
material_phase_maps
```
