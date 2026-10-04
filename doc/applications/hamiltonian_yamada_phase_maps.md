# Hamiltonian Yamada Phase Maps

This workflow accepts an in-memory Hamiltonian or Bloch-vector model and
records extracted topology at each point of a two-parameter sampling grid.

For the additional material-parameter and compact-scaffold plots, continue to
{doc}`material_phase_maps`. That guide starts from saved records and
provides a small compute example. The interactive nodal examples are shown
below.

<div class="kg-hero">
  <p class="kg-lead">Use the interactive result to select a transition and a phase region, then inspect a representative exceptional surface and its simplified spatial-graph skeleton. Use the notebook when you need to regenerate the grid, caches, audits, or figures.</p>
  <div class="kg-link-row">
    <a href="../demos/hamiltonian_yamada_phase_map.html">Open the interactive result full screen</a>
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/06_hamiltonian_yamada_phase_maps.ipynb">Open the reproduction notebook</a>
    <a href="../api/applications.html">Phase-map API</a>
  </div>
</div>

## Interactive result

<iframe
  class="kg-interactive-demo"
  src="../demos/hamiltonian_yamada_phase_map.html"
  title="Interactive Hamiltonian Yamada phase map and representative geometry"
  loading="lazy"
  sandbox="allow-scripts allow-same-origin"
></iframe>

<p class="kg-caption">Choose a classification mode and transition above the phase map. Clicking a stable region updates the representative 3-D exceptional surface, skeleton, graph statistics, and polynomial. The embedded artifact loads Plotly from a pinned CDN URL; use the full-screen link if an iframe is blocked.</p>

## Data flow

For a sampled cell $(\lambda,\Gamma)$, the notebook follows

$$
H(\mathbf{k};\lambda,\Gamma)
\longrightarrow \text{filled exceptional region}
\longrightarrow \text{skeleton}
\longrightarrow G\subset\mathbb{R}^3
\longrightarrow \Upsilon(G;Y).
$$

The reusable `make_yamada_phase_map(...)` API stores one record per cell,
including graph size, connected components, cycle rank, polynomial or error,
and a phase signature. The notebook adds row caches, connected-region
stabilization, endpoint checks, classic phase labels, and a second view that
groups regions up to selected contraction moves.

<a id="interpretation-boundaries"></a>

## Read the phase-map records

- Each colored cell records a computation at the stated parameter values and
  grid resolution; boundaries follow neighboring sampled cells.
- An extraction failure has an error record. An exact zero has a successfully
  evaluated polynomial value.
- The displayed surface gives physical/geometric context. The black skeleton
  is the graph used for topological analysis.
- The classic view groups Yamada signatures. “Up to contraction moves” groups
  representative graphs under the stated contraction convention.
- Inspect resolution, sampling-window contact, endpoint behavior and
  representative geometry when comparing regions.

The full notebook uses up to 60 lambda samples and 50 candidate Gamma samples
with a $120^3$ volume per evaluated cell. Start from the saved interactive
result, then use the notebook's configuration and cache settings for regeneration.
