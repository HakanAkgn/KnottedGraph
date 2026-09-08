# New phase-map integration review

This update incorporates Hakan's `2b2ae6d` phase-map additions and the preceding
`091bfc0` timing-figure data into the user-integration branch. The scope is the
application/user layer. It does not change the input-format Main/S1/S2 figures,
Overleaf files, or protected extraction/projection/invariant implementations.

## Changes to review

1. Moved the reusable material and compact-TPMS engines into the installed
   `knotted_graph.applications.phase_map_examples` package. Kept the old script
   names as compatibility entries and left specialized paper layouts as
   research scripts.
2. Removed developer-machine import paths and temporary-directory assumptions.
   Reproduction tools now accept explicit data and output paths. Geometry
   lookup follows the selected dataset, including the archived `data/` and
   sibling `geometry/` layout. Non-default TPMS domain parameters are recorded
   and checked before representative surfaces are regenerated.
3. Added independent `inspect`, `plot` and `scan` commands. Inspection/replotting
   works in a base installation; scans require the optional stack. A quick
   profile and no-computation `--dry-run` precede dense research workloads.
   Guided scans refuse to overwrite a non-empty result directory.
4. Made scientific status explicit. Raw signatures, exact computations,
   adaptive fills, large-core fallbacks, errors, display merges and contraction
   groups are not interchangeable. New guided scans disable adaptive fills,
   island smoothing and manual signature merges. Historical presentation
   processing remains available only through explicit research options.
5. Added raw categorical PNG/PDF output with Computer Modern mathematics and a
   JSON key containing category IDs, colors, signatures, polynomials, coordinate
   axes and the grid. Missing/error/adaptive-fill cells are distinguished.
6. Connected a guided walkthrough to the root README, User Guide index,
   applications navigation, feature-status matrix and API reference.
7. Added two website viewers without replacing the accepted nodal-only demo.
   Regions load on selection; all original payload values and geometry are
   preserved. Sphinx creates the split assets from pinned reference HTML, so
   another large generated copy is not committed. Mobile sizing and loading
   failures have explicit handling.
8. Added application/portability/regression tests and included the base-safe
   cases and new research scripts in the existing CI checks.

## Acceptance scope

Final local verification on 2026-09-08 used Linux / Python 3.13.7: **376 tests
passed**, with **69.52%** package coverage (the existing 61% gate is unchanged).
Ruff, the configured Pyright checks, notebook normalization/portability,
repository consistency and the protected-code audit passed. Sphinx built with
warnings treated as errors, and generated links/assets passed validation.
A clean base wheel passed Quick Start and the base-safe phase-map tests;
optional-stack-only cases were explicitly skipped in that base environment.
The source comparison confirmed 93 moved scientific functions/classes were
unchanged; the reviewed differences concern orchestration, parameters, paths,
presentation metadata and explicit processing defaults.

Chromium checks covered seven transitions / two modes in the material/nodal
viewer and three transitions / four modes in the TPMS viewer. Each initially
requested only one region attachment. Switching, clicks, narrow-screen layout
and the blocked-CDN message passed without page errors or unexpected failed
requests. The original dense TPMS audit also passed for all 1,323 records.

Validation is performed on compute nodes, not the cluster login node. It covers:

- base and full-optional-stack tests, native wheel construction and clean-wheel
  installation/Quick Start;
- real coarse material and TPMS scans with record round trips, rather than only
  mocked command success;
- saved TPMS grid, geometry, compactness, endpoint and region audits;
- material and TPMS representative-viewer generation on coarse records, and
  the saved-data publication figure renderers;
- Sphinx HTML with warnings as errors and built-site link/asset checks;
- Chromium switching/clicking across material and TPMS modes, one-region initial
  loading, mobile layout, and blocked-CDN error handling;
- source-level comparison of moved scientific functions/classes, unchanged
  reference assets, and lossless reconstruction of all 108 material/nodal plus
  444 TPMS region payloads.

The original dense scans are retained, **not rerun in full**. Passing these
checks is not a proof of resolution convergence, all-parameter topological
equivalence, or cross-platform scientific reproducibility. The existing GitHub
workflow matrix still provides separate platform/version coverage. Local
validation is not the same as a deployed website or a green remote CI run.

## Scientific feedback requested

- Are the historical manual material signature merges and selected-component
  geometry conventions described accurately? They are preserved but no longer
  silently enabled in the introductory route.
- Are the distinctions between exact Yamada, structural signatures, adaptive
  classification and contraction-based display groups clear enough?
- Should the saved material viewer keep its five earlier nodal transitions
  alongside TiB2/Co2MnGa? The current version preserves them and opens TiB2 first.
- Which additional convergence or full-resolution reproduction checks should
  accompany publication? Those require a separate compute plan and review.
