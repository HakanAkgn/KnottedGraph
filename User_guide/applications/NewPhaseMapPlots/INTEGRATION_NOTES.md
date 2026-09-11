# New phase-map integration review

## Boundary-resolved update — 2026-09-11

Hakan's `643fef8` supplies the four helpers required by `56bfbab`:
`boundary_filling_groups`, `enclosed_void_masks`, `resolve_volume_mask` and
`volume_topology`. The exact upstream snapshot imported successfully and passed
all eight core phase-map tests on a PBS compute node. This integration adopts
that implementation and the `56bfbab` material records/figures/HTML, retaining
the installed-package engines and portable source-script entry points.

The material reference contains 46,020 records: 14,634 directly classified,
31,184 adaptive fills and 202 resolution calibrations. TiB2 has 35 raw
signatures; Co2MnGa has 53. The material viewer contains 88 region entries
across seven transitions, including its five earlier nodal transitions. TPMS
retains its 1,323 records and 444 region entries. Region counts are not counts
of distinct polynomials. CSV/JSON signatures, polynomial strings, source/status
metadata and all split-viewer payloads have been checked against the accepted
assets; no dense reference scan was rerun.

The reader and raw plot now distinguish calibrated records from adaptive
fills, including an anchor energy of zero. Guided scans keep display merges,
adaptive fill and resolution calibration off. Research CLI flags enable C6
display grouping and resolution calibration separately. In the relocated
energy sampler, two unconditional calibration calls have moved to the
explicit main/reuse/extend handling; the sampling logic and the teacher's
calibration algorithm remain unchanged. A targeted audit compares this exact
orchestration difference rather than broadly exempting the function.

Representative material geometry now follows the teacher's dominant resolved
body and outer/nested boundary fillings. Regenerating selected material
transitions preserves other transitions in the supplied viewer. The old
component-fraction option produces an explicit migration error.

Folder READMEs explain the route from saved records to a small scan. A base-only
CSV-to-Yamada example completes input loading, embedding validation, projection
and exact evaluation using existing APIs. Notebook export saves execution
copies, HTML and status files outside the source tree. Mathematics and analytic
fields are removed from automatic full execution in notebook CI because their
saved defaults include large galleries or paper-mode scans; source checks and
the full scientific notebook contents are preserved.

Validation/evidence for this update is recorded in the continuation handoff.
The first complete gate passed 384 tests (70.67% coverage), style/types,
portability, repository consistency, Sphinx with warnings as errors, generated
links, reference audits, coarse representative geometry, saved-panel renderers,
clean base-wheel installation/examples and Chromium checks. A subsequent
calibration-routing correction passed the complete suite of 386 tests (71.16%
coverage), the exact source audit and an updated clean base wheel (29 base-safe
phase-map tests passed, seven optional-stack cases skipped). Browser checks cover transition/mode switching,
one-region initial loading, mobile layout, blocked-CDN messaging, failed region
fetch recovery and the unchanged standalone nodal viewer.

Notebook coverage is deliberately explicit: introductory notebooks, offline
protein and Yamada sanity checks run with saved defaults; mathematics covers
original cells 0–16; analytic fields runs all cells with its existing fast
branch selected in a scratch copy. The excluded scan sections, research
galleries and native Repulsor solver are not claimed as executed. Visual review
samples rendered outputs and separately checks the source-page layouts.

This is a local integration review, not a website/PyPI release. The independent
Input Adapters Guide, publication experiments and deployment retain their
separate confirmation boundaries. Main/S1/S2, paper sources, protected core
implementations, benchmark scientific cells and tracked scientific assets are
preserved.

## Earlier integration — 2026-09-08

The earlier update incorporated Hakan's `2b2ae6d` phase-map additions and the preceding
`091bfc0` timing-figure data into the user-integration branch. The scope is the
application/user layer. It does not change the input-format Main/S1/S2 figures,
Overleaf files, or protected extraction/projection/invariant implementations.

### Earlier changes

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

### Earlier acceptance scope

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

- Are the current C6/manual display groups, calibration provenance and
  boundary-resolved geometry described accurately? These remain distinct from
  raw directly classified results in the introductory route.
- Are the distinctions between exact Yamada, structural signatures, adaptive
  classification and contraction-based display groups clear enough?
- Should the saved material viewer keep its five earlier nodal transitions
  alongside TiB2/Co2MnGa? The current version preserves them and opens TiB2 first.
- Which additional convergence or full-resolution reproduction checks should
  accompany publication? Those require a separate compute plan and review.
