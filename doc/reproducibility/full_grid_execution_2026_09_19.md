# Completed original-grid reruns — 19 September 2026

This report concerns the requested `research/solid-spine-maps-41dc87c` branch.
The all-word formulas, existing cavity-processing implementation and original
manuscript/artwork were not revised by these numerical runs.

## Full original Hamiltonian grid

Run: https://github.com/HakanAkgn/KnottedGraph/actions/runs/35411673566

Numerical source commit: `9be4e508f2ca7a48f580253be5de5d5392f055d0`.
The collector job `105819282617` completed successfully. All 20 numerical
shards completed. The collector checked exact, disjoint, same-commit coverage
of the requested grids and the archived collapse/geometry hashes for completed
evaluations. No small-region display filter or abstract polynomial fallback
was used.

All **7,140 original Hamiltonian parameter points** were executed at **120 cubed**.
The exact arrays are `np.linspace(0, 1, 60)` and the retained prefixes of
`np.linspace(0.30, 5.25, 50)`. Thus the interpolation spacing is 1/59, including
both endpoints. The five energy-prefix lengths remain 16, 24, 16, 24 and 39.
The original Bloch functions and finite box are retained; opposite faces are
not periodically identified. Sample-mask and reconstruction scope remain
separate from analytic-source topology.

| Family | Requested points | Spatial/isolated-vertex evaluations | Fixed-diagram evaluations | Time budget | Existing cavity route not revised |
|---|---:|---:|---:|---:|---:|
| Hopf to trefoil | 960 | 960 | 0 | 0 | 0 |
| Hopf to Solomon | 1,440 | 1,429 | 1 | 10 | 0 |
| Unknot to trefoil | 960 | 960 | 0 | 0 | 0 |
| Unknot to Solomon | 1,440 | 1,435 | 2 | 2 | 1 |
| Trefoil to cinquefoil | 2,340 | 2,332 | 2 | 6 | 0 |
| **Total** | **7,140** | **7,116** | **5** | **18** | **1** |

The 7,116 successful entries include the separately audited isolated-vertex
case type. Nontrivial subcubic entries require two distinct generic projection
evaluations with exact equality after normalization. These are two geometric
views using the existing evaluator, not two independently implemented
polynomial algorithms. Higher-valence results remain fixed-diagram values.
Timeouts retain no successful polynomial label. The time budget is 90 seconds
per cell, including its sampling/reconstruction/evaluation work.

This completes execution of the full original grid, not every polynomial
calculation and not a proof of all source-solid phase classes.

## Full TPMS parameter grid at 120 cubed

The same workflow executed all **1,323** points of the original three
21-by-21 parameter grids at an explicitly recorded spatial resolution of
120 cubed. This is a higher-resolution rerun of the parameter grid, not a
claim that the old 64-cubed source masks were identical. The original fields,
spherical domain, interpolation values and level values are retained.

| Family | Requested points | Spatial/isolated-vertex evaluations | Fixed-diagram evaluations | Time budget | Existing cavity route not revised |
|---|---:|---:|---:|---:|---:|
| Schwarz-P to Diamond | 441 | 261 | 39 | 70 | 71 |
| Gyroid to Schwarz-P | 441 | 381 | 0 | 60 | 0 |
| Gyroid to Diamond | 441 | 359 | 34 | 48 | 0 |
| **Total** | **1,323** | **1,001** | **73** | **178** | **71** |

The cavity entries were recognized and left to the excluded existing route;
this run does not claim that those values were recomputed through that route.
Its code was not changed. No failure or timeout was converted into a phase.

## Reconstruction and arithmetic validation

The faster standalone C++ runner implements the existing elementary cubical
free-face collapses using dense cell storage. It separately rebuilds the
input complex and replays incidence, free-face and maximality conditions for
every recorded move. It then compares the complete terminal complex.

Thirty randomized and six structured differential tests match the complete
move sequence and endpoint against the existing Python implementation; the
Python verifier independently checks those test certificates. Full-resolution
map cells use the native replay and an additional independent digital
component/cycle check. This report does not claim that every full-resolution
certificate was additionally expanded and replayed by the Python verifier.

Original graph polylines are retained. The witness proves a deformation
retraction of the specified closed voxel-cell complex when a graph endpoint
is reached. It does not, by itself, prove analytic-to-voxel correspondence,
manifoldness, a regular-neighborhood correspondence, or a complete invariant.
All such stronger claims remain false in the output metadata.

The native-runner full-suite validation passed 618 tests. The exact-plan and
six-case 120-cubed preflight passed 619 tests. The later critical-event module
validation passed 645 full-repository tests, as well as its scoped lint and
repository-consistency check. These runs are separate from proof of a
continuous-source correspondence.

## Retained data

The complete collector artifact is `complete-full-map-summary`, ID
`10575332109`, containing both full record arrays and `full_map_manifest.json`.
Its ZIP size is 7,816,202 bytes and SHA-256 is:

`f3f0c8a8bb866a2a92f49192df162126351b218afc38b63c86dddc0ebb966e44`.

The 20 `full-map-shard-N` artifacts from the same run retain sampled masks,
physical grid axes, compressed full collapse sequences, terminal geometry,
per-cell records and provenance. Their download hashes were verified by the
collector. Artifact retention is 30 days, not a permanent research-data
archive. The collected records include artifact-relative paths and per-cell
hashes. A copy of the collector ZIP was exported to the conversation as
`KnottedGraph_complete_7140_Hamiltonian_1323_TPMS_records.zip`.

Reproduction from the numerical source commit uses the frozen lockfile and
`dev/run_full_resolution_maps.py`; the workflow records all 20 disjoint shard
commands. The native kernel is built with:

```sh
g++ -O3 -std=c++17 dev/full_map_collapse.cpp -o /tmp/kg-collapse
uv run --no-sync python dev/run_full_resolution_maps.py \
  --mode hamiltonian --dimension 120 --kernel /tmp/kg-collapse \
  --out /tmp/kg-hamiltonian-shard-0 --shard 0 --shards 20 --case-seconds 90
```

Use new output directories and execute each shard index exactly once. The
collector rejects incomplete or mixed grids.

## TPMS continuation and final manuscript scope

The adaptive continuation calculation is distinct from these polynomial maps:
it certifies regular analytic source-solid parameter rectangles using both
interval backends. The local critical-event method additionally checks whether
unresolved cells contain actual stationary points. The completed first-shard
probe examined all 28 unresolved leaves: 24 contained certified stationary
points and four remained unknown. Its accepted points were rechecked by both
interval backends with the full critical-value enclosure inside the queried
level interval. This counts cells containing certified events, not 24 distinct
physical transitions or 24 distinct critical curves.

No new final manuscript or replacement main figure is part of this numerical
release. The current 42-page source archive is not in the retained code branch.
A proposed figure-source write was rejected by the tool and was not applied;
it was not retried by an alternative route. The subsequent proposed
same-geometry Hamiltonian-timeout retry commit was also rejected and was not
applied or executed. These proposed edits are excluded from the results above.
