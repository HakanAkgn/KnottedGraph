# Completed computation and restored figures

This package reports executed numerical results. It is not the unavailable current 42-page manuscript.

## Grid outcomes

| Outcome | Hamiltonian | TPMS |
|---|---:|---:|
| evaluated | 7134 | 1168 |
| fixed_diagram_evaluated | 5 | 77 |
| time_budget | 0 | 7 |
| existing_cavity_route_not_revised | 1 | 71 |

## Analytic TPMS rectangles

| Family | Regular area (%) | Event-containing-cell area (%) | Unknown area (%) |
|---|---:|---:|---:|
| gyroid_to_diamond | 89.765625 | 10.187500 | 0.046875 |
| gyroid_to_schwarz_p | 93.203125 | 5.546875 | 1.250000 |
| schwarz_p_to_diamond | 93.921875 | 5.671875 | 0.406250 |

Event outcomes: {"certified_level_resolved_event": 249, "certified_multiseed_level_event": 7, "certified_stationary_event_in_cell": 1063, "exact_degenerate_gyroid_schwarz_event": 1, "exact_stationary_branch_event_in_cell": 50, "unknown": 109}.

These are regularity and event-existence certificates, not a complete classification into distinct ambient-isotopy classes.

## Closed-voxel manifold tests

| Input | Tested | PL manifold | Nonmanifold |
|---|---:|---:|---:|
| hamiltonian | 7140 | 4483 | 2657 |
| tpms | 1323 | 36 | 1287 |

The nonmanifold count refers to the represented closed voxel complexes, not to all analytic level sets.

## Explicit source-to-digital comparison witnesses

### digital_homology_conflict_examples

12 examples retained (a capped example set, not an exhaustive total).

Family `gyroid_to_diamond`, regular analytic rectangle `gyroid_to_diamond_lambda000_c000`.
First sample: lambda=0.0, c=0.0; digital Betti tuple [1, 17, 0].
Second sample: lambda=0.0, c=0.015; digital Betti tuple [1, 14, 0].
Both original parameter points lie inside the same freshly dual-replayed regular analytic rectangle.
Certificate: `comparison_witnesses/gyroid_to_diamond_lambda000_c000.json.gz`.

Family `gyroid_to_diamond`, regular analytic rectangle `gyroid_to_diamond_lambda000_c000`.
First sample: lambda=0.0, c=0.0; digital Betti tuple [1, 17, 0].
Second sample: lambda=0.05, c=0.015; digital Betti tuple [1, 23, 0].
Both original parameter points lie inside the same freshly dual-replayed regular analytic rectangle.
Certificate: `comparison_witnesses/gyroid_to_diamond_lambda000_c000.json.gz`.

Family `gyroid_to_diamond`, regular analytic rectangle `gyroid_to_diamond_lambda000_c000`.
First sample: lambda=0.05, c=0.0; digital Betti tuple [1, 17, 0].
Second sample: lambda=0.0, c=0.015; digital Betti tuple [1, 14, 0].
Both original parameter points lie inside the same freshly dual-replayed regular analytic rectangle.
Certificate: `comparison_witnesses/gyroid_to_diamond_lambda000_c000.json.gz`.

### spine_signature_variation_examples

12 examples retained (a capped example set, not an exhaustive total).

Family `gyroid_to_diamond`, regular analytic rectangle `gyroid_to_diamond_lambda000_c008`.
First sample: lambda=0.0, c=0.12; digital Betti tuple [1, 8, 0].
Second sample: lambda=0.05, c=0.135; digital Betti tuple [1, 8, 0].
Both original parameter points lie inside the same freshly dual-replayed regular analytic rectangle.
Certificate: `comparison_witnesses/gyroid_to_diamond_lambda000_c008.json.gz`.

Family `gyroid_to_diamond`, regular analytic rectangle `gyroid_to_diamond_lambda000_c008`.
First sample: lambda=0.0, c=0.135; digital Betti tuple [1, 8, 0].
Second sample: lambda=0.05, c=0.135; digital Betti tuple [1, 8, 0].
Both original parameter points lie inside the same freshly dual-replayed regular analytic rectangle.
Certificate: `comparison_witnesses/gyroid_to_diamond_lambda000_c008.json.gz`.

Family `gyroid_to_diamond`, regular analytic rectangle `gyroid_to_diamond_lambda000_c009`.
First sample: lambda=0.0, c=0.135; digital Betti tuple [1, 8, 0].
Second sample: lambda=0.05, c=0.135; digital Betti tuple [1, 8, 0].
Both original parameter points lie inside the same freshly dual-replayed regular analytic rectangle.
Certificate: `comparison_witnesses/gyroid_to_diamond_lambda000_c009.json.gz`.

## Remaining scope

Unavailable polynomial values, higher-valence fixed-diagram values and excluded cavity-route cases remain separate.
Unknown analytic rectangles are not relabeled as inequivalent. A selected-spine polynomial is not a complete source-solid invariant.
The original manuscript source, author metadata and all-word formula work were not replaced.
The existing cavity implementation and protected integration branch were not modified.

The figures retain the original grids and representative geometry, with no small-region relabeling or abstract polynomial fallback.
All figure arrays, exact coefficient legends, source hashes, local tests and example witnesses are provided in this package.
