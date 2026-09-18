# Corrected dense TPMS scan: scientific scope and revision evidence

This audit concerns the corrected 1,323-cell scan at 64³ voxels, 21 λ values and 21 c values, exact graph-edge limit 18, and **no display filtering**. Records SHA-256: `b4a8dc6c4feb8905f2c17e8b0567f5a9112135e537f85f764e5af45c0bf09eb2`.

| Family | Subcubic polynomial cells / values | Diagram polynomial cells / values | Abstract graph summaries | Unsupported cavity cells | Non-error raw signatures |
|---|---:|---:|---:|---:|---:|
| Schwarz-P → Diamond | 190 / 1 | 0 / 0 | 179 | 72 | 14 |
| Gyroid → Schwarz-P | 250 / 2 | 7 / 3 | 183 | 1 | 27 |
| Gyroid → Diamond | 122 / 3 | 44 / 19 | 274 | 1 | 77 |

Totals: 562 subcubic evaluations, 51 higher-valence fixed-diagram evaluations, 636 abstract graph summaries, 74 unsupported cells. The 1,249 evaluated graph cells have **zero β0 or β1 mismatches**. All 74 numerical classification errors are explicitly cavity-bearing masks; there are no other classification errors. These masks contain 145 enclosed voids in total. Graph fields are empty in the failed records because the cavity check prevents reconstruction; their apparent β0 mismatches and the two apparent β1 mismatches must not be counted as failed *completed* reconstructions. All 1,323 surface meshes are closed and free of reported nonmanifold/open edges; no sampled volume touches the outer box. Closed surface alone is not a handlebody certificate.

The raw signatures mix different evidence scopes. Their numbers are **not counts of topology or embedding-equivalence classes**. Equal normalized Yamada polynomials are not a complete spatial-graph classifier; higher-valence diagram values have a different invariance scope; graph summaries do not encode embedding; and Yamada on one chosen spine is not automatically invariant under handlebody spine changes. Experimental spatial contraction is separate from this scan.

## Filtering evidence

The archived original map used abstract contraction groups followed by a sequential 4-neighbor component filter, minimum component size 4 and at most 3 passes. Reproducing it exactly from the archived HTML's raw and displayed grids changes **171 distinct cells**: 18 P→D, 67 G→P, 86 G→D. Operation counters 18/71/94 include repeated edits. `tpms_archived_filter_changed_cells.json` retains every raw/final label and changed-cell mask, with source hash; `audit_archived_tpms_filter.py` can replay it from this compact extracted JSON without the large HTML.

Applying that frozen legacy algorithm to the **corrected** raw grids is a separate sensitivity diagnostic:

| Family | Non-error labels at min size 1 / 2 / 4 / 8 | Unique final changed cells at size 4 | Different cells after reversing arbitrary label IDs at size 4 |
|---|---:|---:|---:|
| P→D | 14 / 2 / 2 / 2 | 12 | 7 |
| G→P | 27 / 15 / 11 / 7 | 63 | 51 |
| G→D | 77 / 46 / 20 / 10 | 132 | 82 |

The two isolated cavity errors in the Gyroid-based scans receive non-error display labels under the old filter. This is evidence for retaining unresolved cells and avoiding physical interpretation of neighborhood relabeling. The corrected main figure is unfiltered. The historical 3/11/14 displayed classes and five-regime cut are removed, rather than defended as numerical-noise removal. Root separately improved the optional filter; these frozen historical diagnostics must not be relabeled as measurements of that new rule.

## Resolution and sampling

The 135-cell staged audit consists of three voxel dimensions (24/64/96), three families, five λ values (0/.25/.5/.75/1), and three c values (0/.105/.3). Among 45 matched family/parameter entries, **13 change mask Betti tuple** across the resolutions: 12 between 24 and 64; 8 between 64 and 96. Shared endpoints occur in multiple family entries, so these are not 45 independent shapes. For the Gyroid endpoint at c=.105, β1 is 5/14/11. At c=0, its tuple is (1,5,0)/(4,5,0)/(1,11,0). These changes precede extraction. They show that requiring Betti agreement cannot itself establish convergence of the sampled volume.

In the constrained rerun, all 129 non-cavity samples preserve β0 and β1 through extraction and postprocessing; six β2-positive samples remain explicitly unsupported. Zero-radius extraction matched β0 and β1 in all 135 samples in the original staged probe, but that empirical result is not a general graph-spine guarantee.

Nested parameter subsets measure sampling sensitivity, not independently recomputed convergence:

| Family | Raw non-error signatures at 21 / 11 / 6 points per axis | Corresponding 4-neighbor raw signature regions |
|---|---:|---:|
| P→D | 14 / 5 / 2 | 14 / 5 / 2 |
| G→P | 27 / 18 / 12 | 64 / 25 / 13 |
| G→D | 77 / 41 / 19 | 142 / 59 / 21 |

## Minimal coherent figure and text revision

`PorousMaterialPhaseMap_audited.pdf` is a six-panel replacement for the existing main figure, with the same family order. Top panels show unfiltered voxel β1, numerals for β0>1, and hatching for β2>0. Bottom panels distinguish normalized subcubic polynomial evaluations, high-valence diagram values, abstract graph summaries, and unsupported cavities. The source JSON keeps every plotted value. No appendix plot is needed.

The staged TPMS text replaces topology-class counts with operational scope and explicit unresolved counts, supplies filtering and sampling details in prose, and makes the Morse/physical interpretation prospective. The critical-point and neck-area equations remain. Analytic critical branches, continuum convergence, preserved spatial reconstruction, and mechanical/fluid/thermal calculations have not been established by this scan; measured response changes and guaranteed topology-preserving corridors are therefore not claimed.

## Separate experimental contraction result

The conservative exact rational PL certificate accepts 20 eligible local moves in 12 archived rounded graphs. Nine selected neighboring pairs sharing old abstract-contraction labels give eight unknowns and one already-identical singleton export. The singleton pair is in a cavity-bearing region; it cannot certify source-volume equivalence. No nontrivial cross-parameter equivalence was demonstrated, and no merge is applied to the corrected maps. The implementation refuses unsupported geometry or exhausted budgets instead of treating lack of a witness as inequivalence. Its 14 adversarial tests include equal abstract graphs with differently knotted/linked embeddings.
