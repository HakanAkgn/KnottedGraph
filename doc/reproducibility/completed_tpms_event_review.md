# Completed TPMS event review and full-grid execution

## Tested revision and execution

The tested code is `2e85ec55840b60735830027f190219cf7d495115` on
`research/solid-spine-maps-41dc87c`. This report is a documentation-only
follow-up. GitHub Actions run `35415827400`, job `105824182831`, completed
successfully on 19 September 2026. The actual job log records:

- **650 pytest tests passed**, with no tests removed.
- Scoped Ruff checks passed.
- The repository consistency audit passed.
- Every retained positive local critical-point certificate passed fresh
  primary and rational-trigonometric replay.
- The exact stationary-branch identities and rational-domain inclusion were
  reverified before supplying any additional event witness.
- All 1,479 unresolved adaptive leaves were accounted for exactly once.
- Tracked source files remained unchanged during execution.

The existing all-word formulas, cavity-processing implementation, original
manuscript and original artwork were not revised. The protected integration
branch remains separate; none of these results was merged into it.

## Recovery of the completed event search

The original event search, run `35414418156` at
`24bc40c4b75c5d20ac1af5a685f061c6682d59b9`, completed all 1,479 selected cells.
It obtained 1,063 certified stationary-event cells and 416 unknown cells.
Its final aggregation failed because the final shard had zero selected
cells and the driver did not serialize an empty `event_records.json` file.
This was an output-contract bug, not a failed mathematical certificate.

The driver now writes an empty result array before entering its processing
loop. Regression tests check that zero-row cohorts produce a real empty file,
and that a missing file in a nonempty cohort is never accepted as empty.

The recovery did not rerun root searches, modify retained artifacts or invent
missing positive results. It re-audited the original adaptive partition and
matched every shard's exact selected IDs, source-record hash and source plan.
The absent legacy file was accepted only for shard 19, whose plan and completed
summary both explicitly recorded an empty selection. Every positive local
certificate was then freshly checked using both interval backends against the
caller's reconstructed original analytic field and the queried parameter cell.

## Final event outcomes

The known exact Schwarz-P to Diamond branches supply positive event witnesses
for 50 additional cells that the isolated-root search had left unknown. The
intersection calculation uses exact rational arithmetic on the stored binary
parameter endpoints. It identifies an actual parameter pair and stationary
point inside the requested closed cell. It does not infer membership from a
rounded plotted line or a tolerance.

| Family | Local critical-point certificates | Additional exact-branch witnesses | Still unknown | Total initially unresolved leaves |
|---|---:|---:|---:|---:|
| Schwarz-P to Diamond | 283 | 50 | 56 | 389 |
| Gyroid to Schwarz-P | 265 | 0 | 170 | 435 |
| Gyroid to Diamond | 515 | 0 | 140 | 655 |
| **Total** | **1,063** | **50** | **366** | **1,479** |

Thus **1,113** formerly unresolved leaf rectangles are now proved to contain
at least one stationary event. This is not a count of 1,113 distinct critical
points, distinct curves, or topology-changing transitions. A stationary curve
can cross several parameter rectangles, and adjacent closed rectangles can
share an event on their boundary.

The exact branches are verified from the original scalar field:

$$c=1-2\lambda\quad\text{at }(\pi,0,0),\qquad
c=2\lambda-1\quad\text{at }(\pi,\pi,0),$$

with their coordinate/sign images. At $\lambda=1/2,c=0$, the entire line
$(\pi,t,0)$ inside the ball is stationary. A verifier requiring an isolated
unique root cannot certify a box containing a segment of that line. The exact
identities address that limitation without increasing a search budget or
changing the field. Details and Hessians are in
`exact_tpms_stationary_branches.md`.

## Area accounting of the original parameter domains

The original domains are $0\le\lambda\le1$ and $0\le c\le0.3$, with their
stored binary endpoints retained. The adaptive solver started from the
original 1,200 adjacent parameter rectangles and refined unresolved cells
twice. Every successful regularity leaf has the original two full interval
replays. The complete partition contains 2,538 accepted regularity leaves
and 1,479 unresolved leaves before event analysis.

The following are percentages of each family's parameter-domain area:

| Family | Certified regular rectangles | Rectangles proved to contain an event | Rectangles still unknown |
|---|---:|---:|---:|
| Schwarz-P to Diamond | 93.921875% | 5.203125% | 0.875000% |
| Gyroid to Schwarz-P | 93.203125% | 4.140625% | 2.656250% |
| Gyroid to Diamond | 89.765625% | 8.046875% | 2.187500% |

The three columns partition each domain exactly in the recorded dyadic cell
accounting. The event column is the area of finite rectangles containing
proved events; it is **not** the area of the critical set itself and does not
assert that every point of those rectangles is singular. Boundaries shared
by closed rectangles have zero area in this accounting.

The combined unknown area across the three equal-area parameter domains is
1.90625%. This is **not** a claim that 98.09375% of the domain has been assigned
to distinct ambient-isotopy classes. Regular blocks provide positive
continuation equivalence; event blocks provide existence of an obstruction to
certifying the entire block as regular. Neither result proves inequivalence
between disconnected regular-cover components or completes the global phase
classification.

## Full original-resolution numerical grids

The separate full-grid run `35411673566` at
`9be4e508f2ca7a48f580253be5de5d5392f055d0` completed every requested point.
The exact 7,140-point Hamiltonian grid and the 1,323-point TPMS grid were both
executed at 120 cubed. The original model functions and finite-domain
conventions were retained. No display filtering or abstract fallback was used.

| Outcome | Hamiltonian | TPMS |
|---|---:|---:|
| Requested and recorded parameter points | 7,140 | 1,323 |
| Normalized subcubic spatial or isolated-vertex evaluations | 7,116 | 1,001 |
| Higher-valence fixed-diagram evaluations | 5 | 73 |
| Per-cell time budget reached | 18 | 178 |
| Existing cavity route recognized but not revised/recomputed here | 1 | 71 |

All nontrivial subcubic successful evaluations require exact polynomial
agreement between two distinct generic views of the retained embedded graph.
These are two views using the existing polynomial implementation, not two
independent polynomial algorithms. The full records and source identities
remain unchanged after the event analysis. Detailed family counts are in
`full_grid_execution_2026_09_19.md` and
`completed_adaptive_tpms_summary.json`.

## Artifacts

The complete event-replay artifact from run `35415827400` is:

- Name: `complete-replayed-tpms-events`.
- Artifact ID: `10576140399`.
- ZIP bytes: `1285545`.
- SHA-256: `ee6d60545d409c5ef127a654fa9a4e4d650a680d061b65b77fff3aff3ac54bf1`.
- Exported attachment: `KnottedGraph_complete_TPMS_event_certificates.zip`.

It contains the 650-test JUnit output, the full event records with retained
local certificates and fresh replay results, exact-branch witnesses, the
re-audited adaptive leaf partition, and input/output hashes. Full original
regularity proof trees remain in the 20 `refined-tpms-shard-N` artifacts from
run `35411858561`; they are not all duplicated in this event ZIP.

Other exported results are:

- `KnottedGraph_complete_numerical_review.zip`: complete original-grid
  records, CSV tables, adaptive-area summaries, provenance and data hashes.
- `KnottedGraph_complete_7140_Hamiltonian_1323_TPMS_records.zip`: both complete
  original-grid record arrays and their exact-coverage manifest.
- `KnottedGraph_exact_TPMS_critical_branches.zip`: executed symbolic identities
  and the associated test output.

The Actions artifacts have 30-day retention and are not a permanent data
repository. The source code, result reports and immutable revision identifiers
are committed; the larger masks, graphs and full collapse witnesses are in
the 20 full-map shard artifacts.

## Remaining scope

The complete numerical scans are finished, but some polynomial evaluations
remain unavailable as reported above. The proposed separate Hamiltonian
polynomial-retry commit was rejected by the tool and was not applied or run.
The attempted figure-restoration write was also rejected; no replacement
main figure was generated, and those edits were not bypassed by another route.

No current 42-page manuscript source archive was available in this branch or
the accessible attachments. No final manuscript or figure PDF was rebuilt.
Local execution/stat requests continued to return transport timeouts; the
validated numerical work ran through the recorded GitHub Actions workflows.

The voxel collapse witnesses prove retractions of the represented voxel
complexes. They do not by themselves prove analytic-source-to-voxel
correspondence, and the new event certificates do not supply that separate
correspondence. The remaining 366 parameter rectangles are unknown, not
inequivalent. No complete classification into distinct source-solid isotopy
classes is claimed.
