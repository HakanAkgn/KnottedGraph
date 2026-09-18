# Prior-attempt ledger

Base: fab2673f665796c65a32ec42efdc1bf9b6c37678.
Separate branch: research/continuum-isotopy-fab2673.
The integration/arbitrary-knot-fields-final-audit branch and original manuscript
are not modified. The suggested performance branch points to older commit
18c220d7604b73736258f709a3f55bb792d04667 and is not overwritten.

## Sources reviewed

The supplied Codex conversation and full ten-point arXiv audit were read. The
original manuscript's TPMS definitions, spherical restriction, critical-point
and constrained-boundary equations and both original phase-map figures were
examined. The audited repository's extraction, projection, contraction,
Hamiltonian, timing and formula-discovery records were consulted. The actual
TPMS constructor and nine-pair selection were checked against the implementation.
These sources summarize prior work; unavailable internal Codex execution logs,
historical model conversations and the separate revised 42-page local paper
were not recovered or invented.

## Earlier approaches and the decision made here

| Earlier attempt | Limitation documented by the audit | Decision |
|---|---|---|
| Small-region label filtering | 171 historical cells changed, order dependence, hidden errors | Never use filtering as topology evidence. |
| Abstract contraction comparison | Embedding information discarded | Compare the defining spatial solids instead. |
| Increasing junction-zone persistence | Repeated wrong cycle ranks, including 11 becoming 8 | Retain the audited homology guard; do not weaken it. |
| Dropping small components | Trees and isolates lost | Retain the audited component-preservation fixes. |
| Treating a cavity-bearing volume as a graph | A shell could become a point | Test the whole sublevel solid, with no graph-spine assumption. |
| Spatial-to-abstract fallback | A different polynomial substituted silently | Retain distinct outputs and explicit errors. |
| Additional old camera samples | The directions were confined to a plane | Retain the corrected hemisphere sampler and incidence fix. |
| Exact straight-edge contractions | Terminal coordinate equality too restrictive; eight unknown pairs | Do not just increase search budgets; use analytic continuation. |
| Larger voxel grids alone | 8/45 Betti tuples changed from 64 to 96 cubed | Certify analytic regularity throughout space and parameters. |
| Hamiltonian diagnostic replacement | Coarser finite-window homology, not the historical Yamada map | Do not claim this experiment repairs that separate application. |
| Removing original representative geometries | Scientific presentation was lost | Leave all originals intact; new evidence is a separate figure. |
| Archived benchmark replay | 96 seed associations unresolved | Do not invent seeds or call replay a fresh reconstruction. |
| Finite word checks | Missing family-specific all-word derivation | Do not relabel tests as a universal proof. |

## New method and falsification tests

The manuscript already gives the two possible event systems for a smooth
sublevel set in a fixed ball: a bulk zero of the spatial gradient on the level,
and a constrained critical point on the spherical wall. We turn the absence
of both events into a whole-domain interval certificate. This checks an
analytic parameter path directly rather than comparing independently extracted
spines. Standard regular-level continuation and interval topology are not new
mathematics; the implemented, replayable application is the new work here.

Before accepting TPMS results, tests cover interior births, wall tangencies,
box-corner events, a singular path with equal endpoint topology, incomplete
coverings, tampered certificates, unsupported expressions, budget exhaustion,
cavity-bearing shells and independently evaluated interval bounds.

## Execution history

A first complete run produced nine successful corridors and 914/1200 successful
parameter rectangles. An additional CLI nine-case rerun was interrupted after
eight cases and was never counted as complete. The working runtime later
reinitialized before the package was saved; the module tree already existed on
GitHub, while the prior figure preview and tool logs remained available.
The two modules and runner were reconstructed byte-for-byte, matching their
recorded SHA-256 hashes. The delivered results were then recomputed in the
current runtime, rather than fabricating the lost proof trees. The delivered
plans identify their actual environment. The preview is retained only as an
illustration of the first run, not as evidence replacing the fresh certificates.
