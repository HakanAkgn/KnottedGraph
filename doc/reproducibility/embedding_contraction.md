# Embedding-sensitive contraction: conservative design

## What can be established

Yes: a contraction can retain the embedding information relevant to the **regular neighborhood** of a spatial graph. It changes the chosen graph spine, so it is not an ambient isotopy between the original and quotient abstract graphs. A witnessed local contraction, followed by ambient isotopies, preserves the represented handlebody-link. Conversely, abstract graph contraction alone contains no information about how the tubes are embedded.

The primary reference is [Ishii–Iwakiri, Theorem 2.3 and Figure 2.4](https://doi.org/10.4153/CJM-2011-035-0). Their contraction is a local spatial move; equality of abstract graph quotients is not their equivalence criterion. Their Proposition 2.4 separately requires an invariant to be contraction-invariant before using it as a handlebody invariant.

## Why the current plot comparison is insufficient

The archived `tpms_plotly_region_geometry.py` plot generator discards node positions and edge polylines in `graph_from_geometry` (lines 86–93). `contract_edge` (152–170) retains incidence only. `can_contract_to` accepts plain NetworkX isomorphism (203, 222), deduplicates states by abstract isomorphism (234–238), and represents both failed bounded search and proved obstruction with `False`. `contraction_phase_map` then takes transitive closure of those positive abstract matches. Thus its classes are **abstract contraction matches**, not certified spatial/handlebody equivalence classes.

The simplest counterexample is a circle and a trefoil, each stored as one vertex plus one self-loop. Their incidence is identical, so the existing comparison immediately returns true. Their solid-torus neighborhoods are differently embedded. An analogous multi-component counterexample is an unlink versus Hopf link stored as two self-loops. Correct code must retain geometry even when the abstract graphs agree.

## Implementable conservative move certificate

Initial prototype: contract only an already straight, non-loop edge `u--v`, with `v` fixed. Every polyline segment must be finite, nonzero and belong to a valid embedded PL graph. Interpret binary floating-point input values as exact rationals for the geometric predicates. This certifies the supplied PL graph, not uncertain physical measurements.

For each segment `u--p_i` incident to `u` other than the contracted edge, form the closed swept triangle `D_i=conv(u,v,p_i)`. This includes both half-edges of a loop at `u`. Require:

1. Every `D_i` is nondegenerate. Degenerate/collinear fans are rejected by this prototype, even when a more sophisticated move could work.
2. The fixed graph intersects `D_i` only at the prescribed attachments `p_i` and/or `v`. A contact is permitted only when the segment incidence identifies that point as its legitimate endpoint; geometric coincidence is insufficient.
3. Distinct fan triangles intersect only along the contracted segment `[u,v]`.
4. The input and resulting graph have no intersections except prescribed shared endpoints. Check the entire PL geometry, including parallel graph edges, self-loops, and different components.

Use the explicit motion `u(t)=(1-t)u+t v`; each moved first segment is `[u(t),p_i]`. For `0<=t<1`, the empty-fan conditions prevent strand crossing. At `t=1`, the sole allowed edge collapse is `[u,v]`; all other half-edges attach to `v` without overlap. This is a geometric witness of a neighborhood-preserving spatial edge contraction. Store all unchanged polylines and all first-segment replacements; never redraw the remaining graph with a layout algorithm. A parent implementation can later extend the accepted moves through independently certified empty-triangle simplifications of bent contraction edges.

Exact intersection predicates are needed for a mathematical certificate. Merely sampling the motion at finite time steps or checking the output for intersections does not prove the path is safe. With ordinary tolerance-based geometry, label the result `numerically_audited`, not certified, and classify near contacts as unknown. A clearance threshold is useful to protect physical uncertainty but is separate from exact PL topology; incident segments necessarily have zero distance at their legitimate common vertex.

## Comparing two reconstructed spines

Each search state must preserve the embedded geometry. Safe contractions are only the first half of a comparison. To establish equivalence, two descendants must additionally be connected by a **witnessed ambient isotopy**: exact coincidence under a proven incidence mapping, an orientation-preserving rigid motion followed by exact matching, matching complete oriented planar-diagram structure, or a verified sequence of allowed spatial-graph moves. Matching only the adjacency matrix is never a terminal success condition.

A minimal first implementation may accept only exact geometry coincidence after contractions; it will deliberately leave many truly equivalent pairs unknown. A later robust diagram-isomorphism/diagram-move backend would be much more useful for TPMS data. Matching PD records requires the cyclic order at each crossing, the over/under pairs, vertex records and component embedding data; simply sorting PD rows is not a proof of equivalent planar diagrams. A reflection must not be silently accepted as an ambient isotopy.

Contractions often create vertices of valence above three. A selected-diagram Yamada polynomial can still be computed, but the manuscript's projection-independent subcubic interpretation then no longer applies. Either preserve the relevant ribbon/rotation data and equivalence convention, or return the polynomial as a diagram diagnostic only. An IH-move search on trivalent spines is a possible future alternative.

## What Yamada can and cannot do here

- Different valid normalized Yamada polynomials distinguish subcubic **spatial graphs of the stipulated equivalence type**.
- Equal polynomials do not prove equivalent spatial embeddings.
- Yamada is not, in general, invariant under spine contraction. Therefore different Yamada values of two independently chosen spines do **not** prove their handlebodies inequivalent. Do not use this as a global negative filter on contraction-equivalence search.
- Equal descendant Yamada values are an inexpensive candidate filter, never a reason to union phase classes.
- A genuine handlebody invariant, e.g. a justified quandle/cocycle invariant, can provide a valid negative test. Different component counts and different handlebody-genus multisets also rule out equivalence.

## Bounded search and output semantics

Use a three-way result with a separate reason:

- `equivalent`: a replayable sequence of certified contractions and terminal embedding equivalence is present;
- `inequivalent`: a valid invariant of the **intended equivalence relation** obstructs equivalence;
- `unknown`: no witness within the chosen move set, depth/state/time limit, invalid input, unsupported bent edge or numerical ambiguity.

A fully exhausted contraction-only search still gives `unknown`: equivalence can require expansions or isotopies before contraction. Cap pair attempts, state count, depth and geometric predicate count explicitly, and retain those caps in machine-readable output. Deduplicate by exact serialized embedded states (or a separately certified isotopy), not by abstract isomorphism or Yamada equality. Union-find may combine `equivalent` edges because witnessed equivalence is transitive, but every merge must retain its witness chain. An absence of merges is not evidence of distinct topological phases.

## TPMS experiment before adopting this method

1. Load original reconstructed `pos/pts` geometry. The existing export rounds to six decimals and downsamples polylines (`graph_geometry`, lines 752–789 of the compact scan script). Compare `original_point_count` with stored counts and record any losses. A certificate on those exported graphs only applies to those exported PL graphs unless a separate certificate links them to the originals.
2. Report input validity, component/genus signatures, number of eligible straight edges, accepted local contractions, rejected/unknown reasons, and resulting valences.
3. Start with representative neighboring parameter cells and same-old-class pairs. Keep the unmerged raw labels visible. Present certified merges separately from abstract matches and invariant matches.
4. Record all witnesses and unresolved comparisons. A conservative prototype may find few or no certified cross-cell matches; that is an honest result and a basis for deciding whether a diagram-level backend is worth implementing.

## Required meaningful tests

- Accept a clean two-junction tree contraction, preserve edge/node provenance and loops/multiedges, and replay witness.
- Reject a foreign strand piercing a swept triangle even when both input and proposed final graphs are intersection-free.
- Reject fan overlap, unintended endpoint contact, loop contraction, nonstraight selected edge, invalid input intersections and budget exhaustion.
- Check a foreign strand just off the swept plane is accepted by exact predicates; a strand on the plane is rejected.
- Same abstract graph embedded as unknot/trefoil and unlink/Hopf must never receive `equivalent` solely from abstract isomorphism or equal cycle rank.
- A bent/subdivided edge must remain unknown in this first move set; no opportunistic simplification is allowed without its own certificate.
- Check contracted subcubic endpoints may produce a four-valent vertex and tag the invariant scope appropriately.

The proposed exact-fan certificate is intentionally sufficient rather than necessary. It is a small, auditable move primitive, not a complete handlebody equivalence solver and not a replacement for validating that a voxel reconstruction is itself a spine of the intended volume.

## Prototype findings (18 September 2026)

The self-contained experimental module `embedded_contraction.py` implements the exact-fan certificate and a bounded two-sided comparison whose only successful terminal condition is identical embedded PL geometry. It is available in the experimental module `knotted_graph.core.embedded_contraction` and is not used to merge the revised phase map. Fourteen meaningful tests pass in the revision checkout's Python environment, including the crossing/near-contact and different-knot/link adversarial cases above. Rational coordinates that would be rounded when emitting float polylines are refused.

All 1,323 archived TPMS geometries were inventoried. Of these, 304 contain at least one eligible straight non-loop edge, with 526 such edges total. Exactly one file reports downsampled geometry: `gyroid_to_schwarz_p_lambda008_c014.json.gz` (one edge). Every export still uses coordinates rounded to six decimal places, so these certificates are explicitly about the archived display polylines.

Twenty contraction attempts were selected across twelve geometries: four representatives per interpolation family, spread across the eligible geometries' segment-count range. All twenty passed the exact empty-fan predicate and full input/output embedding checks. The twenty attempts took approximately 3.15 seconds total; the slowest took 0.76 seconds. Their resulting maximum valences ranged from four to nine. Exact Yamada evaluation of these high-valence descendants must therefore retain its diagram/rotation convention and cannot simply inherit the manuscript's subcubic projection-independent interpretation. The replayable witnesses are in `embedding_contraction_tpms_probe.json`, and the full inventory is in `embedding_contraction_geometry_inventory.json`.

Nine neighboring parameter pairs at `c=0.105` were then chosen from the same old stabilized contraction class, three per interpolation family. All nine pairs have abstract-isomorphic graphs. With depth one, at most eight processed states, at most sixteen move attempts and 100,000 exact predicates per move, the result was:

| Outcome | Pairs | Interpretation |
|---|---:|---|
| Exact embedded equality | 1 | Two already identical one-vertex/no-edge spines in Schwarz-P→Diamond at λ=0.50 and 0.55; zero contraction moves needed. |
| Unknown | 8 | No certified common embedded descendant within this deliberately limited move set; this does not prove inequivalence. |

The pair outcomes and rejection reasons are in `embedding_contraction_pair_probe.json`. No nontrivial cross-parameter merge was certified. The useful result is therefore twofold: the local embedding-preserving contraction primitive is feasible on real TPMS outputs, while the existing nontrivial phase-class mergers remain unverified without a stronger witnessed isotopy/diagram comparison stage.

The module is installed at `src/knotted_graph/core/embedded_contraction.py`. It is intentionally not used by the production signature-grouping path. The revised figure shows unfiltered voxel homology and computational status; no nontrivial spatial equivalence is inferred from these experiments.
