# Scientific computation revision

This revision starts from collaborator branch `codex/arbitrary-knot-user-integration` at `43620e2`. The requested integration branch at `643fef8` is an ancestor. The changes address discrepancies found when checking the manuscript against retained code and data. Historical artifacts are not silently overwritten.

## Changes in behavior

- Spatial Yamada evaluation no longer substitutes an abstract graph polynomial when projection fails. Cell records preserve the error, selected projection, normalization and valence-dependent interpretation.
- Subcubic normalized spatial-graph polynomials, higher-valence fixed-diagram polynomials, abstract graph summaries and unavailable computations have separate labels. Equal signatures are not equivalence certificates.
- Extraction can require an expected component count and cycle rank. Candidates must meet both constraints, including any zero-radius fallback. Tiny acyclic voxel components and isolated graph components survive tracing and pruning.
- TPMS volume checks include enclosed cavities and the second Betti number. A cavity-bearing volume is unsupported by a graph-spine reduction; failure no longer becomes an artificial point.
- Projection tangents skip numerically duplicate samples near endpoints. Camera sampling now spans a hemisphere rather than repeatedly sampling one plane.
- The optional TPMS display filter is disabled by default. When requested, it uses synchronous passes, only large-component recipients, unique boundary-contact votes, ties left unchanged and protected error labels. It stores the raw grid and changed-cell mask. Filtering is not a measured classification.
- Benchmark checkpoints bind both input graphs and the source/environment identity. Unavailable comparisons are reported separately from matches and mismatches.

## New measurements

The unfiltered TPMS run contains 1,323 cells on 64-cubed spatial grids. Of these, 74 have nonzero second Betti number and remain unsupported. All 1,249 other extracted graphs match their voxel component and cycle counts. There are 562 normalized subcubic evaluations, 51 higher-valence diagram evaluations and 636 abstract summaries. Their combined signature counts are bookkeeping rather than topology classes.

Independent checks of 45 family/parameter entries at 24, 64 and 96 cubed show changes in 13 Betti tuples; eight still differ between 64 and 96. This prevents a claim that the continuum topology map has converged. Parameter subsampling and the archived filter sensitivity are recorded separately.

The deterministic family-law audit provides exact coefficients for 324 pure-braid words, complete 15-state matrices, a finite-block rank certificate, cubic identities and noncommutation. These finite calculations do not prove the proposed all-word graph-family correspondence. All 459 mixed-family tests and two length-101 pure-braid extrapolation cases also pass exact comparisons. The other 18 planned long-word cases were not completed. Original LLM conversations and historical freeze timestamps were not recovered.

## Experimental spatial contractions

`knotted_graph.core.embedded_contraction` implements a deliberately limited straight-edge contraction using exact rational intersection predicates on the supplied floating-point PL geometry. It records a replayable local motion witness. The bounded comparator requires identical embedded terminal geometry, not abstract isomorphism, for a successful match. Exhausted or unsupported searches return `unknown`.

Twenty tested moves on archived TPMS polylines passed. Nine neighboring-pair searches produced eight unknowns and one already-identical point pair. That point pair was from a legacy cavity-bearing reconstruction, so it does not validate the original volume reduction. No nontrivial phase merger was established. The prototype is not used to merge the revised TPMS map.

## Reproduction

Install the checked lockfile with `uv sync --all-extras --group dev --group docs --frozen`. Run `uv run pytest -q` and `uv run python dev/check_repository_consistency.py`. Exact numerical commands and hashes are stored next to each result under `User_guide/applications/results/arxiv_revision` and `User_guide/benchmarks/results/arxiv_revision`.

The manuscript working copy lives in the ignored local `paper/` directory so that reviewing the new scientific conclusions does not overwrite the original Overleaf draft or publish the manuscript with a code push. Its build uses `latexmk -pdf -no-shell-escape Third_draft.tex` from that directory.

The final code validation comprises 465 passing pytest tests, the configured Ruff gate, the repository consistency audit, both invariant-sanity/application-regression notebooks, and 12 extracted-function Hamiltonian checks. The dated revision manifest records data hashes and the audit environment.
