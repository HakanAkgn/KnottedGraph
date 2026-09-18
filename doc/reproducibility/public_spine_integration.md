# Public solid-spine integration: executed results and scope

## Version and boundaries

The tested source revision is `8e5637c073e938d7099677168ebb064bc61b46bc` on
`research/solid-spine-maps-41dc87c`. This note is a documentation-only follow-up.
Work started from `782a745b7b79a7c41dc993de5cffb9ab8bba0fc5`, which already
contained the opt-in cubical collapse implementation. The protected
`integration/arbitrary-knot-fields-final-audit` branch remains at
`fab2673f665796c65a32ec42efdc1bf9b6c37678`.

Only public reconstruction, scan reporting, positive continuum equivalence,
and their tests/execution records were changed. The all-word formulas,
existing cavity-processing implementation, original manuscript and artwork
were not edited. This is not a finished manuscript or a regenerated full
Hamiltonian/TPMS phase-map release.

The prior-attempt ledger, supplied arXiv audit and Codex figure-comparison
history were reviewed before this implementation. The correction is not a
larger abstract-contraction search, a small-region filter, or a replacement of
failed spatial calculations by abstract polynomials.

## Public reconstruction is now connected to the witnessed path

`NodalSkeleton.skeleton_graph` delegates to the shared `_reconstruction.py`.
The ordinary guarded mode supplies expected component and cycle counts,
checks mask-to-thinning agreement, preserves components through cleanup, and
checks the final graph. Its default smoothing tolerance is zero. Explicit
nonzero smoothing in guarded mode remains unproved geometrically and is
labelled accordingly.

The new `reconstruction="cubical"` route calls `cubical_retract` on the actual
source mask and replays every elementary collapse before returning a graph.
It refuses external skeleton substitution and nonzero geometric smoothing.
The full collapse witness is available as `model.spine_certificate`; graph
metadata records its hash, method, coordinate convention, boundary contact,
and exact scope. Degree-two suppression concatenates original grid points
rather than straightening the embedded paths. Component/cycle agreement is
checked in addition to, not instead of, witness replay.

`NodalPhaseScan` now defaults to this cubical route and calls the existing
audited spatial evaluator. Scan records retain certificate identity and
reconstruction/evaluation metadata. They do not embed the potentially large
full collapse witness; use the reconstruction object or the supplied explicit
source-mask audit runner when archiving those witnesses.

The represented input solid is the union of closed dual voxel cells clipped
at the original sample endpoints. The public graph retains index coordinates;
the existing positive affine conversion to physical coordinates remains
unchanged. No periodic face identification is introduced.

This establishes a replayed deformation retraction of the specified voxel
complex when a graph endpoint is reached. It does not establish the missing
analytic-field-to-voxel correspondence, manifoldness or a regular-neighborhood
certificate. The source metadata explicitly keeps those stronger claims false.

## Failed calculations are no longer measured phases

`NodalPhaseScanResult.phase_grid()` reserves label `-1` for unavailable
calculations and excludes them from its label dictionary. A missing polynomial
is unavailable even when a legacy record lacks an error message. Duplicate,
missing and off-grid records are rejected rather than silently overwritten.
Input parameter arrays must be finite, distinct and increasing.

`transition_intervals()` neither bridges unavailable cells nor treats them as
transitions. It compares only compatible evaluation/normalization/valence
scopes. A returned change is explicitly an `evaluated_spine_signature_change`,
not a certified transition of the analytic source solid. Equal polynomials are
not used as embedding-equivalence certificates.

## Cache provenance

A reusable reconstruction is bound to its source-mask bytes, source shape,
reconstruction settings, actual node/edge geometry, reconstruction metadata
and complete witness. Mutating a returned graph, its claims, or its witness
invalidates the cache. Identical bytes with different array shapes cannot
share a cached reconstruction. Guarded reconstruction recomputes thinning
from the current source mask rather than trusting a stale cached thinning.
Clearing the model cache also removes the cached planar diagram and witness.

## Positive TPMS equivalence paths

`core/continuation_atlas.py` accepts a parameter block only after both the
primary and rational-trigonometric full certificate replays succeed for a
caller-supplied trusted `FieldProblem`. It rejects a changed analytic source
or spatial domain. Closed blocks are connected only by actual intersection,
without tolerance-based gap bridging. Positive comparisons include an explicit
piecewise-linear parameter path through verified blocks and identify the
certificate supporting every segment.

A disconnected or uncovered comparison returns `unknown`, never `inequivalent`.
Cover-component labels are not counts of distinct topological classes. The
atlas never propagates a chosen graph spine's polynomial to other points.

All original nine TPMS corridors at the retained `c=0.105` grid entry were
recomputed and accepted after both replays in the final run. Three belong to
each family. This reruns existing source-solid results and exercises the new
atlas integration; it is not an expansion of the earlier 914/1200 full-grid
coverage, which was not recomputed here.

## Actual execution

The complete final workflow is:

https://github.com/HakanAkgn/KnottedGraph/actions/runs/35383919300

Job `105726248033` completed successfully on the tested source revision using
Ubuntu 24.04.5, CPython 3.13.15 and the frozen repository lockfile. The log records:

- 582 pytest tests passed in 47.66 seconds; no tests were removed.
- The repository consistency check passed.
- Scoped Ruff checks passed for all newly integrated production modules and
  the new regression tests/runner.
- Nine freshly generated TPMS corridor certificates passed primary replay,
  rational replay and positive-path comparison.
- Four actual public `NodalPhaseScan` evaluations at 120 cubed completed with
  spatial or isolated-vertex results; none used an abstract fallback.
- The checked-out tracked source remained unchanged during validation.

The new tests cover unavailable-cell labels, invalid grids, incomparable
invariant scopes, component preservation, boundary-contact interpretation,
forbidden smoothing of certified graphs, stale masks, altered geometry and
witnesses, corrupted continuum certificates, disconnected covers and exact
contact versus a real parameter gap.

An intermediate run (`35380968665`) had all 566 then-present tests passing,
but failed the lint gate on two ambiguous variable names. Those names were
corrected; the gate was not disabled. The native-class edit was committed only
after the subsequent complete suite and lint passed (`35381621077`). Its
commit is `fdf4907685f570a10442a3dfb2e0495a6868a66e`.

### Public-API Hamiltonian checks

These are four explicitly selected regression cases, not 7,140 map samples.
The finite box and original model functions are retained.

| Family | lambda | Energy | Replayed collapses | Boundary contact | Result |
|---|---:|---:|---:|---|---|
| Hopf to trefoil | 0 | 0.3 | 99,809 | No | `-A**4-A**3-2*A**2-A-1` |
| Hopf to trefoil | 0.5 | 0.3 | 163,141 | No | `-A**4-A**3-2*A**2-A-1` |
| Hopf to trefoil | 0.5 | 1.8153061224489797 | 1,310,566 | Yes | `-1` |
| Trefoil to cinquefoil | 1 | 0.3 | 156,147 | Yes | `-A**8+A**7-3*A**6+2*A**5-4*A**4+2*A**3-3*A**2+A-1` |

The high-energy singleton is explicitly labelled `isolated-vertices`; the
other three results are normalized subcubic spatial-graph evaluations. Equal
values in the first two rows do not prove source-solid equivalence. Nongeneric
projection views were rejected while valid alternative views were evaluated;
these warnings are retained in the run log.

### Separate six-case 120-cubed source-mask audit

Run `35381621077` also executed `dev/run_solid_spine_audit.py` at dimension120.
All six selected voxel complexes reached graph endpoints and passed complete
collapse replay and independent component/cycle checks. These cases overlap
with three public-API cases above and must not be counted as ten independent
source geometries.

| Case | Source (beta0,beta1) | Graph vertices/edges | Collapses |
|---|---|---|---:|
| P_to_D_start | (1,5) | 8/12 | 837,404 |
| Hopf_start_low | (1,2) | 2/3 | 99,809 |
| G_to_P_middle | (1,5) | 8/12 | 866,880 |
| Hopf_trefoil_middle_high | (1,0) | 1/0 | 1,310,566 |
| G_to_D_start | (1,8) | 14/21 | 751,489 |
| Cinquefoil_endpoint | (1,4) | 6/9 | 156,147 |

The Gyroid endpoint at c=0.105 had beta1=14 in the prior 64-cubed audit and
beta1=8 here. The low-energy cinquefoil case changed from (beta0,beta1)=(9,1)
to (1,4). Both resolutions have valid retractions of their own digital input.
These differences are evidence that extraction correctness alone does not
establish convergence or analytic-source correspondence. They are not hidden
by filtering or retrospectively reclassified as extraction failures.

## Data artifacts

Final public integration/TPMS artifact: `10563875223` (21 files, 63,140 bytes),
SHA-256 `83dca263cf142be068e971b2e5ed3fc83848b361526d1643378c6a8f6b95ac44`.
It contains the JUnit results, nine complete compressed continuum certificates,
plans, source hashes, equivalence paths and four public scan records.

Full six-case source-mask/graph/collapse-witness artifact: `10563101722`,
SHA-256 `f02048e7e1a3459ea90f1474998f75bc83f428d611ecd75d07253719b94f62b4`.
It includes the actual sampled masks and grids, complete graph geometry and
compressed collapse witnesses for the 120-cubed audit. These artifacts were
also exported as review ZIP attachments. Actions retention is 30 days, not a
permanent research-data archive; retain the downloaded copies or regenerate
using the commands below before relying on them for a paper release.

## Reproduction

From a clean checkout of the tested code with native build prerequisites:

```sh
uv sync --all-extras --group dev --group docs --frozen
uv run --no-sync pytest -q
uv run --no-sync python dev/check_repository_consistency.py
uv run --no-sync python dev/certify_tpms_continuation.py \
  --mode historical-nine --jobs 2 --out /tmp/kg_tpms_corridors_new
uv run --no-sync python dev/verify_public_spine_and_atlas.py \
  --tpms /tmp/kg_tpms_corridors_new --out /tmp/kg_public_spine_new \
  --dimension 120 --case-seconds 180
uv run --no-sync python dev/run_solid_spine_audit.py \
  --out /tmp/kg_six_spines_new --dimension 120 --seconds 850 --case-seconds 140
```

Use fresh output directories. Computation times and timestamps vary with the
machine; replay conditions and stored source identities are the relevant
checks. A budgeted run may report unfinished cases rather than fabricated
completed records.

## What this release does not claim

No new analytic-field-to-discrete-solid certificate was integrated or tested.
No full 7,140-cell Hamiltonian rerun or full TPMS equivalence classification was
completed. No original main figure, representative-geometry strip, caption or
manuscript build was regenerated. The local file/execution environment remained
unavailable; repository builds and computations were instead executed through
the recorded branch-local CI runs. The tested API repairs are not a substitute
for those remaining data and manuscript deliverables.
