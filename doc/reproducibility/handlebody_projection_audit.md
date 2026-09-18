# Projection and handlebody benchmark diagnosis

The historical row 4923 (`random_4920`) mismatch is a projection serialization error on the saved recovered graph. It is not evidence that voxel skeleton recovery changed this embedding's Yamada polynomial. This diagnosis uses the exact recovered payload whose SHA256 is `0adffaddbfeee0e46f126e358da01dd19f77a7b3c60c79527d24fd814bd359f6`.

## Two concrete defects and corrections

1. `PDCode._update_incidences` previously took the next coordinate after a Shapely substring endpoint. When the cut falls on an existing sample, Shapely can emit two coordinates differing only by roundoff. In the failing six-crossing view, recovered crossing 3's outgoing arc 10 had a displacement around 5.6e-17. Its computed angle became 0 instead of the actual outgoing tangent, turning `X[10,7,9,6]` into `X[10,9,6,7]`. The latter joins the wrong opposite half-edges. The fix consistently measures endpoint directions from each arc's own endpoint and skips samples within 64 machine epsilons of its coordinate scale. It still rejects an arc with no resolved direction; it does not skip a genuinely distinct segment of length 1e-10. Height/half-edge matching uses the same direction helper.

2. `generate_isotopy_angles` emitted `(yaw,pitch,0)`, but the existing uppercase ZYX rotation matrix is `Rz(yaw) Ry(pitch)`. Its last row is `(-sin(pitch),0,cos(pitch))`: every view lies in the x-z plane, while yaw merely spins the picture. Thus the advertised hemisphere sample was actually a restricted set of coplanar directions. Vertical planar portions of the generated curves can overlap in every such view. The fix constructs a proper camera rotation whose last row is each Fibonacci hemisphere direction and decomposes that matrix into the requested Euler convention. All 24 proper intrinsic/extrinsic conventions are checked. Explicit rotation-matrix behavior is unchanged; documentation is corrected to its existing standard convention (uppercase intrinsic, lowercase extrinsic).

Neither fix chooses a projection by agreement with an expected polynomial. Selection still uses the valid sampled diagram with the fewest crossings. Nongeneric views remain rejected.

## Recorded mismatch reproduced and corrected

Before the endpoint fix, default six-sample selection reproduced the archived incorrect degree 13 polynomial, while explicit XY and two generic views returned the seed's degree 14 polynomial on the *same* saved recovered graph. After correction, default selection and all four documented explicit views agree with the seed. The correct normalized value is

`A**14 - A**13 + 2*A**12 + 5*A**11 - 10*A**10 + 25*A**9 - 29*A**8 + 40*A**7 - 34*A**6 + 35*A**5 - 20*A**4 + 15*A**3 - 2*A**2 + A + 2`.

The exact recovered fixture, source revision, environment, graph hashes, view angles, and results are retained in `handlebody_failure_final.jsonl` and `tests/projection/data/handlebody_random_4920_recovered.json`.

## Historical data and provenance

The earlier checkout's CSV has 5,000 rows: 4,472 matches, 1 mismatch, 527 blank comparisons. All 527 blanks have recorded projection-failure errors:311 seed-only failures,204 both sides,12 recovered-only. All 527 were still marked `overall_pass=True` because unavailable comparisons were omitted from the final condition. A missing invariant comparison does not establish invariant preservation.

The earlier recovered archive declares commit `49e34e47ca5c182f55ef5f0ea0906220df59befb`, 5,000 cases. The new checkout's saved CSV/recovered archive has 100 cases and commit `62b4ac1ea81c7b51cd71b934d1e32f82bca23f72`; it is a different saved run. These original artifacts have not been overwritten.

The earlier seed archive and CSV are **not a uniformly aligned batch**. Of the first 100 cases, 96 disagree in stored seed structure or certification metadata; all remaining 4,900 agree. Among the 96, 94 differ structurally, and `random_0027` and `random_0036` have matching abstract-graph metadata but different radius/spacing/occupied-voxel certification. The three controls and `random_0000` are the only compatible rows in the first 100. Hence 96 comparisons are withheld; a case name is insufficient to associate inputs. Per-row checks and exact archived seed hashes are saved in `handlebody_seed_alignment.json`.

Historical rows did not store exact seed payload hashes. Agreement of metadata is recorded as compatibility, not independent proof of the original seed coordinates. Every new replay records exact hashes for the actual seed and recovered payload it evaluated.

## Reruns

- All 527 historically unavailable comparisons were replayed using corrected projections:518metadata-compatible seed/recovered pairs match; 9 have incompatible seed metadata and their comparisons are withheld. Every recovered graph is now evaluable. No new mismatches occurred in the 518 compatible pairs. See `handlebody_527_corrected_projection.jsonl`.
- A fresh recovery from the first 20 immutable archived seeds, including voxelization, volume Betti checks, skeletonization, extraction, cleanup, and both invariant evaluations, passes 20/20. This is a new run of those exact archived seeds; it is not a repair or reinterpretation of the incompatible historical first 100 rows. Total measured computation 21.0 seconds, maximum 1.81 seconds per case. Every new recovered embedding and both hashes are in `handlebody_fresh20.jsonl`.
- The complete 5,000-row guarded archived-pair replay finished: all 4,904 metadata-compatible pairs match; 96 conflicting seed associations are withheld. All 5,000 recovered payload hashes agree with the CSV and all 5,000 recovered graphs evaluate successfully. No compatible seed evaluation fails. Its 527 cached entries were reused only after verifying input files, graph payload hashes, projection/invariant source hashes, and dependency versions. Exact per-case results and terminal counts are in `handlebody_all5000_corrected_projection.jsonl`; compact checks are in `handlebody_all5000_summary.json`.

## Reproducibility changes

`dev/audit_handlebody_recovery.py` accepts explicit `--seed-archive`, `--recovered-archive`, `--csv`, and `--case` inputs. It refuses existing output files, stores source/input/environment evidence, withholds unpairable comparisons, and supports at most four workers for archived polynomial replay. `--rebuild --notebook ... --limit 20` performs a serial fresh recovery using only the named audited function definitions from notebook 04, without executing notebook run cells.

Notebook 04 now reuses a checkpoint only when seed-payload hash, exact recovered-payload hash, and source/dependency fingerprint all agree. Legacy rows without these identities are recomputed. It distinguishes `match`, `mismatch`, and `unavailable`; the combined validation cannot pass when Yamada is unavailable. This prevents stale rows from silently becoming evidence for newly generated inputs.

Focused tests cover numerical endpoint duplicates, preservation of small but resolved directions, degenerate arcs, the archived failure, all 24 rotation conventions, invariance across sampled views for small subcubic examples, and stale checkpoint invalidation. Latest consolidated check: 73 tests passed (projection, spatial-PD regression, and notebook checkpoint tests).

Two final documentation corrections changed projection source-file bytes after the run but did not change executable code. `projection_documentation_provenance.json` records both file hashes, exact diffs, and matching AST hashes after docstring removal. No numerical replay was required for these wording corrections.

Portable inputs and run commands are in `User_guide/benchmarks/results/arxiv_revision/handlebody_inputs/README.md`; exact result evidence is in `User_guide/benchmarks/results/arxiv_revision/validation/`. Absolute paths inside original provenance records identify historical inputs and do not need to exist when running the portable commands.
