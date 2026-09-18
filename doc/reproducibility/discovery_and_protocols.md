# Reproducibility audit for the arXiv revision

This audit records a new exact reproduction on 18 September 2026. It does not reconstruct the historical discovery process. The original manuscript, repository release and author metadata were not edited by this audit.

## Exact transfer reconstruction

The supplied `audit_discovery_artifacts.py` reads the pure-braid constructors and raw invariant evaluator from cells 68, 69, 71 and 72 of `User_guide/applications/05_yamada_formula_discovery.ipynb`. The script checks that the imported Python package belongs to the requested checkout. Its provenance records the Git commit, dirty-source status, notebook and selected-cell hashes, source hashes, Python and package versions, and native-backend availability.

Run from the repository with its installed uv environment:

```sh
uv run python dev/audit_discovery_artifacts.py \
  --repo . \
  --out User_guide/applications/results/arxiv_revision/discovery \
  --symbolic
```

Use a fresh output directory to preserve an earlier audit. The data and algebra are deterministic; timestamps and execution timings are not. The script currently addresses notebook cells by index and saves their hashes, so a reorganized notebook requires a reviewed update to the cell indices.

The completed audit recomputed 324 distinct short-word spatial graphs and their raw integer Laurent coefficients:

- 255 words in the union of `u+v`, `u+A+v`, `u+B+v`, where `u,v` range over the manuscript's 15-word Hankel basis. Their maximum length is nine.
- All 255 binary words of lengths zero through seven. These are a different set of 255 words; 69 are additional to the Hankel set.

Exact arithmetic verifies:

1. `det(H(2)) != 0`, hence rank 15 over the rational-function field. The exact rational determinant is retained in `certificate.json`.
2. `H H_inverse = I`, `H T_A = H_A`, `H T_B = H_B`.
3. Both operator identities `(T-Y^2 I)(T-Y^-2 I)(T-Y^-4 I)=0`.
4. `T_A T_B != T_B T_A`, with a nonzero exact entry of the commutator recorded.
5. `L_transpose T_w R` equals direct raw evaluation for all 324 audited words.
6. For all binary words through length seven, raw polynomials agree under reversal; within a fixed A/B count class, no additional equalities between distinct words occur.

The machine-readable outputs are:

| File | Contents |
|---|---|
| `short_exact_records.json` | Direct graph/PD evaluations, integer coefficient dictionaries, hashes, timings |
| `word_plans.json` | 15-word basis, all short-word lists, exact 20-word long-test plan |
| `hankel_input_matrices.json` | Exact H, H_A and H_B entries |
| `exact_transfer_matrices.json` | H, H_A, H_B, H inverse, T_A, T_B, L transpose and R over Q(Y) |
| `certificate.json` | Verification results and nonzero rank/commutator witnesses |
| `exhaustive_short_word_checks.json` | Reversal/collision checks for all words through length seven |
| `provenance.json` | Actual environment and source identity for the new run |
| `new_audit_freeze.json` | Timestamped snapshot of this audit's formula and planned long words only |

The formula payload's canonical SHA256 (excluding its hash field) is
`72ac75fe0ca14632e888207d7c39ec71b05dac0c7e5a44eafb37a27c246e63a3`.
The 255-word and expanded 324-word runs independently produced this same matrix artifact.

### Proof status

These are exact finite-data algebraic checks. They establish the stated finite-block rank, matrix identities and agreement on the explicitly listed words. They do not establish that every graph closure factors through the proposed 15-state representation. A proof of that closure/action correspondence is still needed before inferring the all-word identity or the all-word upper bound and minimality. A 15-dimensional known skein algebra alone does not supply the missing family-specific map.

A subsequent bounded audit completed **two full exact length-101 comparisons**: `A^51 B^50` and `(AB)^50 A`. In each case the directly evaluated raw Laurent coefficient dictionary equals the transfer prediction. The complete 20-word plan and the two successful records are retained in `User_guide/applications/results/arxiv_revision/discovery_long_words/`; the other 18 cases were not completed within the 240-second run. These are retrospective extrapolation checks against an already known candidate, not evidence of historical freezing or data withholding.

To repeat the bounded check into a fresh output directory:

```sh
uv run python dev/audit_discovery_long_words.py --repo . \
  --artifact User_guide/applications/results/arxiv_revision/discovery \
  --out /tmp/new_discovery_long_words --seconds 240
```

Runtime depends on hardware; the budget does not guarantee a fixed completion count. The certificate explicitly records completed and unfinished cases.

For the homogeneous and commuting mixed formulas, use the normalized formulas only for nonempty constructions (m >= 1). At m=0, the common matrix contraction is `-Y^2(1+Y+2Y^2+Y^3+Y^4)`, whose minimum degree is two. Its normalized expression is `-(1+Y+2Y^2+Y^3+Y^4)`. The zero-polynomial normalization convention must also be explicit.

## Historical records and word-plan discrepancy

The discovery notebook has no retained execution outputs. None of its expected 05a, 05b or 05c coefficient/verification ledgers were found in the checkout's application results. The full notebook does contain executable candidate expressions, a prospective order of operations, and deterministic plans. This is sufficient to rerun a procedure; it does not demonstrate when a particular candidate was proposed or what a model saw.

No original prompts, model conversation transcripts or independently timestamped historical candidate snapshots were located in the supplied project material. A Boolean named `prediction_frozen_before_holdout` written by a later notebook execution does not recover historical provenance. The text should not identify a specific model, access date, quoted prompt, or successful historically held-out run unless the corresponding original records are provided.

There is a historical plan/count discrepancy: notebook cell 60 generates **459** mixed-theta words, comprising 363 exhaustive nonempty words of lengths one through five and 96 additional words of lengths six and seven, rather than the original manuscript's 500. No historical 500-row ledger was found. A new audit has now completed **all 459 exact normalized polynomial comparisons successfully** after the constructor/contact/projection preflights. Word lists, per-word coefficient results and the certificate are retained in `User_guide/applications/results/arxiv_revision/mixed_family/`. The separate mixed-family long-word plan has not been newly completed.

```sh
uv run python dev/audit_mixed_family.py --repo . --out /tmp/new_mixed_family_audit
```

Suggested availability wording after these files are included in the release:

> The revision archive includes the family constructors, executable candidate expressions, deterministic word lists, raw coefficients for 324 distinct short pure-braid words and two length-101 words, all 459 mixed-family comparisons, the exact Hankel and transfer matrices, and exact-arithmetic verification scripts. These artifacts document a retrospective computational audit. Original model conversations and independently timestamped discovery/freezing records are not included, and no historical data-withholding chronology is inferred from the reconstructed artifacts.

The revised manuscript distinguishes exact completed finite checks, uncompleted test plans and missing historical records. None of the new finite comparisons supplies an all-family graph-closure proof.

## Hamiltonian sampling recovered from code and retained HTML

The read-only `audit_reproducibility_inventory.py` exports current notebook controls, exact rational sample arrays and retained HTML metadata. It does not execute a phase scan or certify the provenance of the manuscript PDF.

The historical notebook controls recovered before revision were:

- Lambda: 60 inclusive samples, `lambda_j=j/59`, j=0,...,59. Thus delta lambda is 1/59, approximately 0.0169491525424.
- Candidate Gamma/energy: 50 inclusive samples, `Gamma_k=0.30+(99/980)k`, k=0,...,49, ending at 5.25. Step approximately 0.101020408163.
- Spatial grid: 120 cubed.
- `N_JOBS=2` was declared, but the actual row/cache scan loop is sequential; it does not dispatch two workers. OMP/MKL/OpenBLAS/NumExpr thread counts use `setdefault`, so preexisting overrides remain possible. This differs from the new pilot, which explicitly launches two processes.
- Smoothing retries: 4.0, 1.0, 0.25, 0.0.
- Projection retry samples: 32, 96, 256.

The retained HTML `NewPhaseMapPlots/RealMaterials/html/07_hamiltonian_yamada_plotly_region_geometry_with_materials.html` includes all five analytic transitions. Their retained Gamma ranges are prefixes of the common candidate array:

| Transition | Gamma samples | Exact last candidate index | Last Gamma |
|---|---:|---:|---:|
| Hopf to trefoil | 16 | 15 | 1.8153061224489797 |
| Hopf to Solomon | 24 | 23 | 2.623469387755102 |
| Unknot to trefoil | 16 | 15 | 1.8153061224489797 |
| Unknot to Solomon | 24 | 23 | 2.623469387755102 |
| Trefoil to cinquefoil | 39 | 38 | 4.138775510204082 |

This gives 7,140 retained parameter samples. The historical notebook stopped/cropped at the first row classified as a one-vertex result at every lambda. It is an algorithmic terminal-row rule, not an analytic proof of a physical transition. The HTML rounds its displayed floats to six decimal places; exact candidate values follow the notebook formulas above.

### Hamiltonian interpretation and topology limitations

The model source uses m=2. Most knot/link endpoints use c0=0.5, but the Solomon endpoint uses **c0=0.333**, not exactly 1/3. The interpolation to Solomon therefore changes the endpoint coordinate-map parameter as well as p,q. A manuscript statement that one common c0 is fixed for every endpoint is inaccurate.

The numerical model includes `i Gamma sigma_y`; its interior mask is `spectrum.real == 0`, the exceptional-volume region. For real d1,d3 and positive Gamma this is geometrically `d1^2+d3^2 <= Gamma^2`, the same sublevel set as the Hermitian two-band positive-energy region E <= Gamma. A Hermitian Fermi-volume interpretation can therefore be stated through this equality of geometric sublevel sets, but should not describe the computational Hamiltonian itself as Hermitian.

The default sampled domain is [-pi,pi] x [-pi,pi] x [0,pi], with displayed coordinate scale (1,1,2). The audit found no periodic face-identification operation in this route. The manuscript must specify this finite sampled box and what topology, if any, is intended at its boundary, rather than silently identifying it with a full periodic Brillouin torus.

The historical Hamiltonian notebook also performed substantial phase-label postprocessing: four-connected parameter components smaller than `max(5,ceil(0.01*number_of_displayed_grid_cells))` are iteratively assigned a neighboring large component's label. Votes use the eight neighboring pixels, with total label area and then label index as tie breakers. A cell-area threshold is not evidence that a small topological region is numerical noise. Raw and reassigned maps must be distinguished and the reassigned-cell mask retained.

Retained HTML metadata:

| Transition | Raw classes | Filtered classic classes/regions | Contraction-grouped classes/regions |
|---|---:|---:|---:|
| Hopf to trefoil | 9 | 4 / 4 | 4 / 4 |
| Hopf to Solomon | 28 | 7 / 8 | 4 / 4 |
| Unknot to trefoil | 8 | 5 / 6 | 4 / 5 |
| Unknot to Solomon | 40 | 5 / 5 | 5 / 5 |
| Trefoil to cinquefoil | 19 | 6 / 8 | 5 / 7 |

These counts describe legacy display metadata only. The historical evaluator silently mixed spatial and abstract graph results, and complete raw per-cell arrays were not retained. The old display cannot be reconstructed as a validated Yamada map and is excluded from the corrected manuscript evidence.

## Corrected Hamiltonian pilot and notebook controls

The replacement figure is a new finite-window volume-homology/extraction pilot, **not a Yamada map** and not a corrected rerun of all 7,140 historical cells. It evaluates 1,309 distinct points: 11 interpolation values `lambda_j=j/10` and the same energy-prefix lengths 16,24,16,24,39 listed above, with `E_k=0.30+(99/980)k`. Resolution is 64 cubed, with two explicit worker processes. No parameter-space label filter, periodic face identification, continuum convergence result or embedding certificate is claimed.

The raw sampled volume uses 26-connected foreground and 6-connected complement. Its first Betti number is `b1=b0+b2-chi`. Empty masks, positive-b2 masks and boundary-touching masks are unsupported for the present closed-graph interpretation. Other masks are checked for agreement of the voxel skeleton's and extracted graph's b0/b1 through cleanup, pruning, simplification and accepted smoothing. Homology agreement is necessary but does not prove a deformation retraction or ambient isotopy.

Final outcomes:

| Transition | Betti-consistent extraction | Unsupported boundary |
|---|---:|---:|
| Hopf to trefoil | 83 | 93 |
| Hopf to Solomon | 84 | 180 |
| Unknot to trefoil | 34 | 142 |
| Unknot to Solomon | 43 | 221 |
| Trefoil to cinquefoil | 1 | 428 |
| Total | **245** | **1,064** |

All 1,309 masks have b2=0. Across the finite-window masks, max b0=9 and max b1=7. The first pass took 386.51 seconds on the audit machine. Seven first-pass extraction errors were traced to cleanup deleting already existing isolated components. After preserving those components, the seven affected cells were reevaluated in 2.22 seconds and passed. `first_pass_records.json`, `first_pass_summary.json` and `isolated_component_rerun.json` preserve this history separately from final `records.json` and `summary.json`. The first-pass source snapshot and final corrected notebook hashes are both retained.

The updated notebook06 defaults turn off legacy classification reports, display filtering and abstract-contraction diagnostics. Spatial evaluation failures remain unavailable; high-valence diagram results carry their scope/projection metadata. Empty, failed or topologically inconsistent extractions never become a one-vertex result. A missing validated terminal row retains the full requested energy range. Notebook06 and its normalizer preserve these controls.

The supplied function regression audit now passes **12 checks**, including isolated-component preservation; the repository's notebook-hygiene suite passes **11 tests**. The helper loads only the relevant notebook functions and does not trigger a historical full scan.

Regenerate the final raw pilot and figure from the corrected notebook:

```sh
uv run python dev/run_hamiltonian_homology_pilot.py --repo . \
  --notebook User_guide/applications/06_hamiltonian_yamada_phase_maps.ipynb \
  --out /tmp/new_hamiltonian_homology --seconds 600
uv run python dev/plot_hamiltonian_homology_pilot.py \
  --data /tmp/new_hamiltonian_homology --out /tmp/new_hamiltonian_homology/figure
uv run python dev/audit_hamiltonian_notebook.py \
  --notebook User_guide/applications/06_hamiltonian_yamada_phase_maps.ipynb \
  --out /tmp/new_hamiltonian_function_checks
```

A fresh run uses the corrected cleanup for all points. Figure generation also works directly from the archived final `records.json` and `plan.json`, without rerunning geometry. `figure_arrays.npz` retains each plotted b0, b1, b2 and status array; `figure_metadata.json` records the raw-source hash and absence of filtering. The figure displays finite-window homology even at boundary-touching points and keeps unsupported extraction status separate. The archived PDF's SHA256 is `2dc4e2ca39a67492c6ce8d7d176516f25ca127f591f876a552948ddeb6144db1`.

## Benchmark protocol: known versus absent

The current benchmark notebook and `dev/benchmark_topoly_paper_scaling.py` establish:

- Exactly one timed invariant evaluation for each requested graph/framework pair, in a fresh spawned worker; no explicit warmup or repeated timing distribution.
- Worker `perf_counter` starts just before framework evaluation and stops after conversion to the comparison coefficient dictionary.
- Geometric construction/projection, parent/worker startup and result serialization are outside that worker timing. Separate wall/projection columns exist.
- KnottedGraph: raw invariant, `normalize=False`, `n_jobs=1`, `method="negami"`.
- Topoly: its process-global Yamada memo table is cleared inside the timed function; `point(max_cross=5000)` is used.
- Current notebook limits: KnottedGraph 600 seconds, Topoly 100,000 seconds.
- Projection preprocessing may be reused from a geometry/PD cache. Notebook code recomputes KnottedGraph invariant calls, reuses Topoly results from the ledger, and skips later Topoly family cases after a non-pass.
- The retained comparison ledger has 112 rows. Its columns do not record the historical CPU model, RAM, OS, dependency versions, compiler flags, measured thread-environment values or per-row source commit.

Do not substitute the current audit machine's arm64/macOS environment for missing historical benchmark hardware. Recover the original execution environment if available; otherwise report the known protocol and explicitly identify unavailable metadata. A new benchmark rerun would need a separate dated ledger.

## Integration map

Stage these files under the new repository without changing shared phase/projection code:

| Staged file or directory | Repository destination |
|---|---|
| `audit_discovery_artifacts.py` | `dev/audit_discovery_artifacts.py` |
| `audit_reproducibility_inventory.py` | `dev/audit_reproducibility_inventory.py` |
| `audit_mixed_family.py`, `audit_discovery_long_words.py` | `dev/` |
| `audit_hamiltonian_notebook.py`, `run_hamiltonian_homology_pilot.py`, `plot_hamiltonian_homology_pilot.py` | `dev/` |
| `mixed_family_audit/` | `User_guide/applications/results/arxiv_revision/mixed_family/` |
| `discovery_long_word_audit/` | `User_guide/applications/results/arxiv_revision/discovery_long_words/` |
| `hamiltonian_homology_pilot/` | `User_guide/applications/results/arxiv_revision/hamiltonian_homology/` |
| `hamiltonian_integrity_final_audit/` | `User_guide/applications/results/arxiv_revision/hamiltonian_integrity/` |
| `discovery_audit/` | `User_guide/applications/results/arxiv_revision/discovery/` |
| `reproducibility_inventory/` | `User_guide/applications/results/arxiv_revision/inventory/` |
| this document | `doc/reproducibility/discovery_and_protocols.md` |

The recorded commit has working-tree modifications, listed in provenance. The final release should retain those source hashes and ideally rerun the exact audit at the final immutable release commit. This is source-identification work, not a reason to relabel these retrospective artifacts as historical evidence.
