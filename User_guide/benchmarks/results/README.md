# Saved benchmark data

These are existing records from source snapshot
`9607c1c8a35affdc0dbe90e507c2612b7e464c7d`. No benchmark was rerun to build this
index. That snapshot identifies the archived files; the original execution
revision and environment are recorded in the timing rows where available.

Read the [benchmark interpretation guide](../../../doc/benchmarks.md) alongside
the CSVs. The [paper guide](../../../doc/paper_results.md) maps figures to records.

| Saved file | Rows | What it supports |
| --- | ---: | --- |
| [Structured-family scaling](03_knottedgraph_vs_topoly_scaling_rows.csv) | 112 | Figure 3 and Table 1: published-polynomial checks and fixed-PD timings, through 500 crossings |
| [Handlebody topology and Yamada preservation](handlebody_ground_truth/synthetic_ground_truth_yamada_preservation.csv) | 4,400 | Supplementary Figure 6: input acceptance and recovery diagnostics |
| [Stage timings](04_thick_handlebody_pipeline_timings.csv) | 4,400 | Per-case pipeline runtimes and execution metadata |
| [Stage timing summaries](04_thick_handlebody_pipeline_timing_summary.csv) | 21 | Mean, median, standard deviation, quartiles, 95th percentile and range |
| [Timing-figure source](04_thick_handlebody_time_distribution_plot_source.csv) | 4,400 | Supplementary Figure 7: selected stage timings and Topoly completion status |
| [Certified input manifest](handlebody_ground_truth/certified_thick_handlebody_embeddings_manifest.csv) | 4,400 | Seeds, geometry and input-topology acceptance |
| [Recovered graph manifest](handlebody_ground_truth/synthetic_ground_truth_recovered_graphs_manifest.csv) | 4,400 | Exact graph hashes, polynomial strings and comparisons |
| [Recovery checkpoint](handlebody_ground_truth/synthetic_ground_truth_results_checkpoint.csv) | 4,400 | Saved recovery ledger |

## Check the archive without recomputing

```bash
uv run --no-project python scripts/inspect_paper_data.py
uv run --no-project python scripts/inspect_paper_data.py --json
```

The command uses only Python's standard library. It verifies SHA-256 hashes,
columns and row counts from [manifest.json](manifest.json), checks unique case
identifiers and joins the six handlebody tables by `case`. It then summarizes
the **saved** pass/fail and timing-status fields. A changed or missing file
produces a nonzero exit status. This is an archive-integrity check, not a new
independent evaluation of the polynomials.

## Key fields and comparison rules

- Scaling: `case_key` is the unique case identifier; `crossings` describes the
  fixed PD input. `knottedgraph_s` and `topoly_s` time invariant evaluation,
  while `projection_s` is separate. Compare timings only when both results are
  `PASS` and both execution statuses are `ok`.
- `topoly_result` is the comparison outcome; `topoly_status` says whether the
  call completed, errored or was skipped. A row marked
  `skipped_after_family_error` or `skipped_after_family_nonpass` has no new
  Topoly timing. Its `ERROR` result label must not be counted as a new
  execution failure or assigned a fabricated runtime.
- Handlebody: join by `case`, not row position. `input_topology_pass`,
  `primary_topology_pass`, `yamada_match` and `overall_pass` test different
  stages. `abstract_isomorphic` is a separate diagnostic: it is true in 456
  cases, while all 4,400 saved Yamada comparisons pass.
- Timing-figure source: `topoly_status=ok` selects the 3,968 completed Topoly
  cases. The 411 timeouts and 21 errors remain recorded separately. The
  completed-time distribution is censored by a 10-second limit; the marginal
  distributions do not provide a paired speedup over all 4,400 cases.
- Environment: `git_commit`, `python_version`, `platform` and `processor` in
  the stage-timing table describe the original run. Empty metadata fields
  mean the value was not recorded, not that all machines are equivalent.

## Geometry and cached diagrams

The recovered graph archive is
[synthetic_ground_truth_recovered_graphs.json.gz](handlebody_ground_truth/synthetic_ground_truth_recovered_graphs.json.gz).
The 108,611,015-byte certified input archive is stored with Git LFS. To obtain
that optional geometry from a clone with Git LFS installed:

```bash
git lfs pull --include="User_guide/benchmarks/results/handlebody_ground_truth/certified_thick_handlebody_embeddings.json.gz"
```

All CSVs can be inspected without it. The `03_pdcode_cache/` directory contains
the saved projection caches referenced by the scaling rows; these Python
pickle caches are method intermediates, whereas CSVs are the portable first
entry point for inspecting reported results.

## Formula-discovery data availability

The [formula-discovery notebook](../../applications/05_yamada_formula_discovery.ipynb)
contains the constructors and procedures for Figure 4 and Supplementary
Figure 10. Its generated `05a_*`, `05b_*` and `05c_*` CSV/JSON records are **not
included in this snapshot**, and the notebook contains no saved execution
outputs. Their intended location is
`User_guide/applications/results/05_yamada_formula_discovery/`.
The separate 87-row structured-graph demonstration CSV is not a substitute
for these held-out validation records. The paper guide labels this gap
explicitly; no replacement data have been generated.
