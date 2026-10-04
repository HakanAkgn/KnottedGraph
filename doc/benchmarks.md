# Benchmarks and saved results

The results below come from the committed CSVs associated with the
[software paper](https://arxiv.org/abs/2609.31152v1). You can download and
inspect them without running a notebook. The full
[data inventory](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/results/README.md)
describes all files; {doc}`paper_results` links them to the published figures.

## Inspect the records first

From a source checkout:

```bash
uv run --no-project python scripts/inspect_paper_data.py
uv run --no-project python scripts/inspect_paper_data.py --json
```

This standard-library command checks eight CSV files against their saved
hashes, columns and row counts, checks unique identifiers, aligns the six
handlebody tables by `case`, and summarizes the recorded outcomes. It checks
archive integrity by reading the saved files; use {doc}`sanity_checks` for
independent small correctness calculations.

The archived records are identified by source snapshot
`9607c1c8a35affdc0dbe90e507c2612b7e464c7d`. Per-run revisions are recorded in
the CSV execution metadata where available.

## Figure 3 and Table 1: fixed-PD invariant evaluation

**Download:** {download}`112-row scaling CSV <../User_guide/benchmarks/results/03_knottedgraph_vs_topoly_scaling_rows.csv>`.
**Method:** [benchmark notebook 03](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/03_knottedgraph_vs_topoly_scaling.ipynb).

The saved cases cover the Dobrynin–Vesnin family and the Li–Lei–Li–Vesnin
cycle/theta edge-replacement families. All 112 KnottedGraph results are marked
`PASS` against the corresponding published formula. The largest cases in
each of the three plotted families contain 500 crossings.

| Family | KnottedGraph at 7 crossings (s) | Topoly at 7 crossings (s) | Recorded ratio | KnottedGraph at 500 crossings (s) |
| --- | ---: | ---: | ---: | ---: |
| Dobrynin–Vesnin Theta(n) | 0.012525 | 0.475815 | 37.9905 | 0.117355 |
| Cycle edge replacement | 0.015948 | 3.188000 | 199.9028 | 2.376074 |
| Theta edge replacement | 0.016153 | 15.439276 | 955.8345 | 4.998060 |

The seven-crossing rows have `PASS` results for both implementations and use
the same PD input. Rounded to the precision used in the paper, their ratios
are approximately 38, 200 and 960. The ratio is `topoly_s / knottedgraph_s`.
Invariant timings exclude the separately recorded `projection_s`.
The table reports runtime for these three structured families under the
recorded environment. Performance for a new graph depends on its diagram and
connectivity.

### Completion status matters

Of 112 Topoly rows, 38 calls completed (`topoly_status=ok`): 37 match the
reference and one has a `FAIL` comparison. A further 46 rows record execution
errors. The remaining 28 were explicitly skipped after an earlier family
error/non-pass (13 and 15 respectively). They have no new Topoly runtime.
The three 500-crossing rows contain completed KnottedGraph evaluations and
non-completed Topoly statuses. Use completed paired rows with matching reference
results for runtime ratios; skipped rows have no recorded runtime. Read
`topoly_status` together with `topoly_result`.

`case_key` uniquely identifies a row. `pd_hash`, `embedding_hash`,
`projection_cache_key` and `settings_hash` retain input/protocol identifiers;
the saved projection caches are under `User_guide/benchmarks/results/03_pdcode_cache/`.

## Supplementary Figure 6: volume-to-graph recovery

**Download:** {download}`4,400-case preservation CSV <../User_guide/benchmarks/results/handlebody_ground_truth/synthetic_ground_truth_yamada_preservation.csv>`,
{download}`input manifest <../User_guide/benchmarks/results/handlebody_ground_truth/certified_thick_handlebody_embeddings_manifest.csv>`,
{download}`recovered graph manifest <../User_guide/benchmarks/results/handlebody_ground_truth/synthetic_ground_truth_recovered_graphs_manifest.csv>`.
**Method:** [benchmark notebook 04](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/04_thick_handlebody_validation.ipynb).

The inverse benchmark starts with a known embedded graph, constructs and
voxelizes a thick neighborhood, validates the sampled input, and recovers a
graph for comparison. All 4,400 saved cases pass the input-topology,
primary-topology and normalized Yamada comparisons, and have
`overall_pass=True`. The table contains the 4,400 inputs accepted by the
sampled-topology checks. Candidate generation precedes acceptance. Compare
sampling resolutions separately to assess stability.

The benchmark measures topology and invariant preservation. Abstract graph
isomorphism is a separate diagnostic, with 456 matching and 3,944 nonmatching
recovered graphs. The topology and invariant checks and the isomorphism
diagnostic measure different properties of recovery.

The CSVs are enough to inspect these counts. Exact input geometry is optional
and uses Git LFS; the
[data inventory](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/results/README.md#geometry-and-cached-diagrams)
provides the retrieval command and links to the saved recovered-graph archive.

## Supplementary Figure 7: pipeline runtime distributions

**Download:** {download}`timing-figure source, 4,400 rows <../User_guide/benchmarks/results/04_thick_handlebody_time_distribution_plot_source.csv>`,
{download}`full stage timings <../User_guide/benchmarks/results/04_thick_handlebody_pipeline_timings.csv>`,
{download}`21-stage summary <../User_guide/benchmarks/results/04_thick_handlebody_pipeline_timing_summary.csv>`.

| Figure panel | CSV column in the timing-figure source |
| --- | --- |
| a: Skeletonization | `skeletonization_seconds` |
| b: Graph extraction | `graph_extraction_seconds` |
| c: Short-edge contraction | `contract_short_edges_seconds` |
| d: Leaf removal | `remove_leaf_nodes_seconds` |
| e: Degree-two simplification | `simplify_edges_seconds` |
| f: Geometric smoothing | `smooth_edges_seconds` |
| g: Projection selection | `projection_seconds` |
| h: KnottedGraph Yamada evaluation | `knottedgraph_yamada_seconds` |
| i: Topoly Yamada evaluation | `topoly_yamada_seconds`, only where `topoly_status=ok` |

All 4,400 cases have the KnottedGraph projection/invariant records. Topoly
completed 3,968 cases; 411 timed out under the 10-second limit and 21 errored.
Only completed evaluations enter its runtime distribution. The other statuses
remain in the downloadable data. The histograms summarize each
implementation's completed cases; use paired completed rows with matching
polynomial results to calculate per-case speedups.

The full stage-timing file also records `git_commit`, `python_version`,
`platform` and `processor`, along with geometry and acceptance metadata.
The summary covers 21 recorded stages, including aggregates of constituent
operations. Use one level of aggregation when adding stage timings.
Candidate search/generation is outside the measured recovery workload.

## Where the remaining evidence lives

- {doc}`sanity_checks`: independent algebraic checks and projection-invariance procedures.
- {doc}`paper_results`: all main/supplementary figure originals and their method links.
- {doc}`applications/yamada_formula_discovery`: saved Figure 4 / Supplementary Figure 10 data, exact coefficients and formula-discovery procedures.
- {doc}`applications/material_phase_maps`: existing material data and viewers, with raw classifications distinguished from adaptive fills and display processing.
