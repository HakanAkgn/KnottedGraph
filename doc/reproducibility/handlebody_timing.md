# Timing-source audit

The quoted manuscript means and 95th percentiles are numerically correct for the archived timing CSV. The issue is incomplete and unequal cohorts, and the distinction between historical timings and the corrected code.

Source: `User_guide/benchmarks/results/04_thick_handlebody_timing_distributions_source.csv`, introduced by commit `091bfc09a7727d829653e02532be8e1133e28171` on 2026-09-04. SHA256: `13dd29bf4a87632115ea858e9437f672601628bd10ebbe494a5c1b15ffda3516`. It contains 5,000 unique case rows. All six reconstruction/cleanup timing columns match the earlier 5,000-row recovery CSV exactly, including missing entries. They do not come from the new checkout's saved 100-case run.

The repository PDF and manuscript PDF have different byte hashes (`5b21d9...` versus `3f58babc...`) but identical extracted annotations and pixel-identical 110dpi renderings at 1793x1364 pixels. The source CSV exactly reproduces every displayed mean and percentile, resolving the earlier uncertainty about which timing dataset supports this figure.

| Stage | Available n | Mean (s) | P95 (s) |
|---|---:|---:|---:|
| Skeletonization | 5,000 | 0.05227444585 | 0.08565950774 |
| Extraction | 5,000 | 0.008936261936 | 0.01371189166 |
| Contraction | 4,900 | 0.02421401499 | 0.03931169582 |
| Leaf removal | 4,900 | 0.00006572524713 | 0.0001268371852 |
| Simplification | 4,900 | 0.01264999040 | 0.02042661808 |
| Smoothing | 4,900 | 0.01343247545 | 0.02131140651 |
| Projection attempts | 5,000 | 0.03099758730 | 0.06307608990 |
| KnottedGraph, all available | 4,784 | 0.001144902095 | 0.003268609802 |
| KnottedGraph, matched completed subset | 3,491 | 0.0005715492375 | 0.001459166999 |
| Topoly, same matched completed subset | 3,491 | 3.634864190 | 10.79445710 |

The candidate replacement figure uses all available positive finite times in panels a-g and the **same 3,491 cases** for both invariant panels h-i. Both engines must have a duration, Topoly status must be `ok`, Topoly agreement up to a unit must be recorded as true, and recomputed KG agreement with the stored recovered value must be true. The gate does not claim agreement of the seed and recovered embedding; it only defines a fair external-engine timing comparison on the selected diagrams.

There are 1,241 Topoly timeouts, 52 errors, and 216 skipped rows with no selected PD. Five timeout rows retain a stale true match flag but have no duration; a truth flag alone is therefore not an adequate inclusion rule. The strict status gate excludes them. Also, 184 successful Topoly durations exceed their recorded `topoly_timeout_s=10`, up to 449.282 seconds. The CSV alone does not establish a uniform enforced timeout. Do not describe this dataset as uniformly censored at 10 seconds without the run/retry logs.

The ratio of mean engine times on the matched subset is approximately 6,360; the median per-case ratio is approximately 58.3. These answer different questions and must not be substituted for one another. The figure reports durations, not a universal speedup factor.

All timings are historical and were not rerun under the corrected hemisphere sampler/endpoint-incidence code. They should not be presented as current-code performance or as evidence that all 5,000 inputs passed invariant validation. The corrected invariant replay is separate evidence. A current performance claim would require a new controlled timing run with the revised code and matched selected diagrams.

Deliverables: `Time_distributions_verified.pdf`, its unchanged `_source.csv`, exact `_plotted_records.csv` (one panel/case/time per row), `_summary.csv`, `_provenance.json`, and portable `build_timing_figure.py`. No original PDF or source data was overwritten.
