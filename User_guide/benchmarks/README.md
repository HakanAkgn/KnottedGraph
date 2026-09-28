# Benchmarks and correctness checks

Start with the [published results and CSV guide](../../doc/benchmarks.md),
or the [sanity-check instructions](../../doc/sanity_checks.md) for a small local check.
All commands below run from the repository root.

| Entry | Purpose | Does it compute new results? |
| --- | --- | --- |
| [Saved data inventory](results/README.md) | Download CSVs, identify fields and inspect recorded results | No |
| `uv run --no-project python scripts/inspect_paper_data.py` | Check hashes, case alignment and recorded status counts | No |
| `uv run python dev/run_yamada_sanity_checks.py` | Independent algebraic identities and small public-API checks | Yes, small correctness checks |
| [01: Yamada sanity checks](01_yamada_sanity_checks.ipynb) | The script above plus 25 deterministic trivalent embeddings under multiple projections | Yes |
| [02: Application regression checks](02_application_regression_checks.ipynb) | Compare current application output with a historical revision in temporary worktrees | Yes; needs Git history and application dependencies |
| [03: KnottedGraph vs Topoly](03_knottedgraph_vs_topoly_scaling.ipynb) | Structured-family invariant-evaluation benchmark | Yes; inspect the saved CSV first |
| [04: Thick handlebody validation](04_thick_handlebody_validation.ipynb) | Certified volume-to-graph construction and topology preservation | Yes; publication workload includes 4,400 volumes |

The saved results are already in this checkout. Do not use **Run All** merely to
view a paper figure. The [figure guide](../../doc/paper_results.md) includes the
published PDF originals. Benchmark notebooks can update their result ledgers;
use a separate checkout/output location for a new experiment.
