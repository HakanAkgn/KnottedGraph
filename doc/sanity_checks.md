# Sanity checks

Use the small checks below before interpreting a new invariant calculation.
To inspect the paper's existing results without running any calculations, use
{doc}`benchmarks` and its saved-record checker instead.

## Independent algebraic and public-API checks

From the repository root, install the base development environment and run:

```bash
uv sync --group dev
uv run python dev/run_yamada_sanity_checks.py
```

The [standalone script](https://github.com/HakanAkgn/KnottedGraph/blob/main/dev/run_yamada_sanity_checks.py)
checks trees, cycles, bouquets, theta graphs, an isthmus, a one-point union
and planar K4 against explicit formulas. It also compares native-dispatched
results with the exact Python compact evaluator, checks the two public
methods on planar and one-crossing theta embeddings, and checks the mirror
relation. These small correctness calculations use the base package and
test tools.

Success ends with:

```text
PASS: all published/independent Yamada sanity checks succeeded.
```

If an assertion fails, retain its full message, source commit and environment
details. These checks compare exact symbolic expressions.

## Projection invariance on fixed spatial embeddings

Open [benchmark notebook 01](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/01_yamada_sanity_checks.ipynb):

```bash
uv sync --group dev --extra notebook
uv run jupyter lab User_guide/benchmarks/01_yamada_sanity_checks.ipynb
```

Its first stage runs the standalone checks above. The next stage tests 25
deterministic connected trivalent embeddings: five each of theta, K4,
triangular-prism, K3,3 and cube graphs. It holds each embedding fixed, samples
projection directions, checks the returned PD code and crossing count, and
requires the same normalized Laurent polynomial across all valid views.
The final message reports projection invariance for all 25 embeddings.
Each family contributes five spatial embeddings of the same abstract graph.

## Further checks and their scope

| Check | Entry point | Coverage |
| --- | --- | --- |
| Five-minute installation check | `uv run python examples/quickstart.py` | 3D trefoil input, three-crossing projection and expected normalized Yamada polynomial |
| Small CSV input check | `uv run python examples/input_to_yamada.py` | Planar theta graph loaded from node/edge files agrees with its explicit crossing-free formula |
| Saved paper-data integrity | `uv run --no-project python scripts/inspect_paper_data.py` | File hashes, CSV columns, case alignment and saved statuses; no scientific calculation |
| Application regressions | [Notebook 02](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/02_application_regression_checks.ipynb) | Current vs historical application outputs; requires the reference Git revision and application dependencies |
| Published structured families | [Notebook 03](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/03_knottedgraph_vs_topoly_scaling.ipynb) | Published formulas and fixed-PD timing protocol; see {doc}`benchmarks` before executing |
| Volume-to-graph recovery | [Notebook 04](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/04_thick_handlebody_validation.ipynb) | Topology/invariant preservation for certified sampled inputs; publication-scale workload |

These checks cover the listed algebraic identities and projection examples.
The {doc}`paper_results` guide links each scientific dataset to its validation
procedure.
