# Troubleshooting

Use the checks below to diagnose installation problems, validate embedded
geometry and inspect the computational cost of topology calculations before
changing scientific parameters.

## Confirm the interpreter and package version

KnottedGraph 0.2 requires Python 3.11 or newer. In the recommended uv checkout,
inspect the managed environment without activating it:

```bash
uv run python --version
uv run python -c "import knotted_graph; print(knotted_graph.__version__); print(knotted_graph.__file__)"
uv pip check
```

In an activated conventional venv, use `python --version`, `python -m pip
--version`, the same `python -c ...` command, and `python -m pip check`. Run
all checks with the interpreter from your selected environment.

The current 0.2 development API is newer than the legacy 0.1 release on PyPI.
If the version command prints `0.1.x`, or imports such as
`knotted_graph.core` and `knotted_graph.inputs` are missing, install the
current main-branch source version using the [installation guide](installation.md).
Use an environment containing the selected source version and its dependencies.

## Install the optional dependency for the workflow

The base install contains the graph, projection, and invariant code. Optional
imports need the corresponding extra:

| Symptom | Install |
| --- | --- |
| `No module named 'pyvista'` while loading a surface | `knotted_graph[surface]` |
| Missing Plotly, Kaleido, Pillow, or PyMuPDF | `knotted_graph[viz]` |
| Missing `poly2graph` or scikit-image in a nodal workflow | `knotted_graph[nodal]` |
| Missing JupyterLab | `knotted_graph[notebook]` |
| Missing Biopython or Plotly in repulsive layout | `knotted_graph[repulsion]` |

For a source checkout, replace `pip install` with the matching `uv sync
--extra ...` command. The `all` extra installs every Python optional dependency.
For native Repulsor, install the C++ source and system libraries with the
[Repulsive Layout guide](user_guide/repulsive_layout.md).

## Distinguish a source checkout from an installed package

The beginner notebooks accept either:

1. a checkout containing `src/knotted_graph`; or
2. an environment in which the current package is installed.

If a notebook reports that KnottedGraph cannot be found, launch Jupyter through
the environment used to install the package:

```bash
uv run jupyter lab
```

Avoid adding a guessed `/src` directory to `sys.path`. Print
`knotted_graph.__file__` when you are unsure which copy Python imported.

## Validate an embedded graph before projection

Projection expects an undirected `networkx.MultiGraph` with finite three-
dimensional node positions. Each node stores `pos` with shape `(3,)`; each edge
stores a polyline `pts` with shape `(N, 3)` whose endpoints match its incident
nodes. Missing edge polylines can be materialized as straight segments by
`ensure_embedding`.

```python
from knotted_graph.core import ensure_embedding, validate_embedding

issues = validate_embedding(graph)
if issues:
    raise ValueError("; ".join(issues))

graph = ensure_embedding(graph)
```

Typical failures identify the exact node or edge. Fix missing coordinates,
non-finite values, endpoint mismatches, empty graphs, and directed graphs before
sampling projections.

## Fix projection failures

Rotation angles are three values in **degrees** by default. A rotation order is
a three-character Euler sequence such as uppercase `ZYX` (extrinsic) or
lowercase `xyz` (intrinsic). Keep the sequence entirely uppercase or lowercase.
To reproduce one view, pass its angles explicitly. To search for a cleaner
regular projection, increase the sample count gradually:

```python
from knotted_graph.projection import sample_projections, select_projection

projection = select_projection(graph, num_rotation_samples=10)
print(projection.rotation_angles, projection.num_crossings)

alternatives = sample_projections(graph, num_rotation_samples=20)
```

An overlapping-collinear-segment error identifies a degenerate view. Try
another rotation to obtain a regular diagram. If some sampled views fail,
KnottedGraph reports a warning and keeps the valid views; if every view fails,
the raised error includes the failed sample details.

## Understand slow Yamada calculations

A diagram with `c` crossings has `3**c` resolved states before memoization and
structural shortcuts. The warning emitted at ten crossings indicates
computational cost. Inspect the selected crossing count
before computing and start with one worker:

```python
import sympy as sp
from knotted_graph.projection import compute_yamada_polynomial

Y = sp.Symbol("Y")
result = compute_yamada_polynomial(
    graph,
    Y,
    num_rotation_samples=10,
    n_jobs=1,
    return_result=True,
)
print(result.projection.num_crossings)
print(result.polynomial)
```

Use `n_jobs=-1` only when you have intentionally allocated all available CPU
cores. A result of zero can be mathematically correct: in particular, the
Yamada polynomial vanishes for a graph containing a bridge (cut edge).

## Run surface workflows without a display

On a server, set PyVista to off-screen mode before rendering:

```bash
export PYVISTA_OFF_SCREEN=true
```

Surface loading returns a `PolyData` mesh and reports open boundary edges in
`result.issues`. Choose a skeleton-extraction route using the
[Workflow Overview](user_guide/workflow_overview.md). If PyVista is missing,
install the `surface` extra.

## Diagnose Repulsor setup

With uv, install the `repulsion` Python dependencies first, then bootstrap
the native source with `--skip-python-install`:

```bash
uv sync --extra repulsion
uv run python scripts/bootstrap_repulsion.py --skip-python-install
export REPULSOR_ROOT="$PWD/external/Repulsor"
uv run kg-repulsive-layout --help
```

Inside an activated conventional venv with pip, the bootstrap may install the
extra itself: `python scripts/bootstrap_repulsion.py`.

Compiler, BLAS/LAPACK, `fmt`, or SuiteSparse errors are native dependency
problems; see the [Repulsive Layout guide](user_guide/repulsive_layout.md) for
the required libraries and reproducibility checks.

## Report a reproducible issue

Include the following when opening a [GitHub issue](https://github.com/HakanAkgn/KnottedGraph/issues):

- the output of the environment/version commands at the top of this page;
- the smallest input or graph that reproduces the problem;
- `validate_embedding(graph)` or adapter `result.issues`;
- the explicit rotation angles, order, and crossing count for projection issues;
- the selected optional extras and operating system; and
- the complete traceback as text.

Remove private paths or unpublished scientific data before posting.
