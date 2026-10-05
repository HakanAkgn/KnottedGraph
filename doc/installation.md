# Installation

## Choose the API version first

KnottedGraph requires Python 3.11 or newer. This documentation describes the
0.2.0 development API on `main`. The website is built from that branch.

The [PyPI release](https://pypi.org/project/knotted_graph/) uses the legacy
0.1.2 nodal-only API. Install this source checkout for the examples in these
documents. The software and application papers are listed in {doc}`citing`.

## Recommended source setup with uv

Install [uv](https://docs.astral.sh/uv/), then clone the documented branch:

```bash
git clone --branch main --single-branch \
  https://github.com/HakanAkgn/KnottedGraph.git
cd KnottedGraph
```

For the generic library and the test tools, run:

```bash
uv sync --group dev
```

This command installs the base dependencies and the `dev` dependency group.
To include every optional Python workflow, use:

```bash
uv sync --group dev --all-extras
```

To add only selected workflows, name each extra explicitly. For example:

```bash
uv sync --group dev --extra surface --extra notebook
```

Use `uv run` so commands execute inside the managed environment:

```bash
uv run python examples/quickstart.py
```

The coordinate-curve example is the base-install smoke test. The complete test suite imports
optional workflows; install all extras before running it:

```bash
uv sync --group dev --all-extras
uv run pytest
```

## Source setup with pip

If uv is unavailable, create and activate a Python 3.11-or-newer virtual
environment, clone `main` as above, and install the checkout:

```bash
python -m pip install --upgrade pip
python -m pip install -e .
```

Add one or more extras with standard pip syntax:

```bash
python -m pip install -e ".[surface,notebook]"
```

Use `python -m pip install -e ".[all]"` only when all optional Python workflows
are required.

## Optional-feature matrix

The base installation contains the graph data structures, embedded-graph
utilities, projection and PD-code pipeline, Yamada polynomial backends, and
their required numerical dependencies.

| Install target | Additional capability | Setup notes |
| --- | --- | --- |
| `knot-fields` | Analytic knot-field level sets, tubular diagnostics, and spatial-graph extraction | Field construction is available in the base install; sampled level-set extraction adds scikit-image |
| `nodal` | Non-Hermitian nodal skeleton extraction | Adds `poly2graph`, scikit-image, PyVista, minorminer, and tabulate |
| `surface` | Surface-mesh workflows | Adds PyVista |
| `viz` | Interactive plots and publication-image export | Adds Plotly, Kaleido, PDF/image conversion packages |
| `repulsion` | Python-side repulsive-layout I/O and visualization | Install the native C++ solver separately; see below |
| `notebook` | Interactive notebooks | Adds JupyterLab |
| `all` | Every optional Python workflow | Also adds `igraph`; native system dependencies remain separate |

If you know your starting data or intended result but not the relevant extra,
use the {doc}`feature_status` matrix. It links each documented route to its
public call, return object, next step and computational scaling.

Install the `dev` and `docs` dependency groups with `--group`:

```bash
uv sync --group dev       # test and lint tools
uv sync --group docs      # documentation build tools
```

Groups and extras are independent. For example, a documentation environment
with every optional feature is:

```bash
uv sync --group docs --all-extras
```

## Repulsive-layout native dependency

The `repulsion` extra supplies the Python dependencies. Install the Repulsor
C++ source and its build libraries separately. From a source checkout,
bootstrap the pinned upstream revision with:

```bash
uv sync --extra repulsion
uv run python scripts/bootstrap_repulsion.py --skip-python-install
export REPULSOR_ROOT="$PWD/external/Repulsor"
```

The reference Linux/WSL build requires a C++20 compiler, OpenBLAS,
LAPACK/LAPACKE, `fmt`, and AMD/SuiteSparse. Continue with the
{doc}`user_guide/repulsive_layout` guide before running that backend.

## Verify the installation

Run the small 3D trefoil example to check the coordinate-input, projection and
Yamada workflow:

```bash
uv run python examples/quickstart.py
```

The final two lines should be:

```text
Projection crossings: 3
Normalized Yamada: -Y**11 + Y**9 + Y**8 + Y**7 - Y**4 - Y**3 - Y**2 - Y - 1
```

Continue with the {doc}`quickstart` after these results appear.
If they do not, use the symptom-based checks in {doc}`troubleshooting`.
