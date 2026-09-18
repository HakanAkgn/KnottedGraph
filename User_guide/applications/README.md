# Choose an application

Complete the [Quick Start](../../doc/quickstart.md) first. These examples apply
the graph/projection/Yamada tools to specific scientific models; each notebook
states its inputs, dependencies and interpretation limits.

| Notebook | Use it for | Execution boundary |
| --- | --- | --- |
| [Physics](01_physics_applications.ipynb) | Nodal and material skeleton examples | Includes large 3-D sampling examples; review resource notes before running |
| [Mathematics](02_mathematics_applications.ipynb) | Mathematical graph families and invariants | Early graph examples are small; later galleries sample 200³/300³ grids and write results |
| [Proteins](03_protein_applications.ipynb) | Inspect an ordered protein backbone and its input contract | Base installation; creates a temporary offline PDB example; does not infer a protein interaction graph |
| [Analytic knot fields](04_analytic_knot_fields.ipynb) | Knot functions, level sets and deformations | Needs the extraction stack; the saved configuration enables paper mode and extensive scans |
| [Formula discovery](05_yamada_formula_discovery.ipynb) | Research searches for graph-family formulas | Large searches are reproduction workloads; preserve the stated branch restrictions |
| [Hamiltonian/Yamada phase maps](06_hamiltonian_yamada_phase_maps.ipynb) | Unified finite-grid phase-map workflows | Scan cost depends on resolution and graph/projection complexity |

For a first material or compact-TPMS exercise, use the
[saved-result walkthrough](../../doc/applications/material_and_tpms_phase_maps.md)
and [research directory map](NewPhaseMapPlots/README.md). Reading or replotting
saved records uses the base installation; a small new scan needs `nodal` and
`viz`. The source checkout supplies the reference assets, which are not bundled
in the wheel.

Follow the sequence **browse → inspect records → replot → small scan → research
reproduction**. A colored region is a recorded classification/display group,
not a proof that every represented graph is homeomorphic. Read the polynomial,
structural fallback, error, calibration and grouping information separately.

Execute and render on a compute node when using a cluster. Preserve notebook
sources and accepted results, and keep execution copies in a separate output
directory. Execution success and visual review should be recorded separately.

For a reduced analytic-fields run, make an execution copy and set the existing
`PAPER_FIGURE_MODE=False` and `SAVE_PAPER_FIGURES=False` controls **before**
running its setup cell. This selects the existing `FAST_INTERACTIVE` branch;
`QUICK=True` alone does not disable paper mode. Its optional scan sections then
print their skip messages. A successful reduced run does not validate those
scans or guarantee a certified knot at the reduced grid resolution.

The automatic notebook CI executes the introductory notebooks, the offline
protein example and the two small benchmark notebooks. The mathematics,
analytic-fields and other research notebooks receive source/portability checks;
their complete research workloads require a separate compute run.

From a checkout, preserve execution output with:

```bash
MPLBACKEND=module://matplotlib_inline.backend_inline \
uv run --no-sync python dev/execute_notebook.py \
  User_guide/applications/03_protein_applications.ipynb \
  --output-dir /path/to/fresh/notebook-review
```

Install the `docs` dependency group for HTML export (`uv sync --group docs`),
along with the extras needed by the selected notebook. The runner saves a
notebook, HTML and a JSON execution status even after a cell failure; inspect
the rendered HTML separately. It refuses to overwrite an existing export.
