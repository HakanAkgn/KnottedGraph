# Material and compact-TPMS phase-map examples

**Start with the [guided walkthrough](../../../doc/applications/material_and_tpms_phase_maps.md).**
It separates reading saved records, making a raw plot, a coarse new scan and
publication reproduction. Run commands from the repository root after installing
this development branch with `uv`; the old PyPI API is not sufficient.

For the integration scope, validation results and scientific feedback requested,
see [the review notes](INTEGRATION_NOTES.md).

```bash
uv run python -m knotted_graph.applications.phase_map_examples --help
uv run python -m knotted_graph.applications.phase_map_examples inspect \
  User_guide/applications/NewPhaseMapPlots/TPMS/data/tpms_parameter_phase_map_records.csv
```

## Directory map

| Location | Purpose | Edit / regenerate? |
| --- | --- | --- |
| `RealMaterials/data/`, `TPMS/data/` | accepted record tables, metadata and audit reports | preserve; write new runs elsewhere |
| `TPMS/geometry/` | per-cell compressed geometry for the accepted scan | preserve |
| `*/figures/`, `*/html/` | original paper plots and interactive artifacts | preserve pending scientific review |
| `*/scripts/` | compatibility entries and specialized research presentation tools | portable input/output arguments below |
| `src/knotted_graph/applications/phase_map_examples/` (repository root) | installed application helpers and private compute engines | maintained Python implementation |
| `dev/build_phase_map_demos.py` (repository root) | lossless, lazy-loading website packaging | no scan; generated assets ignored by Git |

The material reference assets came from Hakan's `56bfbab` update, with the
required boundary-topology implementation completed in `643fef8`. The TPMS
assets remain those from `2b2ae6d`. Recorded absolute paths
inside original data/HTML describe its provenance; they are not required on your
machine. The separate timing source CSV and `assets/paper/Time_distributions.pdf`
from `091bfc0` are also retained without alteration.

## Specialized reproduction commands

Use a compute node for rendering, geometry regeneration and validation. New
output belongs under `_build/new_phase_maps/` or an explicit scratch directory,
not beside the reference records. The examples below do not overwrite accepted
files. They require `uv sync --extra nodal --extra viz` except raw record reading.

| Script under this folder | Purpose and key options |
| --- | --- |
| `RealMaterials/scripts/material_parameter_phase_maps.py` | detailed material scan; `--output-dir`, `--dimension`, `--lambda-count`, `--energy-count`, `--only`, `--workers` |
| `RealMaterials/scripts/reproduce_multiband_material_surfaces.py` | original surface-gallery reproduction; `--output-dir`, `--only` |
| `RealMaterials/scripts/generate_material_phase_map_3panel.py` | reproduce the two material panels from saved region HTML; `--html`, `--output-dir` (historical filename retained) |
| `RealMaterials/scripts/integrate_material_maps_into_region_geometry.py` | replace selected material transitions while preserving other base-viewer transitions; `--input-html`, `--data-dir`, `--output`, mandatory `--historical-display-processing` acknowledgement |
| `TPMS/scripts/tpms_compact_scaffold_phase_maps.py` | detailed compact-scaffold scan; `--output-dir`, `--dimension`, `--threshold-count`, `--threshold-min`, `--threshold-max`, `--domain`, `--only` |
| `TPMS/scripts/tpms_plotly_region_geometry.py` | build region geometry and bounded contraction groups from records; `--scan-dir`, `--output`, `--max-contraction-states`, `--max-contraction-depth`, `--max-contraction-edges` |
| `TPMS/scripts/make_tpms_compact_1row3col_figure.py` | render the specialized publication layout from saved HTML; `--html`, `--output`, optional `--include-panel-c-slice`, `--geometry-dir` |
| `TPMS/scripts/verify_tpms_compact_phase_maps.py` | audit saved TPMS outputs; `--scan-dir`, `--html`, `--output`; nonzero exit status on validation failures |

For example, render saved material panels without a dense scan:

```bash
uv run --no-sync python \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/scripts/generate_material_phase_map_3panel.py \
  --output-dir _build/new_phase_maps/material_panels
```

The detailed material engine defaults to no adaptive fill or display smoothing.
For historical display processing, the relevant switches include
`--adaptive-energy-step 0.05 --min-stable-cells 8 --island-merge-strategy below
--apply-signature-merges`. This is an explicit presentation choice, not a claim
that the merged polynomials are equal. Match **all** parameters and saved energy
grids from the reference summaries before claiming full reproduction; the
historical scan also included incremental extensions. The new research switches
`--apply-c6-review` and `--apply-resolution-calibration` separately enable the
teacher's C6 display partition and resolution calibration; both are off by
default. Calibration is recorded with `resolution_calibration_energy` and is
distinct from ordinary adaptive fill. Reusing records applies any explicitly
selected processing and writes to the chosen output directory, so use a copy of
accepted input when exploring these choices.

The material geometry script now uses the dominant resolved volume body and
its outer/nested boundary fillings, matching the latest teacher implementation.
It requires an acknowledgement because it also performs historical display
processing. The former `--primary-component-min-fraction` convention no longer
applies; supplying that option produces a migration error. The geometry and
record classifications must be interpreted with the same version and settings.

The guided CLI's `--dry-run` is the safest way to inspect a new run plan. None of
these tools submit jobs or select a billing account on your behalf.
