# Material phase-map examples

**Start with the [guided walkthrough](../../../doc/applications/material_phase_maps.md).**
It separates reading saved records, making a raw plot, a coarse new scan and
publication reproduction. Run commands from the repository root after installing
the 0.2 development checkout with `uv`.

For the data provenance and interpretation,
see [the provenance notes](INTEGRATION_NOTES.md).

```bash
uv run python -m knotted_graph.applications.phase_map_examples --help
uv run python -m knotted_graph.applications.phase_map_examples inspect \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/data/material_parameter_phase_map_records.csv
```

## Directory map

| Location | Purpose | Edit / regenerate? |
| --- | --- | --- |
| `RealMaterials/data/` | accepted record tables, metadata and audit reports | preserve; write new runs elsewhere |
| `*/figures/`, `*/html/` | original paper plots and interactive artifacts | preserve as reference assets |
| `*/scripts/` | compatibility entries and specialized research presentation tools | portable input/output arguments below |
| `src/knotted_graph/applications/phase_map_examples/` (repository root) | installed application helpers and private compute engines | maintained Python implementation |
| `dev/build_phase_map_demos.py` (repository root) | lossless, lazy-loading website packaging | no scan; generated assets ignored by Git |

The material reference assets came from Hakan's `56bfbab` update, with the
required boundary-topology implementation completed in `643fef8`. Recorded absolute paths
inside original data/HTML record the producing machine. The commands below use
repository-relative paths and explicit output locations. The separate timing source CSV and `assets/paper/Time_distributions.pdf`
from `091bfc0` are also retained without alteration.

## Specialized reproduction commands

Reading the saved CSVs and viewers needs no computation server. New
output belongs under `_build/new_phase_maps/` or an explicit scratch directory,
not beside the reference records. The examples below do not overwrite accepted
files. They require `uv sync --extra nodal --extra viz` except raw record reading.

| Script under this folder | Purpose and key options |
| --- | --- |
| `RealMaterials/scripts/material_parameter_phase_maps.py` | detailed material scan; `--output-dir`, `--dimension`, `--lambda-count`, `--energy-count`, `--only`, `--workers` |
| `RealMaterials/scripts/reproduce_multiband_material_surfaces.py` | original surface-gallery reproduction; `--output-dir`, `--only` |
| `RealMaterials/scripts/generate_material_phase_map_3panel.py` | reproduce the two material panels from saved region HTML; `--html`, `--output-dir` (historical filename retained) |
| `RealMaterials/scripts/integrate_material_maps_into_region_geometry.py` | replace selected material transitions while preserving other base-viewer transitions; `--input-html`, `--data-dir`, `--output`, mandatory `--historical-display-processing` acknowledgement |

For example, render saved material panels without a dense scan:

```bash
uv run --no-sync python \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/scripts/generate_material_phase_map_3panel.py \
  --output-dir _build/new_phase_maps/material_panels
```

## Reference processing settings

To reproduce the saved presentation, match the parameters and energy grids in
the reference summaries, including the recorded incremental extensions. The
detailed material engine classifies cells directly by default. Its optional
processing settings are:

| Setting | Effect on the output |
| --- | --- |
| `--adaptive-energy-step 0.05` | Enables the reference adaptive energy sampling/filling setting. |
| `--min-stable-cells 8 --island-merge-strategy below --apply-signature-merges` | Applies the reference display partition. A displayed group can contain several distinct polynomial signatures; use the raw signatures for polynomial comparisons. |
| `--apply-c6-review` | Applies the recorded C6 display partition. Default: off. |
| `--apply-resolution-calibration` | Applies resolution calibration. Default: off; calibrated records use `resolution_calibration_energy`, a separate status from adaptive fill. |

When reusing records, the selected processing writes to your output directory.
Use a working copy of the input to explore different settings.

The material geometry script uses the dominant resolved volume body and its
outer/nested boundary fillings from the boundary-topology implementation in
`643fef8`. Select `--historical-display-processing` to use its reference display
processing. This replaces the earlier `--primary-component-min-fraction`
convention; the CLI reports a migration error for that earlier option. Keep the
geometry and record classifications together with their version and settings.

Use the guided CLI's `--dry-run` to inspect the scan settings and output paths
before execution.
