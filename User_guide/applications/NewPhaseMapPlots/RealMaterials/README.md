# Material phase maps: start with saved results

These examples scan a Hamiltonian coefficient with lambda and an energy/gap
threshold for TiB2 and Co2MnGa. Start with the
[complete walkthrough](../../../../doc/applications/material_and_tpms_phase_maps.md).
Commands below run from the repository root after installing this development
checkout.

| Directory | What it contains | First use |
| --- | --- | --- |
| [figures](figures/) | Accepted static material panels | Browse the saved figures |
| [html](html/) | Saved interactive phase maps and representative geometry | Download the region-geometry HTML and open in a browser; Plotly requires its declared CDN |
| [data](data/) | CSV/JSON cell records and scan summaries | Inspect provenance and sampling parameters |
| [scripts](scripts/) | Compatibility entries and specialized reproduction tools | Use explicit input/output paths from the [script map](../README.md) |

Inspect both families without recomputing:

```bash
uv run python -m knotted_graph.applications.phase_map_examples inspect \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/data/material_parameter_phase_map_records.csv
```

The reference contains 30,180 TiB2 records (35 raw signatures) and 15,840 Co2MnGa
records (53 raw signatures). It contains 14,634 directly classified, 31,184
adaptive-fill and 202 resolution-calibration records. The 88 viewer region
attachments include both display modes and the retained nodal examples; they
are not 88 distinct material polynomials.

Choose one family and replot the raw signatures:

```bash
uv run python -m knotted_graph.applications.phase_map_examples plot \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/data/material_parameter_phase_map_records.csv \
  --family tib2_d6_F --output _build/new_phase_maps/my_tib2
```

Success produces PNG, PDF and a JSON phase key. The JSON preserves full
signatures, source labels and cell classification/calibration status. Category
IDs are labels, not polynomial values. `source=yamada-set` can include outer
and inner boundary results and structural fallbacks; inspect its contents.

For new computation, install `uv sync --extra nodal --extra viz` on a compute
node, then inspect the plan:

```bash
uv run --no-sync python -m knotted_graph.applications.phase_map_examples scan materials \
  --profile quick --output-dir _build/new_phase_maps/my_material_scan --dry-run
```

Remove `--dry-run` on the compute node to execute. The default writes nine cells
from a coarse scan plus summaries, plots and `run_plan.json`. Use a new or empty
directory. This exercise does not demonstrate high-resolution convergence.

The guided scan classifies each requested cell without adaptive filling,
resolution calibration or display merging. Detailed research options include
`--apply-signature-merges`, `--apply-c6-review` and
`--apply-resolution-calibration`; their scientific/presentation roles are
explained in the walkthrough. Reuse or extend a working copy of a dataset,
keeping its records and summary together. Reference assets remain unchanged.

If imports fail, verify this checkout includes the boundary helpers from
upstream `643fef8` and that the optional dependencies are installed. If a family
is ambiguous, select the key printed by `inspect`. For more failures, see
[Troubleshooting](../../../../doc/troubleshooting.md).
