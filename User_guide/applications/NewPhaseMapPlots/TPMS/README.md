# Compact-TPMS phase maps: start with saved results

These examples interpolate gyroid, Schwarz-P and diamond fields, then cut a
compact scaffold inside a finite enclosing domain. See the
[complete walkthrough](../../../../doc/applications/material_and_tpms_phase_maps.md).
The resulting boundary topology describes the finite specimen; periodic
boundary identifications are not automatically imposed.

| Directory | What it contains | First use |
| --- | --- | --- |
| [figures](figures/) | Accepted static phase maps | Browse the saved panels |
| [html](html/) | Interactive maps with representative geometry | Download the region viewer and open in a browser with access to its declared Plotly CDN |
| [data](data/) | Records, summaries and existing validation reports | Inspect one family and its domain settings |
| [geometry](geometry/) | Compressed geometry referenced by the saved scan | Keep it together with the selected dataset |
| [scripts](scripts/) | Scan, representative-viewer and validation tools | Use the [script map](../README.md) for explicit paths and options |

Run from the repository root after installing this development checkout:

```bash
uv run python -m knotted_graph.applications.phase_map_examples inspect \
  User_guide/applications/NewPhaseMapPlots/TPMS/data/tpms_parameter_phase_map_records.csv
uv run python -m knotted_graph.applications.phase_map_examples plot \
  User_guide/applications/NewPhaseMapPlots/TPMS/data/tpms_parameter_phase_map_records.csv \
  --family gyroid_to_diamond --output _build/new_phase_maps/my_tpms
```

Each of the three families has 441 reference cells. The raw signature counts
are 65 for `gyroid_to_diamond`, 20 for `gyroid_to_schwarz_p`, and 15 for
`schwarz_p_to_diamond`. PNG/PDF and JSON output give a raw categorical plot and
its full phase key. The reference lacks `classification_computed`, so the
inspector reports its classification status as unrecorded.

The viewer offers `classic`, `stable`, `contraction` and `stable_contraction`
modes. Smoothing and bounded contraction grouping have different meanings from
literal polynomial equality. A region's attachment is one representative,
not every geometry in that region.

For a small new example, install `uv sync --extra nodal --extra viz` on a
compute node, then run:

```bash
uv run --no-sync python -m knotted_graph.applications.phase_map_examples scan tpms \
  --profile quick --output-dir _build/new_phase_maps/my_tpms_scan --dry-run
```

Remove `--dry-run` on the compute node to execute. Expect nine cell records,
summary and plot files, `run_plan.json`, and a `geometry/` directory. Use a fresh
output directory. Inspect boundary contact, domain parameters and geometry
before interpreting the classification; this coarse example is not a
convergence result.

Reference CSV/JSON, HTML and geometry are source assets, not wheel contents.
When moving a dataset, preserve either the scan's `geometry/` subdirectory or
the archived sibling `data/` and `geometry/` layout. A missing geometry error
should be fixed using that dataset's assets, not a file from a different run.
