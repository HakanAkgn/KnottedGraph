# Material Phase Maps

Start here for the material-parameter examples supplied in
`User_guide/applications/NewPhaseMapPlots`. You can inspect saved results without
running a scan. The application provides material-specific models and helpers
to inspect, plot and compute record tables. The [phase-map API](../api/applications.md)
documents the general scan interface.

## Choose how much work to run

| Goal | Entry | What runs / prerequisites |
| --- | --- | --- |
| Understand the published examples | previews and interactive links below | browser; Plotly needs access to its pinned CDN |
| Inspect or replot saved cells | `phase_map_examples inspect` / `plot` | base installation; reads CSV/JSON, no extraction |
| Learn the computation | `phase_map_examples scan --profile quick` | optional `nodal` and `viz`; small optional calculation |
| Reproduce high-resolution research plots | source scripts and reference data | advanced; expensive scans, explicit display processing and validation |

Install the `main` source version with the
[source-install instructions](../installation.md) to use these commands. They
assume the repository root; the installed Python helpers also work outside a
checkout when supplied with your own records file.

## Saved research results

The material examples interpolate one Hamiltonian coefficient with $\lambda$
and vary an energy/gap threshold $E$.

<div class="kg-wide-figure">
  <img src="../demos/new_phase_maps/materials/preview.png" alt="Saved material-parameter phase-map panels for TiB2 and Co2MnGa">
</div>

<p><a href="../demos/new_phase_maps/materials/index.html">Open material and Hamiltonian regions</a>.</p>
The viewer includes nodal transitions and material scans. Select TiB2 or
Co2MnGa to explore the material results.

Choose a transition and classification mode, then click a region to see its
representative surface and skeleton from a selected cell in that region.
In “Up to contraction” views, colors follow the reported grouping convention.

The viewer loads geometry **on selection**. If Plotly or a region download
is blocked, use the static figures. Serve downloaded split viewers over HTTP.
The original self-contained research HTML files remain in the source tree;
use their declared browser dependencies.

## First exercise: inspect saved cells

Inspect saved records with the base installation:

```bash
uv run python -m knotted_graph.applications.phase_map_examples inspect \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/data/material_parameter_phase_map_records.csv
```

Expected reference counts:

| Family key | Records | Grid | Distinct raw signatures |
| --- | ---: | --- | ---: |
| `tib2_d6_F` | 30,180 | 60 lambda × 503 energy samples | 35 |
| `co2mnga_t8` | 15,840 | 60 lambda × 264 energy samples | 53 |

The inspector also reports errors, missing cells, source counts, direct
classifications, adaptive fills and resolution calibrations.

## Replot without recomputing

```bash
uv run python -m knotted_graph.applications.phase_map_examples plot \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/data/material_parameter_phase_map_records.csv \
  --family tib2_d6_F --output _build/new_phase_maps/tib2_raw
```

This saves `tib2_raw.png`, `tib2_raw.pdf` and `tib2_raw.json`. The JSON is the
phase key: it maps each category ID to its full signature, color, source and
recorded polynomial, and includes the coordinate axes and the category-ID grid.
White cells are missing, shaded cells are adaptive fills, open circles mark
resolution calibrations, and crosses mark errors. The JSON also records each
cell's classification status and calibration anchor energy when present.
This raw view retains the recorded signatures, with small-island smoothing,
C6 grouping and manual signature merging disabled. Use the JSON key to identify
categories when the plot contains many colors.

The same operation from Python:

```python
import matplotlib.pyplot as plt
from knotted_graph.applications.phase_map_examples import load_phase_map, plot_phase_map

data = load_phase_map("my_records.csv", family="tib2_d6_F")
print(data.summary())
figure = plot_phase_map(data, output="my_results/raw_map")
plt.close(figure)
```

Replace `my_records.csv` with the path to your saved records. Find the
reference CSVs and large research assets in the source checkout.

## Compute a small new example

Install the optional stack, then inspect the plan before doing any computation:

```bash
uv sync --extra nodal --extra viz
uv run --no-sync python -m knotted_graph.applications.phase_map_examples scan materials \
  --profile quick --output-dir _build/new_phase_maps/first_material --dry-run
```

The default quick profile uses one family, a $24^3$ sampling grid and $3\times3$
parameter cells. It uses one process and caps exact-Yamada attempts at eight
core edges. This coarse profile demonstrates the input, extraction and
recording steps. Runtime depends on graph and projection complexity. Use the
research configuration and resolution checks for quantitative comparisons.

To execute the optional scan, remove `--dry-run`:

```bash
uv run --no-sync python -m knotted_graph.applications.phase_map_examples scan materials \
  --profile quick --output-dir _build/new_phase_maps/first_material
```

Each run writes CSV/JSON cell records, summary/plot files and `run_plan.json`. Choose a new or empty output
directory: the guided command refuses to overwrite existing results. Use
`--family` (repeatable), `--dimension`, `--lambda-count`, `--level-count` and
`--max-exact-yamada-edges` to adjust the scan. Material families are `tib2_d6_F`,
`co2mnga_t8`, `ti3al_M2` and `yh3_m1`. Material scans support `--workers`.

For cluster runs, request resources through your site's scheduler. Estimate
the time and memory for larger scans from a coarse run.

<a id="what-a-color-doesand-does-notmean"></a>
<a id="what-a-color-does-and-does-not-mean"></a>

## Read colors and classification status

| Record or display status | Interpretation |
| --- | --- |
| `source=yamada` | a Yamada result for the extracted finite-grid graph; also inspect `classification_computed` |
| `source=vertex` | the engine classified a vertex-only core; inspect components and boundary contact |
| `source=large-core` | a structural signature because an exact attempt was outside the configured limit; not an exact polynomial |
| `source=yamada-set` | an outer/inner boundary classification collection; inspect each entry because structural fallbacks may also occur |
| `source=error`, or nonempty `error` | computation failed; this is not a zero invariant |
| `classification_computed=False`, no calibration anchor | assigned by adaptive energy filling, not independently evaluated at that cell |
| nonempty `resolution_calibration_energy` | assigned using the recorded anchor energy under the upstream TiB2 resolution rule; separate from adaptive filling |
| small-island smoothing or manual signature merges | historical display postprocessing; distinct raw results may share a displayed class |
| TiB2 C6 display grouping | groups audited non-C6 representatives under the stated upstream rule; raw signatures remain available |
| contraction mode | grouping under stated bounded contraction tests; not literal equality of classic Yamada polynomials |

The supplied material records contain 14,634 directly classified, 31,184
adaptive-fill and 202 resolution-calibration cells. A missing classification status is reported as **unrecorded**.
The guided scans classify every requested cell and disable adaptive filling,
resolution calibration, smoothing, C6 display grouping and manual signature
merges. They still use the upstream volume-resolution rule when extracting a
material body; `removed_component_voxels` and `filled_void_voxels` describe that
operation. Compare boundary contact, sampling resolution, components and
the representative skeleton when interpreting classifications.

## Research reproduction and provenance

The dataset guide records the reference assets from upstream `56bfbab` and
the boundary-volume helpers from `643fef8`. Reusable compute engines live in
`knotted_graph.applications.phase_map_examples`; the script paths are
compatibility entry points. The research scripts provide publication layouts
and bounded contraction processing; reusable scan and record-handling functions
are documented in the application API.

Detailed material runs retain explicit research processing options:
`--apply-signature-merges` enables manual display merges, `--apply-c6-review`
enables the TiB2 audit grouping, and `--apply-resolution-calibration` applies
the upstream thin-tube calibration to the records. Adaptive energy scans also
use the upstream per-column stabilization rule; inspect anchor energies to
distinguish its assignments. These options serve different purposes and are
disabled in the guided route.

For `--reuse-records` or `--extend-existing-records`, work on a copy of the
records and matching summary in one output directory. Reuse rebuilds outputs
and records the chosen processing in the summary; extension computes missing
energy rows. The representative-geometry script follows the resolved dominant
body and its outer/nested fillings. The CLI reports a migration error for
`--primary-component-min-fraction`; use the resolved-body selection described above.

`--profile paper` selects denser sampling defaults. Choose families and display
processing explicitly using the options above. For the publication layout, use
the [research-script map](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/NewPhaseMapPlots/README.md)
and the saved metadata. Compare regenerated figures with the reference before
updating a published result.
