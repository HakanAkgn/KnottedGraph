# Material And Compact-TPMS Phase Maps

Start here for the material-parameter and compact-scaffold examples supplied in
`User_guide/applications/NewPhaseMapPlots`. You can inspect saved results without
running a scan. These are application-specific examples, not new generic input
formats or a replacement for the [phase-map API](../api/applications.md).

## Choose how much work to run

| Goal | Entry | What runs / prerequisites |
| --- | --- | --- |
| Understand the published examples | previews and interactive links below | browser; Plotly needs access to its pinned CDN |
| Inspect or replot saved cells | `phase_map_examples inspect` / `plot` | base installation; reads CSV/JSON, no extraction |
| Learn the computation | `phase_map_examples scan --profile quick` | optional `nodal` and `viz`; coarse grid on a compute node |
| Reproduce high-resolution research plots | source scripts and reference data | advanced; expensive scans, explicit display processing and validation |

Follow the [source-install instructions](../installation.md) for this review
branch. The legacy PyPI release does not contain these commands. Commands below
assume the repository root; the installed Python helpers also work outside a
checkout when supplied with your own records file.

## Saved research results

The saved displays below are historical artifacts. The [scientific revision](../reproducibility/arxiv_revision.md) supplies corrected unfiltered TPMS data, cavity and component checks, and explicit evaluation scopes. Its results supersede the historical TPMS class interpretation. [Reproduction commands and filtering audits](../reproducibility/tpms.md) retain the original display records separately.

The material examples interpolate one Hamiltonian coefficient with $\lambda$
and vary an energy/gap threshold $E$. The compact-TPMS examples interpolate
between gyroid, Schwarz-P and diamond scalar fields with $\lambda$, then vary a
threshold $c$ inside a finite enclosing domain. This is a compact scaffold cut
from TPMS fields, not an infinite periodic structure with periodic-boundary
topology automatically identified.

<div class="kg-wide-figure">
  <img src="../demos/new_phase_maps/materials/preview.png" alt="Saved material-parameter phase-map panels for TiB2 and Co2MnGa">
</div>

<p><a href="../demos/new_phase_maps/materials/index.html">Open material and Hamiltonian regions</a>.</p>
The saved material artifact also includes the earlier nodal transitions; select
TiB2 or Co2MnGa to see the new material scans.

<div class="kg-wide-figure">
  <img src="../demos/new_phase_maps/tpms/preview.png" alt="Saved compact-TPMS phase-map overview for the three field interpolations">
</div>

<p><a href="../demos/new_phase_maps/tpms/index.html">Open compact-TPMS regions</a>.</p>
Choose a transition and classification mode, then click a region to see its
representative surface and skeleton. A representative is one selected cell,
not every geometry in that colored region. In “Up to contraction” views,
read the reported grouping convention before comparing colors.

These links load geometry **on selection**, not when you open this guide. If
Plotly or a region download is blocked, the viewer reports the problem; the
static figures remain available. Downloaded split viewers must be served over
HTTP, not opened with `file://`. The original self-contained research HTML files
remain in the source tree, but still require their declared browser dependencies.

## First exercise: inspect saved cells

No optional extraction stack is needed:

```bash
uv run python -m knotted_graph.applications.phase_map_examples inspect \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/data/material_parameter_phase_map_records.csv
```

Expected reference counts:

| Family key | Records | Grid | Distinct raw signatures |
| --- | ---: | --- | ---: |
| `tib2_d6_F` | 30,180 | 60 lambda × 503 energy samples | 35 |
| `co2mnga_t8` | 15,840 | 60 lambda × 264 energy samples | 53 |
| `gyroid_to_diamond` | 441 | 21 lambda × 21 threshold samples | 65 |
| `gyroid_to_schwarz_p` | 441 | 21 lambda × 21 threshold samples | 20 |
| `schwarz_p_to_diamond` | 441 | 21 lambda × 21 threshold samples | 15 |

For TPMS, substitute
`User_guide/applications/NewPhaseMapPlots/TPMS/data/tpms_parameter_phase_map_records.csv`.
The inspector also reports errors, missing cells, source counts, direct
classifications, adaptive fills and resolution calibrations. It never
interprets a signature string as executable Python.

## Replot without recomputing

```bash
uv run python -m knotted_graph.applications.phase_map_examples plot \
  User_guide/applications/NewPhaseMapPlots/RealMaterials/data/material_parameter_phase_map_records.csv \
  --family tib2_d6_F --output _build/new_phase_maps/tib2_raw
```

This saves `tib2_raw.png`, `tib2_raw.pdf` and `tib2_raw.json`. The JSON is the
phase key: it maps each category ID to its full signature, color, source and
recorded polynomial, and includes the coordinate axes and the category-ID grid.
Category IDs are not polynomial values. The plot uses serif labels
and Computer Modern mathematics without changing your global plotting settings.
White cells are missing, shaded cells are adaptive fills, open circles mark
resolution calibrations, and crosses mark errors. The JSON also records each
cell's classification status and calibration anchor energy when present.
No small-island smoothing, C6 grouping or manual signature merging is applied.
With many categories, use the JSON key rather than relying on color alone.

The same operation from Python:

```python
import matplotlib.pyplot as plt
from knotted_graph.applications.phase_map_examples import load_phase_map, plot_phase_map

data = load_phase_map("my_records.csv", family="gyroid_to_diamond")
print(data.summary())
figure = plot_phase_map(data, output="my_results/raw_map")
plt.close(figure)
```

`my_records.csv` is a placeholder for a saved records file, not a file bundled
with the wheel. Reference CSVs and large research assets are in the source
checkout, not installed into `site-packages`.

## Compute a small new example

Install the optional stack, then inspect the plan before doing any computation:

```bash
uv sync --extra nodal --extra viz
uv run --no-sync python -m knotted_graph.applications.phase_map_examples scan tpms \
  --profile quick --output-dir _build/new_phase_maps/first_tpms --dry-run
```

The default quick profile uses one family, a $24^3$ sampling grid and $3\times3$
parameter cells. It uses one process and caps exact-Yamada attempts at eight
core edges. That cap does not guarantee a fixed runtime: projection complexity
also matters. A low-resolution result is a smoke example, not a converged
topological measurement or a substitute for the saved high-resolution figures.

On a compute node, remove `--dry-run`:

```bash
uv run --no-sync python -m knotted_graph.applications.phase_map_examples scan tpms \
  --profile quick --output-dir _build/new_phase_maps/first_tpms
uv run --no-sync python -m knotted_graph.applications.phase_map_examples scan materials \
  --profile quick --output-dir _build/new_phase_maps/first_material
```

Each run writes CSV/JSON cell records, summary/plot files and `run_plan.json`;
TPMS also writes representative per-cell geometry. Choose a new or empty output
directory: the guided command refuses to overwrite existing results. Use
`--family` (repeatable), `--dimension`, `--lambda-count`, `--level-count` and
`--max-exact-yamada-edges` to adjust the scan. Material families are `tib2_d6_F`,
`co2mnga_t8`, `ti3al_M2` and `yh3_m1`; TPMS families are listed above. Material
scans support `--workers`; TPMS is currently serial.

On a cluster, even the quick scan and rendering/build tasks belong in a batch
or interactive compute allocation, not on a login node. Resource choices are
site-specific; no queue, account or scheduler settings are embedded in this
application. Dense scans can be substantially more expensive in both time and
memory. Estimate them from a coarse run before requesting larger resources.

## What a color does—and does not—mean

| Record or display status | Interpretation |
| --- | --- |
| `source=yamada` | in revised TPMS scans, a normalized subcubic spatial-graph polynomial; historical records require their original scope audit |
| `source=diagram-yamada` | a higher-valence fixed-diagram polynomial; retain its projection and local vertex convention |
| `source=vertex` | the engine classified a vertex-only core; inspect components and boundary contact |
| `source=large-core` | a structural signature because an exact attempt was outside the configured limit; not an exact polynomial |
| `source=yamada-set` | an outer/inner boundary classification collection; inspect each entry because structural fallbacks may also occur |
| `source=error`, or nonempty `error` | computation failed; this is not a zero invariant |
| `classification_computed=False`, no calibration anchor | assigned by adaptive energy filling, not independently evaluated at that cell |
| nonempty `resolution_calibration_energy` | assigned using the recorded anchor energy under the upstream TiB2 resolution rule; separate from adaptive filling |
| small-island smoothing or manual signature merges | historical display postprocessing; distinct raw results may share a displayed class |
| TiB2 C6 display grouping | groups audited non-C6 representatives under the stated upstream rule; raw signatures remain available |
| historical contraction mode | abstract graph grouping that discards embedding information; not certified spatial or handlebody equivalence |

The supplied material records contain 14,634 directly classified, 31,184
adaptive-fill and 202 resolution-calibration cells. TPMS records do not carry `classification_computed`; the
inspector reports this as **unrecorded**, rather than silently assuming a value.
The new guided scans classify every requested cell and disable adaptive filling,
resolution calibration, smoothing, C6 display grouping and manual signature
merges. They still use the upstream volume-resolution rule when extracting a
material body; `removed_component_voxels` and `filled_void_voxels` describe that
operation. Inspect boundary contact, sampling resolution,
components and the representative skeleton before interpreting any class.

## Research reproduction and provenance

The latest material reference data, figures and HTML from upstream `56bfbab`
are preserved, together with the TPMS assets from `2b2ae6d`. The boundary-volume
helpers are supplied by upstream `643fef8`. Historical path strings are provenance, not runtime
requirements. Reusable compute engines now live in
`knotted_graph.applications.phase_map_examples`; the old script paths are
compatibility entry points. Specialized publication layouts and bounded
contraction processing remain documented research scripts, not generic public
APIs. Their map mathematics and Hamiltonian definitions were not redesigned.

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
body and its outer/nested fillings. The old
`--primary-component-min-fraction` selection option no longer applies.

`--profile paper` chooses denser sampling defaults, but does **not** silently
enable historical display merges, adaptive fills or all families, and is not
a promise of pixel-identical publication reproduction. For that task use the
[research-script map](https://github.com/sarinstein-yan/KnottedGraph/blob/codex/arbitrary-knot-user-integration/User_guide/applications/NewPhaseMapPlots/README.md),
match the saved metadata, and review regenerated results before replacing an
accepted figure. Input-format figures Main/S1/S2 are unrelated and unchanged.

The website builder separates saved region JSON from the original HTML without
changing grid values or geometry. `dev/build_phase_map_demos.py` checks pinned
source hashes and writes a manifest of region hashes. Sphinx calls it for HTML
builds; generated split assets are ignored by Git, avoiding another checked-in
copy of the large payloads. This branch's preview is not automatically the
published site: deployment still follows the repository's main-branch workflow.
