# Paper figures, results and data

Explore **[KnottedGraph: Scalable knotted-graph topology for scientific and
mathematical discovery](https://arxiv.org/abs/2609.31152v1)** through its figures.
Each entry below connects the published result to the relevant interface,
notebook and available supporting data. Figure numbers refer to arXiv v1.
The PDF originals were copied unchanged from that version's author source;
the PNGs are display previews. Download the
{download}`figure provenance and SHA-256 inventory <assets/paper_figures/provenance.json>`.

**Start with {doc}`benchmarks` to verify numerical statements from the saved
CSVs, or {doc}`sanity_checks` to run a small correctness check.**
Please {doc}`cite our work <citing>` when using the library or its applications.

## Figure 1: scientific inputs meet at a common graph

```{figure} assets/paper_figures/figure-1.png
:alt: Paper Figure 1 showing scientific input representations entering the knotted-graph pipeline

Different input representations enter the workflow at the stage appropriate
to their existing information. {download}`Open the original PDF <assets/paper_figures/figure-1.pdf>`.
```

Follow {doc}`user_guide/input_adapters` for coordinate, molecular, polymer,
spatial-CSV and surface inputs. The {doc}`feature_status` matrix distinguishes
public loaders from application-specific constructions; an input pictured in
the paper does not automatically imply a generic file parser.
{doc}`user_guide/workflow_overview` explains the shared node `pos` / edge `pts`
representation and the next projection step.

## Figure 2: operations on the shared representation

```{figure} assets/paper_figures/figure-2.png
:alt: Paper Figure 2 illustrating graph construction, layout, projection and invariant functionality

Inspect the intermediate graph and diagram as the representation moves through
the pipeline. {download}`Open the original PDF <assets/paper_figures/figure-2.pdf>`.
```

The working entry points are {doc}`quickstart`,
{doc}`user_guide/projection_yamada`, {doc}`user_guide/repulsive_layout` and
{doc}`api/index`. The supplementary figure directory below links the individual
extraction, PD-encoding and evaluation stages to their documentation.

## Figure 3: exact evaluation through 500 crossings

```{figure} assets/paper_figures/figure-3.png
:alt: Paper Figure 3 showing reference graph families and KnottedGraph versus Topoly runtime curves

Structured families are checked against published polynomials before their
timings are interpreted. {download}`Open the original PDF <assets/paper_figures/figure-3.pdf>`.
```

**Evidence:** {download}`112-row scaling CSV <../User_guide/benchmarks/results/03_knottedgraph_vs_topoly_scaling_rows.csv>`.
**Procedure:** [benchmark notebook 03](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/benchmarks/03_knottedgraph_vs_topoly_scaling.ipynb).

All saved KnottedGraph rows pass the reference comparison. The three
500-crossing evaluation times are approximately 0.117, 2.376 and 4.998 seconds.
The {doc}`benchmarks` page also reproduces the seven-crossing Table 1 values
from these rows and explains completed, errored and skipped Topoly calls.
The high-crossing endpoints are not completed paired Topoly comparisons.

## Figure 4: discovering and testing family laws

```{figure} assets/paper_figures/figure-4.png
:alt: Paper Figure 4 showing homogeneous motifs, commuting mixed words and noncommuting ordered braid words

Homogeneous motifs, mixed families and ordered braid words connect exact
polynomial calculations with candidate family formulas.
{download}`Open the original PDF <assets/paper_figures/figure-4.pdf>`.
```

**Procedure:** [application notebook 05](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/05_yamada_formula_discovery.ipynb),
with the three parts explained in {doc}`applications/yamada_formula_discovery`.

**Download the data:**
{download}`Complete data and source ZIP <assets/data/figure4-suppfig10-source-data.zip>` ·
[Dataset guide and provenance](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/results/figure4_suppfig10/README.md).
The supplied exact-check records can be browsed directly, without running a notebook.

| Figure panels | Included records | Download |
| --- | --- | --- |
| 4b,c,e,f,h,i: homogeneous examples | 6 examples across three families, m = 1 and 2 | {download}`CSV <../User_guide/applications/results/figure4_suppfig10/data/figure4_homogeneous_m1_m2.csv>` |
| 4k: mixed families | 459 coefficient comparisons | {download}`CSV <../User_guide/applications/results/figure4_suppfig10/data/figure4_mixed_family_459.csv>` |
| 4l: ordered braid words | 324 short-word records | {download}`CSV <../User_guide/applications/results/figure4_suppfig10/data/figure4_pure_braid_324.csv>` |
| 4l: order-sensitive examples | AAB, ABA and BAA | {download}`CSV <../User_guide/applications/results/figure4_suppfig10/data/figure4_order_sensitive_AAB_ABA_BAA.csv>` |
| Additional long-word comparisons | Two completed length-101 cases | {download}`JSON <../User_guide/applications/results/figure4_suppfig10/audit/discovery_long_words/records.json>` |

Each coefficient table links through the dataset guide to its JSON records,
verification certificate and source information. The small example tables are
subsets of the larger tables. The {doc}`formula-discovery guide <applications/yamada_formula_discovery>`
explains how to read them.

## Supplementary Figure 10: an ordered-word worked example

```{figure} assets/paper_figures/supp-10.png
:alt: Supplementary Figure 10 working through the AAABA ordered-braid example

The AAABA = A³BA example connects the graph, ordered transfer product and
Yamada coefficients. {download}`Open the original PDF <assets/paper_figures/supp-10.pdf>`.
```

**Data:** {download}`AAABA coefficient record (JSON) <../User_guide/applications/results/figure4_suppfig10/data/suppfig10_AAABA.json>` ·
{download}`Exact transfer matrices (JSON) <../User_guide/applications/results/figure4_suppfig10/audit/discovery/exact_transfer_matrices.json>` ·
{download}`Verification certificate (JSON) <../User_guide/applications/results/figure4_suppfig10/audit/discovery/certificate.json>`.

This worked example is included in the 324-word dataset above. Its JSON entry
identifies the graph and the corresponding Hankel matrix input.

## Supplementary Figures 6–7: construction and runtime evidence

```{figure} assets/paper_figures/supp-6.png
:alt: Representative thick-handlebody recovery cases from Supplementary Figure 6

Supplementary Figure 6 illustrates the inverse volume-to-graph benchmark.
{download}`Original PDF <assets/paper_figures/supp-6.pdf>`.
```

The {download}`4,400-case preservation table <../User_guide/benchmarks/results/handlebody_ground_truth/synthetic_ground_truth_yamada_preservation.csv>`
contains the recorded input acceptance, topology diagnostics and normalized
Yamada comparisons. See {doc}`benchmarks` for the distinction between invariant
preservation and abstract isomorphism, and for geometry/manifests.

```{figure} assets/paper_figures/supp-7.png
:alt: Nine pipeline timing distributions in Supplementary Figure 7

Supplementary Figure 7 breaks down the same construction workload by stage.
{download}`Original PDF <assets/paper_figures/supp-7.pdf>`.
```

The {download}`4,400-row timing-figure source <../User_guide/benchmarks/results/04_thick_handlebody_time_distribution_plot_source.csv>`
maps to panels a–i in the {doc}`benchmarks` field guide. Topoly's 3,968 completed
calls, 411 timeouts and 21 errors are recorded separately; its histogram uses
completed calls only.

## Complete supplementary figure directory

| Figure | Original PDF | Follow the method or inspect the data |
| --- | --- | --- |
| S1: Input and Yamada examples | {download}`PDF <assets/paper_figures/supp-1.pdf>` | {doc}`user_guide/input_adapters`; {doc}`sanity_checks` |
| S2: Skeletonization and extraction | {download}`PDF <assets/paper_figures/supp-2.pdf>` | {doc}`user_guide/workflow_overview`; benchmark 04 linked above |
| S3: PD-code generation | {download}`PDF <assets/paper_figures/supp-3.pdf>` | {doc}`user_guide/projection_yamada` |
| S4: PD code to Yamada | {download}`PDF <assets/paper_figures/supp-4.pdf>` | {doc}`user_guide/projection_yamada`; {doc}`api/yamada` |
| S5: Yamada resolutions | {download}`PDF <assets/paper_figures/supp-5.pdf>` | {doc}`api/yamada`; {doc}`sanity_checks` |
| S6: Handlebody recovery | {download}`PDF <assets/paper_figures/supp-6.pdf>` | {doc}`benchmarks`: preservation CSV and geometry manifests |
| S7: Timing distributions | {download}`PDF <assets/paper_figures/supp-7.pdf>` | {doc}`benchmarks`: stage-by-stage CSV column map |
| S8: Parameter sweeps | {download}`PDF <assets/paper_figures/supp-8.pdf>` | {doc}`applications/hamiltonian_yamada_phase_maps`: saved interactive view and notebook |
| S9: Repulsive layout | {download}`PDF <assets/paper_figures/supp-9.pdf>` | {doc}`user_guide/repulsive_layout`; {doc}`applications/analytic_knot_fields` |
| S10: Non-Abelian worked example | {download}`PDF <assets/paper_figures/supp-10.pdf>` | [Worked example and data](#supplementary-figure-10-an-ordered-word-worked-example); notebook 05, Part III |
| S11: Software architecture | {download}`PDF <assets/paper_figures/supp-11.pdf>` | {doc}`user_guide/workflow_overview`; {doc}`api/index`; {doc}`installation` |

The archived PDFs show the submitted presentation. The notebooks expose
scientific procedures; a saved record and a notebook are not necessarily a
pixel-identical figure compositor. The figure provenance identifies the exact
PDF version available here.

## Scientific application data and interactive views

For the separate [scientific application paper](https://arxiv.org/abs/2609.13390),
{doc}`applications/material_phase_maps` exposes saved TiB2/Co2MnGa records and
an interactive region viewer. It distinguishes directly classified data,
adaptive fills, resolution calibration and historical display processing.
{doc}`applications/hamiltonian_yamada_phase_maps` provides the retained nodal
viewer and its provenance. These pages allow exploration of existing data
without starting a new parameter scan.
