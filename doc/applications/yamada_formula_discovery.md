# Yamada Formula Discovery

This is an **advanced publication-reproduction workflow**, not the first place
to learn the Yamada API. Start with {doc}`../quickstart` and
{doc}`../user_guide/projection_yamada` if you want to compute one invariant.

<div class="kg-hero">
  <p class="kg-lead">Explore exact Laurent-polynomial data for homogeneous families, commuting mixed words and order-sensitive braid words. Start with the saved CSV/JSON records, then follow the notebook to inspect the construction and formula procedures.</p>
  <div class="kg-link-row">
    <a href="https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/05_yamada_formula_discovery.ipynb">Open the notebook</a>
    <a href="../paper_results.html#figure-4-discovering-and-testing-family-laws">Figures and data downloads</a>
    <a href="../user_guide/projection_yamada.html">Review projection and Yamada concepts</a>
  </div>
</div>

## What the notebook is testing

The notebook has three independent scientific parts:

1. **Homogeneous theta-derived families** — construct certified spatial
   representatives, compute exact Laurent polynomials, and reserve values that
   are not used during formula proposal.
2. **Abelian mixed theta words** — test whether the closed invariant factors
   through motif counts rather than order.
3. **Non-Abelian ordered pure-braid words** — reconstruct exact symbolic
   transfer data and test complete held-out Laurent polynomials coefficient by
   coefficient.

The held-out calculations are evidence for the displayed identities over the
tested families. They are not, by themselves, an all-parameter mathematical
proof; an analytic proof must still derive the transfer identities from the
Yamada skein algebra.

## Browse the saved data

The {doc}`../paper_results` page places Figure 4 and Supplementary Figure 10
beside their CSV/JSON downloads. The supplied dataset includes **459 mixed-family
comparisons**, **324 short pure-braid words**, the **AAABA worked example** and
**two completed length-101 comparisons**.

{download}`Download all data and source files (ZIP) <../assets/data/figure4-suppfig10-source-data.zip>`
or open the [dataset guide](https://github.com/HakanAkgn/KnottedGraph/blob/main/User_guide/applications/results/figure4_suppfig10/README.md)
for column definitions, panel mappings, certificates and provenance.

Columns ending in `_coefficients_json` store exact integer coefficients keyed
by Laurent exponent. `master_formula_pass` records the mixed-family comparison;
the long-word JSON uses `raw_coefficient_identity`. The example CSVs are subsets
of the larger tables, so they should not be added to the total case count.

Verify the saved hashes and table/JSON agreement from a source checkout:

```bash
uv run --no-project python scripts/inspect_paper_data.py
```

This command and the website build read the archived records without starting
new scientific calculations. Source versions and the distinction between these
records and the notebook's own output schemas are documented in the dataset guide.

## Browse versus regenerate

Opening and reading the notebook no longer depends on a literal Git branch
name. The setup records the audited source commit and warns when the current
checkout differs. This is appropriate for review and exploratory execution.

Strict publication regeneration additionally requires the audited source
ancestry, a clean `src/knotted_graph` tree, and the factorized native backend:

```bash
export KNOTTEDGRAPH_STRICT_PUBLICATION_REGENERATION=1
```

Do not enable strict mode merely to read the derivation. The dataset-generation
and extreme held-out stages can be expensive and are intentionally excluded
from automatic beginner-notebook execution.

## How to read a successful run

- inspect the environment, commit, backend, and dataset hashes before formulas;
- distinguish training/discovery rows from frozen held-out rows;
- require exact symbolic equality, not numerical samples at selected values;
- treat an error or missing backend as a failed regeneration, not as a
  topological result; and
- preserve generated tables and their provenance together.

The notebook uses explicit display-math delimiters so the symbolic summary
renders consistently in Jupyter, GitHub, and the documentation toolchain.
