# Yamada Formula Discovery

Use this **publication-reproduction workflow** to explore graph-family
formulas and their held-out checks. For a first invariant calculation, follow
{doc}`../quickstart` and {doc}`../user_guide/projection_yamada`.

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
   representatives, compute exact Laurent polynomials, and reserve held-out
   values for testing proposed formulas.
2. **Abelian mixed theta words** — test whether the closed invariant factors
   through motif counts rather than order.
3. **Non-Abelian ordered pure-braid words** — reconstruct exact symbolic
   transfer data and test complete held-out Laurent polynomials coefficient by
   coefficient.

The held-out tables test the displayed identities coefficient by coefficient
over the listed families. To establish an all-parameter formula, derive the
general transfer identities from the Yamada skein algebra.

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
of the larger tables. Count cases from the full tables to obtain the totals.

Verify the saved hashes and table/JSON agreement from a source checkout:

```bash
uv run --no-project python scripts/inspect_paper_data.py
```

This command checks the saved records. The dataset guide documents source
versions and the corresponding notebook output schemas.

## Browse versus regenerate

The setup records the reference source commit and reports version differences
for exploratory runs.

Strict publication regeneration additionally requires the audited source
ancestry, a clean `src/knotted_graph` tree, and the factorized native backend:

```bash
export KNOTTEDGRAPH_STRICT_PUBLICATION_REGENERATION=1
```

Reading saved records and exploring the derivation use the standard mode.
Dataset generation and extreme held-out calculations are larger research
workloads; choose their settings in the notebook before running them.

## How to read a successful run

- inspect the environment, commit, backend, and dataset hashes before formulas;
- distinguish training/discovery rows from frozen held-out rows;
- compare full expressions for exact symbolic equality;
- record errors or missing backends as regeneration failures; and
- preserve generated tables and their provenance together.
