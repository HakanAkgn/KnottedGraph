# Figure 4 and Supplementary Figure 10 data

Browse the exact coefficient records accompanying the paper's homogeneous,
mixed-family and ordered-braid examples. This release organizes the data package
provided on 29 September 2026; its recorded checks were run on 18 September 2026.

**[Download the complete data and source ZIP](../../../../doc/assets/data/figure4-suppfig10-source-data.zip)**
or [view the figures beside their downloads](https://hakanakgn.github.io/KnottedGraph/paper_results.html#figure-4-discovering-and-testing-family-laws).

## Choose a dataset

| Figure / example | Records | Download | What to inspect |
| --- | ---: | --- | --- |
| Figure 4b,c,e,f,h,i: homogeneous examples | 6 | [CSV](data/figure4_homogeneous_m1_m2.csv) | Three families at m = 1 and 2; `figure4_panels` identifies each panel. |
| Figure 4k: mixed families | 459 | [CSV](data/figure4_mixed_family_459.csv) | Word, motif counts, actual/predicted normalized coefficients and comparison result. |
| Figure 4l: short pure-braid words | 324 | [CSV](data/figure4_pure_braid_324.csv) | Word, graph size, crossings, raw and normalized coefficients. |
| Figure 4l: order-sensitive examples | 3 | [CSV](data/figure4_order_sensitive_AAB_ABA_BAA.csv) | AAB, ABA and BAA extracted from the short-word table. |
| Supplementary Figure 10: AAABA | 1 | [JSON](data/suppfig10_AAABA.json) | The A³BA example, including its direct coefficient record and Hankel input entry. |
| Additional length-101 comparisons | 2 completed | [JSON](audit/discovery_long_words/records.json) | A⁵¹B⁵⁰ and (AB)⁵⁰A: stored predictions compared with direct results. |

The 6-row and 3-row tables are subsets of the 459-row and 324-row tables,
respectively. The AAABA record also occurs in the 324-row table.

## Read the columns

- `word` is the case identifier. In the mixed-family tables, L means cross-linked,
  D means lower-paired and B means two-sided. In the pure-braid tables, A and B
  are braid blocks.
- Columns ending in `_coefficients_json` contain a JSON object mapping integer
  Laurent exponents to exact integer coefficients. For example,
  `{"-1": 2, "0": -3, "2": 1}` means 2Y⁻¹ − 3 + Y².
- `master_formula_pass` records equality of the mixed-family coefficient
  dictionaries. `status` records whether a calculation completed.
- The long-word JSON records use `raw_coefficient_identity` for the comparison;
  each `actual` object contains the independently evaluated record.

## Check the saved files

From the repository root, run:

```bash
uv run --no-project python scripts/inspect_paper_data.py
```

This reads files, verifies the hashes in [manifest.json](manifest.json), and
checks that the compact tables match the archived JSON records. It also reports
the saved benchmark counts. It does not evaluate a new graph or polynomial.

The included [table exporter](REGENERATE_EXTRACTS.py) converts the JSON records
back to the four CSV files and AAABA JSON. Run it on a copy if you want to
recreate the compact exports; its five outputs were checked byte-for-byte
against the supplied files before publication.

## Source records and verification

| Group | Records and definitions | Verification | Source information |
| --- | --- | --- | --- |
| Mixed families | [Records](audit/mixed_family/records.json) · [Word list](audit/mixed_family/word_plan.json) | [Certificate](audit/mixed_family/certificate.json) | [Provenance](audit/mixed_family/provenance.json) |
| Pure-braid words | [Records](audit/discovery/short_exact_records.json) · [Word lists and basis](audit/discovery/word_plans.json) · [Transfer matrices](audit/discovery/exact_transfer_matrices.json) · [Hankel matrices](audit/discovery/hankel_input_matrices.json) | [Certificate](audit/discovery/certificate.json) · [Short-word comparisons](audit/discovery/exhaustive_short_word_checks.json) | [Provenance](audit/discovery/provenance.json) |
| Long words | [Records](audit/discovery_long_words/records.json) · [Plan](audit/discovery_long_words/word_plan.json) | [Certificate](audit/discovery_long_words/certificate.json) | [Provenance](audit/discovery_long_words/provenance.json) |

The complete ZIP contains the supplied notebook, verification scripts,
manuscript snapshot, protocol notes and original SHA-256 inventory. The source
notebook uses the graph constructions and evaluation procedures described in
[application notebook 05](../../05_yamada_formula_discovery.ipynb). To run the
scripts, use a full KnottedGraph checkout with its dependencies installed; use
the ZIP to browse the source and saved records.

## Record provenance

| Record | Provenance and interpretation |
| --- | --- |
| Exact comparisons | Retrospective checks of the recorded cases. The certificates describe these finite comparisons; an all-family identity requires an analytical derivation. |
| Historical notebook outputs | The original `05a_*`, `05b_*` and `05c_*` ledgers and an independent historical formula-freezing record are absent from this package. The compact tables use their own schemas; use the notebook's output schemas when regenerating its ledgers. |
| AAABA | A worked example included in the Hankel input. It belongs to the construction data rather than an independent held-out set. |
| Length-101 plan | Two completed comparisons from the 20-case plan: A⁵¹B⁵⁰ and (AB)⁵⁰A. The remaining 18 cases have no completed result in the supplied certificate. |

Detailed protocols and the original handoff notes are included in the ZIP.
Browse these alongside the figures and coefficient records to follow the
recorded comparisons.
