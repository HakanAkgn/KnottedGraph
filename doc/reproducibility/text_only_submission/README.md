# KnottedGraph — text-only arXiv revision

This replaces the manuscript text supplied on 19 September 2026, corresponding to `Knotted_Graph_Code_Paper-22.pdf`. It does not replace any figure PDF, internal artwork, numerical plot or scientific implementation. The TPMS section and its unchanged figure move from the main text to the supplementary application note. All seventeen original figure references are retained.

## What to use

`Third_draft.tex` in the generated package is the **complete manuscript source**, including the main text, all supplementary sections, all figure calls and the unchanged inline TikZ drawing. It is not an insertion fragment. The `parts/` directory in the repository is the editing source from which this full file is concatenated.

Place `Third_draft.tex` and `text_only_references.bib` at the Overleaf/project root that contains `GeneralVersion/`. Keep your existing `GeneralVersion/Figures/`, `GeneralVersion/references.bib` and `GeneralVersion/macros.tex`. The revised source also permits the original `macros.tex` at the compile root when the GeneralVersion copy is absent. Select `Third_draft.tex` as the main document. No figure image needs replacement.

The single added bibliography entry is the primary Morse-theory book. The supplied text includes the original citation keys. The original `.bib`, macro file and seventeen standalone figure assets were **not supplied with this request**, so the generated source package is not a self-contained arXiv upload until those existing project files are included. Do not upload the standalone source ZIP as though it contained the complete original artwork.

## Verification scope

The assembly script checks all seventeen figure references and their original size options; requires the known retained TikZ drawing verbatim; checks environment and brace balance, duplicate labels and textual cross-reference targets; records citation keys and author placeholders; and exports the shortened plain-text abstract. It does not claim a fresh LaTeX build or visual approval of a new PDF without the original macro, bibliography and figure dependencies.

The user-provided compiled PDF was consulted as the appearance/reference document. Figures were neither extracted, regenerated nor modified. Caption corrections and movement of the TPMS float are text/layout changes. The bad RII overpass remains in the drawing and is explicitly identified as an invalid cancellation in the adjacent text. The historical maps likewise remain filtered historical illustrations rather than being relabeled as newly verified phase maps.

## Before submission

1. Fill the intentional author, affiliation, correspondence, funding and supplementary-date placeholders.
2. Apply the small bibliography corrections in `BIBLIOGRAPHY_NOTES.md` to your existing `.bib`; its unavailable contents were not silently fabricated.
3. Compile the full existing project with the revised source, then inspect the actual submitted PDF and arXiv preview. A text-only revision cannot certify the old numerical displays or guarantee acceptance.

The scientific formulas are retained with corrected domains, proof-status and verification claims. The cavity-treatment direction is preserved; no new cavity method or numerical rerun is claimed. The main changes are catalogued in `AUDIT_RESPONSE.md`.
