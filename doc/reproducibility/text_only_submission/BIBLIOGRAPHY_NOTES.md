# Bibliography corrections for the original project

The original `GeneralVersion/references.bib` was not attached. The revised main source preserves its keys and adds only `text_only_references.bib` for the standard primary Morse-theory reference. Therefore no claim is made that the original bibliography has already been edited or freshly compiled.

The attached 18 September audit supports the following changes to that original file:

| Key | Action |
|---|---|
| `kg_galeski2022lifshitz` | Add article number `7418` to Nature Communications **13** (2022), DOI `10.1038/s41467-022-35106-7`. |
| `kg_shi2017lifshitz` | Add article number `14988` to Nature Communications **8** (2017), DOI `10.1038/ncomms14988`. |
| `kg_li2019density` | Add DOI `10.1134/S0081543819030076` if absent. |
| `kg_shapely` | Cite the actual version/year from the relevant original environment, not the currently installed audit environment. The historical version has not been recovered here. Do not invent a year or version-specific DOI merely to suppress a warning. |
| `kg_dobrynin1996yamada` | Choose the Russian original or the circulated English translation, then make author order, publication details and translation note consistent. The audit reports Vesnin–Dobrynin for the original and Dobrynin–Vesnin for the English translation. Do not combine the two conventions without explanation. |
| `kg_lundstrom2022transfer` | Retain the master's-thesis type and add a verified institutional record URL when available. No new institutional URL was fabricated for this revision. |

Publication years should not be changed solely to match an online-first date. The internal key `kg_chbili2016` need not match the year of the cited 2012 preprint. The citation keys in the source do not need renaming for cosmetic consistency.

The source prose was corrected for citation scope: Ishii/Iwakiri is not used to justify abstract graph contraction as spatial equivalence; Ramer/Douglas–Peucker is not cited as an unconditional isotopy guarantee; Dziawa is an example of tunable band topology, not direct validation of a Lifshitz line in these maps; and TPMS response papers remain motivation rather than evidence for this paper's historical boundaries.

`apply_available_bib_fixes.py` performs only the three unambiguous numerical/DOI additions above, writing a new output file and a change log. It leaves the unresolved version/translation/institutional-link choices to the authors. It never overwrites the input bibliography.
