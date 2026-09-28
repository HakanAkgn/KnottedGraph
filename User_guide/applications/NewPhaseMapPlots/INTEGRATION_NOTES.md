# Material phase-map provenance

The retained material reference assets were introduced in `56bfbab`; the
boundary-topology helpers were completed in `643fef8`. These are existing
research results, not new calculations from the publication documentation update.

## What the saved records mean

The material CSV contains 46,020 rows: 14,634 directly classified records,
31,184 adaptive fills and 202 resolution calibrations. TiB2 has 35 raw
signatures and Co2MnGa has 53. The saved viewer contains 88 region entries
across seven transitions, including five earlier nodal transitions. A region
count is not a count of distinct polynomials.

The `source`, status and processing metadata distinguish exact computations,
adaptive fills, resolution calibration, display merges and contraction groups.
Historical display grouping does not assert equality of the grouped polynomials.
The guided CLI leaves adaptive filling, manual merging and calibration off for
new scans; specialized scripts expose those choices explicitly.

## Geometry and presentation

Representative geometry uses the dominant resolved volume body and its outer
and nested boundary fillings. `dev/build_phase_map_demos.py` splits the saved
viewer into lazily loaded region attachments without regenerating the scan.
Absolute paths inside the original research HTML record the producing machine;
they are not installation requirements.

The [material walkthrough](../../../doc/applications/material_phase_maps.md)
explains saved-data inspection, raw plotting, and optional new calculations.
The [paper guide](../../../doc/paper_results.md) identifies which figures and
records belong to the software paper and which illustrate the separate
scientific application paper. Earlier integration notes and their historical
validation reports remain available in Git history.
