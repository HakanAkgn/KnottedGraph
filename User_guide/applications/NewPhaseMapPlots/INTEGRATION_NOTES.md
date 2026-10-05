# Material phase-map provenance

The material reference assets originate from `56bfbab`, with the
boundary-topology helpers completed in `643fef8`. The tables and viewers preserve
the recorded research calculations and their processing metadata.

## What the saved records mean

The material CSV contains 46,020 rows: 14,634 directly classified records,
31,184 adaptive fills and 202 resolution calibrations. TiB2 has 35 raw
signatures and Co2MnGa has 53. The saved viewer contains 88 region entries
across seven transitions, including five nodal transitions. Viewer entries
count displayed regions. Compare polynomial values and computation status in
the raw cell records.

The `source`, status and processing metadata distinguish exact computations,
adaptive fills, resolution calibration, display merges and contraction groups.
Historical display groups represent the recorded presentation partition and
can contain several distinct polynomial signatures. The guided CLI classifies
each requested cell; specialized scripts expose adaptive filling, display
merging and calibration as selectable processing settings.

## Geometry and presentation

Representative geometry uses the dominant resolved volume body and its outer
and nested boundary fillings. `dev/build_phase_map_demos.py` splits the saved
viewer into lazily loaded region attachments without regenerating the scan.
Absolute paths inside the original research HTML record the producing machine.
The walkthrough commands use repository-relative inputs and explicit outputs.

The [material walkthrough](../../../doc/applications/material_phase_maps.md)
explains saved-data inspection, raw plotting, and optional new calculations.
The [paper guide](../../../doc/paper_results.md) identifies which figures and
records belong to the software paper and which illustrate the separate
scientific application paper.
