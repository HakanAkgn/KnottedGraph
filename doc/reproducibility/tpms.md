# TPMS revision reproducibility bundle

These scripts and records support the corrected finite-resolution TPMS analysis. They do not assert complete spatial-graph or handlebody equivalence. The experimental contraction probe is separate and performs no automatic map merging.

Use the repository's pinned `uv` environment. All script input and output paths are arguments; none requires an author's machine path. The historical JSON files retain original absolute source locations as provenance, alongside SHA-256 hashes. They are not paths that must exist to rerun the portable scripts.

Run from the repository root. In the examples, `B` is that root and `D` is the corrected dense scan directory. `R` is a fresh output directory.

```sh
B=.
D=User_guide/applications/results/arxiv_revision/tpms_n64
R=_build/tpms_validation
mkdir -p "$R"

# Reproduce the exact earlier display filter (not the revised optional filter).
uv run python "$B/dev/audit_archived_tpms_filter.py" \
  --extracted-input "$B/User_guide/applications/results/arxiv_revision/tpms_audit/tpms_archived_filter_changed_cells.json" \
  --output "$R/archived_filter.json"

# Audit corrected records and thresholds 1/2/4/8 using the frozen old rule.
uv run python "$B/dev/analyze_tpms_dense_scan.py" \
  --records "$D/tpms_parameter_phase_map_records.json" \
  --filter-source "$B/dev/legacy_tpms_filter.py" \
  --output-prefix "$R/corrected"

# Reproduce parameter-subset sensitivity and matched voxel-resolution counts.
uv run python "$B/dev/analyze_tpms_sampling.py" \
  --records "$D/tpms_parameter_phase_map_records.json" \
  --stage-records "$B/User_guide/applications/results/arxiv_revision/tpms_audit/tpms_extraction_stage_diagnosis.json" \
  --output "$R/sampling.json"

# Rebuild the unfiltered six-panel main figure, not an appendix figure.
uv run python "$B/dev/plot_tpms_audited_main.py" \
  --records "$D/tpms_parameter_phase_map_records.json" \
  --output-prefix "$R/PorousMaterialPhaseMap"

# Recompute the 135-cell extraction audit from the installed source.
uv run python "$B/dev/diagnose_tpms_extraction.py" --out-dir "$R"
uv run python "$B/dev/diagnose_tpms_base_fallback.py" --out-dir "$R"
uv run python "$B/dev/diagnose_tpms_extraction_constrained.py" --out-dir "$R"

# Experimental exact PL contraction tests and bounded archived-data probes.
uv run pytest tests/core/test_embedded_contraction.py
uv run python "$B/dev/evaluate_embedded_contraction.py" \
  --geometry-dir User_guide/applications/NewPhaseMapPlots/TPMS/geometry \
  --out-dir "$R"
uv run python "$B/dev/evaluate_contraction_pairs.py" \
  --tpms-dir User_guide/applications/NewPhaseMapPlots/TPMS \
  --out-dir "$R"
```

The corrected dense scan uses dimension 64, 21 λ samples over [0,1], 21 c samples over [0,.3], spherical compactification radius fraction .72 inside a box of half-width 2.25π, an exact evaluation limit of 18 graph edges, and minimum stable-cell size 1 (unfiltered). The parent revision manifest records its command, source commit and environment. Raw records SHA-256 is `b4a8dc6c4feb8905f2c17e8b0567f5a9112135e537f85f764e5af45c0bf09eb2`.

The corrected summary and findings note distinguish successful Betti checks, unsupported cavity masks, subcubic invariant values, fixed-diagram evaluations and abstract summaries. Filter JSON files include per-cell masks. `tpms_sampling_sensitivity.json` treats nested 21/11/6 grids as sampling sensitivity, not independent convergence. The matched 24³/64³/96³ audit measures mask Betti changes separately from extraction agreement.

The historical extraction-stage files document the unconstrained candidate selection before the optional Betti gate; recomputation with current source keeps the unconstrained default for that probe. The constrained probe explicitly passes expected component and cycle counts and refuses cavity-bearing cases. Source-module hashes are saved with the results. All stage-level findings are empirical; Betti equality is not a general proof of deformation retraction or preserved embedding.

The contraction prototype treats binary floating-point coordinates as exact rationals and certifies only a restricted straight-edge move when its swept regions satisfy disjointness predicates. Bounded comparison returns `unknown` when no witness is found. A graph-rank obstruction is reported separately from lack of a witness. Historical geometry is rounded and one exported graph was downsampled; certificates concern the actual input polylines only. The 20 certified moves do not validate the original smooth solids. Of nine neighboring pair probes, eight are unknown; the remaining identical singleton export belongs to a cavity-bearing region. No nontrivial cross-parameter embedding equivalence was demonstrated.

The frozen legacy filter is intentionally retained for audit reproducibility. It may hide unsupported cells and depends on label ordering. It must not be substituted for the revised optional filter or used in the main figure. The final revision manifest records the installed artifact hashes.
