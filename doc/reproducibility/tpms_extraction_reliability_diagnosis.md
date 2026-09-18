# TPMS extraction diagnosis and optional homology guard

## Independent staged experiment

I evaluated 135 parameter cells: three pairwise TPMS families, five values λ∈{0,.25,.5,.75,1}, three offsets c∈{0,.105,.3}, and resolutions 24³, 64³ and 96³. The radius/domain and trigonometric fields match the packaged TPMS implementation. Some family endpoints repeat the same field, so these are 135 evaluation cells rather than 135 unique geometries.

Each cell records the voxel mask's digital `(β0,β1,β2)`, the Lee skeleton's corresponding digital topology, the sparse adjacency graph counts, every junction radius 0–4, the persistence-selected graph, and leaf/simplification/smoothing stages. The data and source hashes are in `tpms_extraction_stage_diagnosis.json`; the executable script is `diagnose_tpms_extraction.py`.

| Resolution | Evaluated cells | Mask→skeleton Betti mismatches | Selected graph β0/β1 mismatches | Positive β2 cells |
|---|---:|---:|---:|---:|
| 24³ | 45 | 0 | 29 | 2 |
| 64³ | 45 | 0 | 0 | 2 |
| 96³ | 45 | 0 | 5 | 2 |

The zero-radius graph matched the mask's β0/β1 in all 135 cells. Therefore the observed loss of cycles is introduced by junction-zone expansion/persistence selection, not by Lee thinning. The sample also shows why simply increasing resolution cannot be assumed to cure the extraction problem.

Concrete 96³ failures:

- Gyroid endpoint λ=0, c=0: mask and skeleton β1=11; the selected graph has β1=8. Zero-radius β1=11 is rejected by a cleanup diagnostic, and hops1–4 repeatedly exhibit the wrong β1=8.
- Gyroid endpoint λ=0, c=.105: mask/skeleton β1=11; hop0 gives11, hop1 gives10, hops2–4 give8, and persistence selects8.
- Gyroid→Diamond λ=.25, c=.105: mask/skeleton β1=8; hop0 gives8 and later hops give5, which is selected.

Both gyroid endpoints appear in two interpolation families, accounting for the five 96³ failures. At24³, expansion often swallows the entire graph and the repeated one-vertex/no-edge result wins persistence despite losing all handles.

## Definite out-of-scope volume topology

Schwarz-P→Diamond at λ=.5 with c=.105 and c=.3 has `(β0,β1,β2)=(1,0,1)` at all three sampled resolutions. These are cavity-containing shell regions, not ordinary handlebodies with graph spines. The skeleton preserves the two-dimensional topology, but graphization collapses it to one vertex. A Betti1-only check would miss this error. The original compact TPMS treatment must not present these cells as ball-like phases. Keep them explicitly unsupported under the graph-spine model or develop the already-discussed boundary-resolved extension; do not fill the void merely to obtain a desired graph.

This also qualifies the earlier experimental contraction pair probe: its exact singleton-graph match at Schwarz-P→Diamond λ=.50/.55 concerned only the archived PL graph. It does not establish equivalence or correctness of the original finite-volume reconstructions. This is precisely why the graph-to-volume scope separation matters.

## Root cause in tracing and selection

`trace_component` collapses each connected expanded junction zone into a point, omitting all edges within the zone, and explicitly drops short loops at positive radius. Neither action checks whether a genuine cycle is erased. The selection rule rewards repeated graph isomorphism and local edge cleanliness, but repetition alone does not imply correct homology. Meanwhile, the abstract sparse-adjacency cycle rank is not the right target either: it can contain lattice clique cycles. For the96³ cavity skeleton, that graph's β1 is8173 while the digital object's β1=0 and β2=1. Use independently supplied volume/skeleton **digital** topology, not the raw adjacency E−V+C, when constraining extraction from a volume.

## Implemented conservative change

The isolated revision checkout now exposes optional keywords:

```python
skeleton_image_to_graph(
    skeleton,
    expected_components=topology.connected_components,
    expected_cycle_rank=topology.handle_rank,
)
```

The same options are available through `topology_aware_skeleton_image_to_graph`, `_optimized.extract`, `persistent_extract`, and its constrained compatibility wrapper. Candidate eligibility includes these counts. A zero-radius fallback may only return if it satisfies them; otherwise an explicit `ValueError` is raised. Both counts must be nonnegative integers. Unconstrained behavior and default call forwarding are unchanged. The change is confined to `extraction/skeleton.py`, `_optimized.py`, `_topology_optimized.py` and the new `tests/extraction/test_expected_skeleton_homology.py`, preserving the parent's prior isolated-node change.

The optional guard does not turn matching Betti numbers into a proof of embedding or handlebody equivalence. It prevents a demonstrated necessary-condition violation. The application must independently reject β2>0 and verify mask-to-skeleton topology; the guard cannot infer those properties from a graph.

## Verification after the change

The extraction tests and topology-aware skeleton tests pass: **26 passed**. They include the actual96³ gyroid regression, impossible expected counts, invalid argument types/ranges, both public aliases, and unchanged unconstrained behavior.

I reran the same135-cell staged experiment with supplied expected counts and explicit β2 rejection. The result is **129 cells processed with zero β0/β1 mismatches and zero postprocessing errors, plus six explicitly unsupported cavity cells**. Artifacts: `diagnose_tpms_extraction_constrained.py` and `tpms_extraction_constrained_diagnosis.json`.

A separate targeted check of all34 originally wrong-selected cells found that the zero-radius graph has nondegenerate edge geometry and survives leaf removal, degree-two simplification and ε=0 smoothing while retaining the target β0/β1 in every case. Those graphs have maximum valence3–5, so higher-valence invariant outputs still require the correct diagram-dependent labeling. This audit is `tpms_zero_radius_fallback_diagnosis.json`.

## Remaining limitations

Global β0/β1 equality does not ensure the right genus is assigned to each component, does not prove a spatial isotopy and does not certify that the mask is a3-manifold or regular neighborhood. Component-wise targets and explicit geometric contraction certificates would strengthen future reconstruction. The existing junction representative relocation and RDP smoothing likewise need geometric safety checks for a rigorous embedding claim. No such stronger result is claimed here, and the old phase-class counts have not been reinterpreted as certified classifications.
