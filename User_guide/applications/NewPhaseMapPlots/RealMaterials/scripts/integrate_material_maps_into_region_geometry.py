#!/usr/bin/env python3
"""Merge material-parameter phase maps into the existing region-geometry HTML."""

from __future__ import annotations

import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/codex-matplotlib")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import networkx as nx
import numpy as np
from skimage.measure import euler_number, label as label_volume

THIS_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(THIS_REPO / "tmp"))

from material_parameter_phase_maps import (  # noqa: E402
    K_SYMBOLS,
    MaterialFermiSurface,
    PHASE_COLORS,
    apply_signature_phase_merges,
    connected_components_for_label,
    material_families,
    remap_labels_contiguous,
    reset_gap_threshold,
    stable_labels,
)
from knotted_graph.extraction import skeletonize_volume  # noqa: E402


TARGET_HTML = Path(
    "/Users/hakanakgun/Desktop/Projects/ProfLeeProjects/Knotted_graph_code_paper/"
    "FigureGeneration/figures/07_hamiltonian_yamada_phase_maps/"
    "07_hamiltonian_yamada_plotly_region_geometry.html"
)
MATERIAL_DIR = THIS_REPO / "tmp" / "material_parameter_phase_maps"
RECORDS_JSON = MATERIAL_DIR / "material_parameter_phase_map_records.json"
SUMMARY_JSON = MATERIAL_DIR / "material_parameter_phase_map_summary.json"
OUT_HTML = MATERIAL_DIR / "07_hamiltonian_yamada_plotly_region_geometry_with_materials.html"

MAX_SURFACE_TRIANGLES = 1200
ROUND_DIGITS = 5
MIN_STABLE_CELLS = int(os.environ.get("MATERIAL_PHASE_MIN_ISLAND_CELLS", "8"))
ISLAND_MERGE_STRATEGY = os.environ.get("MATERIAL_PHASE_ISLAND_MERGE", "below").strip().lower()
PRIMARY_COMPONENT_MIN_FRACTION = float(os.environ.get("MATERIAL_PRIMARY_COMPONENT_MIN_FRACTION", "0.25"))


def rounded(value: Any, digits: int = ROUND_DIGITS) -> Any:
    if value is None:
        return None
    if isinstance(value, (np.floating, float)):
        return round(float(value), digits)
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def rounded_list(values: Any, digits: int = ROUND_DIGITS) -> list[Any]:
    arr = np.asarray(values)
    if arr.size == 0:
        return []
    flat = arr.reshape(-1)
    return [rounded(v, digits) for v in flat]


def parse_payload(html: str) -> tuple[dict[str, Any], str, str]:
    marker = "const payload = "
    start = html.index(marker) + len(marker)
    end = html.index(";\nconst transitions", start)
    return json.loads(html[start:end]), html[:start], html[end:]


def phase_grids(
    records: list[dict[str, Any]],
    lambdas: np.ndarray,
    energies: list[float],
) -> dict[str, Any]:
    signatures: list[str] = []
    seen: set[str] = set()
    for record in records:
        signature = str(record["phase_signature"])
        if signature in seen:
            continue
        seen.add(signature)
        signatures.append(signature)

    old_signature_to_id = {signature: idx + 1 for idx, signature in enumerate(signatures)}
    old_id_to_label = {
        idx + 1: next(
            str(record["phase_label"])
            for record in records
            if record["phase_signature"] == signature
        )
        for idx, signature in enumerate(signatures)
    }
    old_id_to_signature = {idx + 1: signature for idx, signature in enumerate(signatures)}

    lookup = {
        (round(float(record["energy"]), 12), round(float(record["lam"]), 12)): record
        for record in records
    }

    raw_old = np.zeros((len(energies), len(lambdas)), dtype=int)
    for row, energy in enumerate(energies):
        for col, lam in enumerate(lambdas):
            record = lookup[(round(float(energy), 12), round(float(lam), 12))]
            raw_old[row, col] = old_signature_to_id[str(record["phase_signature"])]

    stable_old, changed_cells = stable_labels(
        raw_old,
        min_cells=MIN_STABLE_CELLS,
        merge_strategy=ISLAND_MERGE_STRATEGY,
    )
    stable_old, manual_merges = apply_signature_phase_merges(
        stable_old, old_signature_to_id, str(records[0]["material"])
    )
    stable, stable_label_lookup = remap_labels_contiguous(stable_old, old_id_to_label)
    used_old = sorted(int(value) for value in set(stable_old.ravel()))
    old_to_new = {old: new for new, old in enumerate(used_old, start=1)}
    raw_remapped = np.vectorize(lambda value: old_to_new.get(int(value), 0))(raw_old).astype(int)
    stable_signature_lookup = {
        new: old_id_to_signature[old]
        for old, new in old_to_new.items()
    }

    return {
        "lookup": lookup,
        "raw_old": raw_old,
        "stable": stable,
        "raw_remapped": raw_remapped,
        "stable_label_lookup": stable_label_lookup,
        "stable_signature_lookup": stable_signature_lookup,
        "changed_cells": int(changed_cells),
        "manual_merges": manual_merges,
        "raw_classes": len(signatures),
    }


def choose_representative(
    component: list[tuple[int, int]],
    phase_id: int,
    raw_remapped: np.ndarray,
    lookup: dict[tuple[float, float], dict[str, Any]],
    lambdas: np.ndarray,
    energies: list[float],
) -> dict[str, Any]:
    center = np.asarray(component, dtype=float).mean(axis=0)

    def score(cell: tuple[int, int]) -> tuple[int, int, float, int, int]:
        row, col = cell
        same_raw = int(raw_remapped[row, col]) == int(phase_id)
        record = lookup[(round(float(energies[row]), 12), round(float(lambdas[col]), 12))]
        computed = bool(record.get("classification_computed", True))
        distance = float(np.linalg.norm(np.asarray([row, col], dtype=float) - center))
        return (0 if same_raw else 1, 0 if computed else 1, distance, row, col)

    row, col = min(component, key=score)
    return lookup[(round(float(energies[row]), 12), round(float(lambdas[col]), 12))]


def empty_surface() -> dict[str, Any]:
    return {"x": [], "y": [], "z": [], "i": [], "j": [], "k": [], "point_count": 0, "triangle_count": 0}


def boundary_faces_for_mask(mask: np.ndarray) -> list[str]:
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0 or not mask.any():
        return []
    faces = []
    for label, touched in (
        ("kx_min", mask[0, :, :].any()),
        ("kx_max", mask[-1, :, :].any()),
        ("ky_min", mask[:, 0, :].any()),
        ("ky_max", mask[:, -1, :].any()),
        ("kz_min", mask[:, :, 0].any()),
        ("kz_max", mask[:, :, -1].any()),
    ):
        if bool(touched):
            faces.append(label)
    return faces


def topology_for_mask(mask: np.ndarray) -> dict[str, Any]:
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0 or not mask.any():
        return {
            "interior_voxels": 0,
            "components": 0,
            "euler_characteristic": 0,
            "handle_rank": 0,
            "boundary_faces": [],
            "touches_boundary": False,
            "closed_in_window": True,
        }
    _, component_count = label_volume(mask, connectivity=3, return_num=True)
    euler = int(euler_number(mask, connectivity=3))
    faces = boundary_faces_for_mask(mask)
    return {
        "interior_voxels": int(mask.sum()),
        "components": int(component_count),
        "euler_characteristic": euler,
        "handle_rank": int(max(0, int(component_count) - euler)),
        "boundary_faces": faces,
        "touches_boundary": bool(faces),
        "closed_in_window": not bool(faces),
    }


def component_target_point(span: Any) -> np.ndarray:
    target = []
    for lo, hi in span:
        lo_f = float(lo)
        hi_f = float(hi)
        target.append(0.0 if lo_f <= 0.0 <= hi_f else 0.5 * (lo_f + hi_f))
    return np.asarray(target, dtype=float)


def selected_volume_component(obj: MaterialFermiSurface) -> dict[str, Any]:
    mask = np.asarray(obj._interior_mask, dtype=bool)
    total_voxels = int(mask.sum())
    if total_voxels == 0:
        return {
            "mask": mask,
            "component_id": 0,
            "component_count": 0,
            "kept_voxels": 0,
            "discarded_voxels": 0,
            "discarded_components": 0,
            "centroid": np.asarray([(lo + hi) * 0.5 for lo, hi in obj.span], dtype=float),
            "topology": topology_for_mask(mask),
        }

    labels, component_count = label_volume(mask, connectivity=3, return_num=True)
    counts = np.bincount(labels.ravel(), minlength=int(component_count) + 1)
    largest = int(counts[1:].max()) if component_count else total_voxels
    min_count = max(1, int(math.ceil(float(largest) * PRIMARY_COMPONENT_MIN_FRACTION)))
    target = component_target_point(obj.span)

    candidates = []
    for component_id in range(1, int(component_count) + 1):
        count = int(counts[component_id])
        if count < min_count:
            continue
        pts = np.argwhere(labels == component_id)
        centroid = np.asarray(obj._idx_to_coord(pts.mean(axis=0).reshape(1, 3))[0], dtype=float)
        distance = float(np.linalg.norm(centroid - target))
        candidates.append((distance, -count, component_id, centroid))
    if not candidates:
        for component_id in range(1, int(component_count) + 1):
            pts = np.argwhere(labels == component_id)
            centroid = np.asarray(obj._idx_to_coord(pts.mean(axis=0).reshape(1, 3))[0], dtype=float)
            distance = float(np.linalg.norm(centroid - target))
            candidates.append((distance, -int(counts[component_id]), component_id, centroid))

    _, neg_count, component_id, centroid = min(candidates)
    component_mask = labels == int(component_id)
    return {
        "mask": component_mask,
        "component_id": int(component_id),
        "component_count": int(component_count),
        "kept_voxels": int(-neg_count),
        "discarded_voxels": int(total_voxels - int(-neg_count)),
        "discarded_components": int(max(0, int(component_count) - 1)),
        "centroid": centroid,
        "topology": topology_for_mask(component_mask),
    }


def surface_component_near(mesh: Any, centroid: np.ndarray) -> Any:
    try:
        filtered = mesh.connectivity(
            extraction_mode="closest",
            closest_point=tuple(float(v) for v in centroid),
        )
        if filtered.n_points > 0 and filtered.n_cells > 0:
            return filtered
    except Exception as exc:
        print(f"    closest surface component extraction failed: {type(exc).__name__}: {exc}", flush=True)
    try:
        filtered = mesh.connectivity("largest")
        if filtered.n_points > 0 and filtered.n_cells > 0:
            return filtered
    except Exception as exc:
        print(f"    largest surface component extraction failed: {type(exc).__name__}: {exc}", flush=True)
    return mesh


def surface_payload(obj: MaterialFermiSurface, component_info: dict[str, Any] | None = None) -> dict[str, Any]:
    try:
        mesh = obj.exceptional_surface_pv
        if mesh.n_points == 0 or mesh.n_cells == 0:
            return empty_surface()
        mesh = mesh.extract_surface().triangulate().clean()
        if component_info and int(component_info.get("component_count", 1)) > 1:
            mesh = surface_component_near(mesh, np.asarray(component_info["centroid"], dtype=float))
            mesh = mesh.extract_surface().triangulate().clean()
        if mesh.n_cells > MAX_SURFACE_TRIANGLES:
            reduction = 1.0 - (MAX_SURFACE_TRIANGLES / float(mesh.n_cells))
            reduction = min(0.96, max(0.0, reduction))
            try:
                mesh = mesh.decimate_pro(reduction, preserve_topology=True).triangulate().clean()
            except Exception:
                mesh = mesh.decimate_pro(reduction, preserve_topology=False).triangulate().clean()
        points = np.asarray(mesh.points, dtype=float)
        faces = np.asarray(mesh.faces, dtype=int)
        if points.size == 0 or faces.size == 0:
            return empty_surface()
        faces = faces.reshape((-1, 4))
        faces = faces[faces[:, 0] == 3][:, 1:4]
        return {
            "x": rounded_list(points[:, 0]),
            "y": rounded_list(points[:, 1]),
            "z": rounded_list(points[:, 2]),
            "i": [int(v) for v in faces[:, 0]],
            "j": [int(v) for v in faces[:, 1]],
            "k": [int(v) for v in faces[:, 2]],
            "point_count": int(points.shape[0]),
            "triangle_count": int(faces.shape[0]),
        }
    except Exception as exc:
        print(f"    surface extraction failed: {type(exc).__name__}: {exc}", flush=True)
        return empty_surface()


def vertex_coordinate(obj: MaterialFermiSurface, component_info: dict[str, Any] | None = None) -> np.ndarray:
    try:
        mask = np.asarray(
            component_info["mask"] if component_info is not None else obj._interior_mask,
            dtype=bool,
        )
        if mask.any():
            center_idx = np.argwhere(mask).mean(axis=0)
            return np.asarray(obj._idx_to_coord(center_idx.reshape(1, 3))[0], dtype=float)
    except Exception:
        pass
    return np.asarray([(lo + hi) * 0.5 for lo, hi in obj.span], dtype=float)


def single_node_skeleton(coord: np.ndarray) -> dict[str, Any]:
    return {
        "lines": {"x": [], "y": [], "z": []},
        "nodes": {
            "x": rounded_list([coord[0]]),
            "y": rounded_list([coord[1]]),
            "z": rounded_list([coord[2]]),
        },
        "node_count": 1,
        "edge_count": 0,
        "components": 1,
        "cycle_rank": 0,
    }


def component_vertex_skeleton(
    obj: MaterialFermiSurface,
    component_info: dict[str, Any] | None = None,
) -> dict[str, Any]:
    try:
        mask = np.asarray(
            component_info["mask"] if component_info is not None else obj._interior_mask,
            dtype=bool,
        )
        labels, component_count = label_volume(mask, connectivity=3, return_num=True)
        coords = []
        for component_id in range(1, int(component_count) + 1):
            pts = np.argwhere(labels == component_id)
            if pts.size == 0:
                continue
            coords.append(obj._idx_to_coord(pts.mean(axis=0).reshape(1, 3))[0])
        if coords:
            arr = np.asarray(coords, dtype=float)
            return {
                "lines": {"x": [], "y": [], "z": []},
                "nodes": {
                    "x": rounded_list(arr[:, 0]),
                    "y": rounded_list(arr[:, 1]),
                    "z": rounded_list(arr[:, 2]),
                },
                "node_count": int(arr.shape[0]),
                "edge_count": 0,
                "components": int(arr.shape[0]),
                "cycle_rank": 0,
            }
    except Exception as exc:
        print(f"    component vertex extraction failed: {type(exc).__name__}: {exc}", flush=True)
    return single_node_skeleton(vertex_coordinate(obj, component_info))


def filtered_skeleton_image(
    obj: MaterialFermiSurface,
    component_info: dict[str, Any] | None,
) -> np.ndarray | None:
    if not component_info or int(component_info.get("component_count", 1)) <= 1:
        return None
    mask = np.asarray(component_info["mask"], dtype=bool)
    try:
        skeleton = np.asarray(obj._skeleton_image, dtype=bool) & mask
        if skeleton.any():
            return skeleton
    except Exception:
        pass
    try:
        skeleton = np.asarray(skeletonize_volume(mask), dtype=bool)
        if skeleton.any():
            return skeleton
    except Exception as exc:
        print(f"    component skeletonization failed: {type(exc).__name__}: {exc}", flush=True)
    return None


def skeleton_payload_from_graph(obj: MaterialFermiSurface, graph: nx.MultiGraph) -> dict[str, Any]:
    if graph.number_of_nodes() == 0:
        coord = vertex_coordinate(obj)
        return single_node_skeleton(coord)

    node_positions = np.asarray([data["pos"] for _, data in graph.nodes(data=True)], dtype=float)
    node_coords = np.asarray(obj._idx_to_coord(node_positions), dtype=float)
    line_x: list[Any] = []
    line_y: list[Any] = []
    line_z: list[Any] = []
    for u, v, _, data in graph.edges(keys=True, data=True):
        pts = np.asarray(data.get("pts", []), dtype=float)
        if pts.ndim != 2 or pts.shape[0] < 2 or pts.shape[1] != 3:
            pts = np.asarray([graph.nodes[u]["pos"], graph.nodes[v]["pos"]], dtype=float)
        coords = np.asarray(obj._idx_to_coord(pts), dtype=float)
        if coords.shape[0] > 100:
            keep = np.unique(np.linspace(0, coords.shape[0] - 1, 100).astype(int))
            coords = coords[keep]
        line_x.extend(rounded_list(coords[:, 0]))
        line_y.extend(rounded_list(coords[:, 1]))
        line_z.extend(rounded_list(coords[:, 2]))
        line_x.append(None)
        line_y.append(None)
        line_z.append(None)

    return {
        "lines": {"x": line_x, "y": line_y, "z": line_z},
        "nodes": {
            "x": rounded_list(node_coords[:, 0]),
            "y": rounded_list(node_coords[:, 1]),
            "z": rounded_list(node_coords[:, 2]),
        },
        "node_count": int(graph.number_of_nodes()),
        "edge_count": int(graph.number_of_edges()),
        "components": int(nx.number_connected_components(graph)) if graph.number_of_nodes() else 0,
        "cycle_rank": int(
            graph.number_of_edges()
            - graph.number_of_nodes()
            + (nx.number_connected_components(graph) if graph.number_of_nodes() else 0)
        ),
    }


def skeleton_payload(
    obj: MaterialFermiSurface,
    source: str,
    component_info: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], str]:
    if source == "vertex":
        skeleton = component_vertex_skeleton(obj, component_info)
        core_mode = "vertex-primary-component" if component_info else "vertex-after-prune"
        return (skeleton, core_mode)
    try:
        component_skeleton = filtered_skeleton_image(obj, component_info)
        graph_kwargs: dict[str, Any] = {}
        if component_skeleton is not None:
            graph_kwargs["skeleton_image"] = component_skeleton
        graph = obj.skeleton_graph(
            smooth_epsilon=0,
            simplify=True,
            force_small_edge_contraction=True,
            small_edge_limit=math.pi * 0.1,
            previous_n_edgepoint=20,
            **graph_kwargs,
        )
        suffix = "-primary-component" if component_skeleton is not None else ""
        return skeleton_payload_from_graph(obj, graph), f"material-smooth-0{suffix}"
    except Exception as exc:
        print(f"    skeleton extraction failed: {type(exc).__name__}: {exc}", flush=True)
        coord = vertex_coordinate(obj, component_info)
        return (single_node_skeleton(coord), "skeleton-unavailable")


def geometry_for_record(
    family: Any,
    record: dict[str, Any],
    object_cache: dict[float, MaterialFermiSurface],
) -> tuple[dict[str, Any], dict[str, Any], str, dict[str, Any]]:
    if str(record["source"]) == "boundary-open":
        coord = np.asarray([(lo + hi) * 0.5 for lo, hi in family.span], dtype=float)
        return empty_surface(), single_node_skeleton(coord), "boundary-open-not-skeletonized", {
            "component_id": 0,
            "component_count": int(record.get("interior_components", 0)),
            "kept_voxels": int(record.get("interior_voxels", 0)),
            "discarded_voxels": 0,
            "discarded_components": 0,
            "topology": {
                "interior_voxels": int(record.get("interior_voxels", 0)),
                "components": int(record.get("interior_components", 0)),
                "euler_characteristic": int(record.get("euler_characteristic", 0)),
                "handle_rank": int(record.get("handle_rank", 0)),
                "boundary_faces": list(record.get("boundary_faces", [])),
                "touches_boundary": bool(record.get("touches_boundary", False)),
                "closed_in_window": not bool(record.get("touches_boundary", False)),
            },
        }

    lam = round(float(record["lam"]), 12)
    if lam not in object_cache:
        print(
            f"  building volume {family.key} lambda={float(record['lam']):.4f} "
            f"{family.lambda_param}={family.param_value_at(float(record['lam'])):.6g}",
            flush=True,
        )
        object_cache[lam] = MaterialFermiSurface(
            family.hamiltonian_at(float(record["lam"])),
            k_symbols=K_SYMBOLS,
            span=family.span,
            dimension=family.dimension,
            band_pair=family.band_pair,
            gap_tol=float(record["energy"]),
            check_hermitian=False,
            check_pt_symmetry=False,
            chunk_size=100_000,
            force_small_edge_contraction=True,
            small_edge_limit=math.pi * 0.1,
            previous_n_edgepoint=20,
        )
    obj = object_cache[lam]
    reset_gap_threshold(obj, float(record["energy"]))
    component_info = selected_volume_component(obj)
    surface = surface_payload(obj, component_info)
    skeleton, core_mode = skeleton_payload(obj, str(record["source"]), component_info)
    return surface, skeleton, core_mode, component_info


def span_label(span: list[list[float]] | tuple[tuple[float, float], ...]) -> str:
    labels = []
    for axis, (lo, hi) in zip(("kx", "ky", "kz"), span):
        labels.append(f"{axis}=[{float(lo):.4g},{float(hi):.4g}]")
    return "; ".join(labels)


def attachment_component_payload(component_info: dict[str, Any]) -> dict[str, Any]:
    topology = dict(component_info.get("topology", {}))
    topology["boundary_faces"] = list(topology.get("boundary_faces", []))
    return {
        "selection": "primary connected volume component",
        "component_id": int(component_info.get("component_id", 0)),
        "full_component_count": int(component_info.get("component_count", 0)),
        "kept_voxels": int(component_info.get("kept_voxels", 0)),
        "discarded_components": int(component_info.get("discarded_components", 0)),
        "discarded_voxels": int(component_info.get("discarded_voxels", 0)),
        "topology": topology,
    }


def build_material_transition(
    payload: dict[str, Any],
    family_summary: dict[str, Any],
    family: Any,
    records: list[dict[str, Any]],
    lambdas: np.ndarray,
) -> dict[str, Any]:
    energies = [float(value) for value in family_summary["energies"]]
    grids = phase_grids(records, lambdas, energies)
    stable = grids["stable"]
    raw_remapped = grids["raw_remapped"]
    phase_ids = sorted(int(value) for value in set(stable.ravel()))
    colors = [PHASE_COLORS[(phase_id - 1) % len(PHASE_COLORS)] for phase_id in phase_ids]

    base_region_keys = [["" for _ in range(stable.shape[1])] for _ in range(stable.shape[0])]
    material_regions: dict[str, dict[str, Any]] = {}
    object_cache: dict[float, MaterialFermiSurface] = {}
    region_id = 0

    for phase_id in phase_ids:
        components = connected_components_for_label(stable, phase_id)
        components.sort(key=lambda comp: (min(row for row, _ in comp), min(col for _, col in comp)))
        for component in components:
            region_id += 1
            representative = choose_representative(
                component,
                phase_id,
                raw_remapped,
                grids["lookup"],
                lambdas,
                energies,
            )
            print(
                f"  representative {family.key} region {region_id}: phase={phase_id}, "
                f"lambda={float(representative['lam']):.4f}, E={float(representative['energy']):.6g}, "
                f"source={representative['source']}, cells={len(component)}",
                flush=True,
            )
            surface, skeleton, core_mode, component_info = geometry_for_record(family, representative, object_cache)
            attachment = attachment_component_payload(component_info)
            topology = dict(attachment["topology"])
            color = PHASE_COLORS[(phase_id - 1) % len(PHASE_COLORS)]
            base = {
                "transition": family.key,
                "title": family.title,
                "regionId": int(region_id),
                "phaseId": int(phase_id),
                "stableRegionId": int(region_id),
                "stablePhaseId": int(phase_id),
                "stableSignature": grids["stable_signature_lookup"].get(int(phase_id), str(representative["phase_signature"])),
                "phaseLabel": grids["stable_label_lookup"].get(int(phase_id), str(representative["phase_label"])),
                "lambda": rounded(float(representative["lam"]), 6),
                "Gamma": rounded(float(representative["energy"]), 6),
                "parameterName": "E",
                "parameterLabel": "E",
                "parameterValue": rounded(float(representative["energy"]), 6),
                "lambdaParameter": family.lambda_param,
                "materialParameterValue": rounded(float(representative["parameter_value"]), 6),
                "dimension": int(family.dimension),
                "spanLabel": span_label(family.span),
                "bandPair": list(family.band_pair),
                "cellCount": int(len(component)),
                "source": str(representative["source"]),
                "coreMode": core_mode,
                "nodes": int(skeleton.get("node_count", representative["nodes"])),
                "edges": int(skeleton.get("edge_count", representative["edges"])),
                "components": int(skeleton.get("components", representative["components"])),
                "cycleRank": int(skeleton.get("cycle_rank", representative["cycle_rank"])),
                "polynomial": str(representative["polynomial"] or representative["phase_label"]),
                "classicRegionIds": [int(region_id)],
                "classicPhaseIds": [int(phase_id)],
                "startName": f"{family.lambda_param}={family.param_start:.6g}",
                "endName": f"{family.lambda_param}={family.param_end:.6g}",
                "color": color,
                "surface": surface,
                "skeleton": skeleton,
                "attachmentComponent": attachment,
                "fullVolumeTopology": {
                    "interior_voxels": int(representative["interior_voxels"]),
                    "components": int(representative["interior_components"]),
                    "euler_characteristic": int(representative["euler_characteristic"]),
                    "handle_rank": int(representative["handle_rank"]),
                    "boundary_faces": list(representative["boundary_faces"]),
                    "touches_boundary": bool(representative["touches_boundary"]),
                    "closed_in_window": not bool(representative["touches_boundary"]),
                    "forces_vertex": str(representative["source"]) == "vertex",
                },
                "topology": {
                    "interior_voxels": int(topology.get("interior_voxels", 0)),
                    "components": int(topology.get("components", 0)),
                    "euler_characteristic": int(topology.get("euler_characteristic", 0)),
                    "handle_rank": int(topology.get("handle_rank", 0)),
                    "boundary_faces": list(topology.get("boundary_faces", [])),
                    "touches_boundary": bool(topology.get("touches_boundary", False)),
                    "closed_in_window": bool(topology.get("closed_in_window", True)),
                    "forces_vertex": str(representative["source"]) == "vertex",
                },
            }
            for row, col in component:
                base_region_keys[row][col] = f"classic:{family.key}:{region_id}"
            for mode_key, mode_title in (
                ("classic", "Classic"),
                ("contraction", "Up to contraction moves"),
            ):
                region = dict(base)
                region["mode"] = mode_key
                region["modeTitle"] = mode_title if mode_key == "classic" else "Material scan in contraction view"
                region["key"] = f"{mode_key}:{family.key}:{region_id}"
                material_regions[region["key"]] = region

    modes = {}
    for mode_key in ("classic", "contraction"):
        region_keys = [
            [key.replace("classic:", f"{mode_key}:") for key in row]
            for row in base_region_keys
        ]
        modes[mode_key] = {
            "z": stable.tolist(),
            "regionKeys": region_keys,
            "colors": colors,
            "classes": int(len(phase_ids)),
            "components": int(region_id),
            "rawClasses": int(grids["raw_classes"]),
            "classicStableComponents": int(region_id),
            "smallComponents": int(grids["changed_cells"]),
            "minIslandCells": int(MIN_STABLE_CELLS),
            "islandMergeStrategy": str(ISLAND_MERGE_STRATEGY),
            "manualPhaseMerges": grids["manual_merges"],
            "contractionAccessiblePairs": 0,
            "contractionMergedClassicRegions": int(region_id),
        }

    transition = {
        "key": family.key,
        "title": family.title,
        "displayTitle": (
            f"{family.title}: {family.lambda_param}(λ)="
            f"{family.param_start:.6g}→{family.param_end:.6g}"
        ),
        "lambdas": [rounded(float(value), 6) for value in lambdas],
        "gammas": [rounded(float(value), 6) for value in energies],
        "parameters": [rounded(float(value), 6) for value in energies],
        "parameterName": "E",
        "parameterLabel": "E",
        "terminalGamma": rounded(float(energies[-1]), 6),
        "startName": f"{family.lambda_param}={family.param_start:.6g}",
        "endName": f"{family.lambda_param}={family.param_end:.6g}",
        "lambdaSamples": int(len(lambdas)),
        "materialParameter": family.lambda_param,
        "materialParameterStart": rounded(family.param_start, 6),
        "materialParameterEnd": rounded(family.param_end, 6),
        "dimension": int(family.dimension),
        "bandPair": list(family.band_pair),
        "spanLabel": span_label(family.span),
        "postProcessing": {
            "minIslandCells": int(MIN_STABLE_CELLS),
            "islandMergeStrategy": str(ISLAND_MERGE_STRATEGY),
            "islandRule": "connected lambda-E components smaller than the cutoff are relabeled to the closest large phase below in E, falling back to nearest large phase if no lower phase exists",
            "manualPhaseMerges": grids["manual_merges"],
            "geometryRule": "click attachments keep only the selected connected volume component; disconnected pieces are ignored unless merged in the mask",
            "primaryComponentMinFraction": float(PRIMARY_COMPONENT_MIN_FRACTION),
        },
        "modes": modes,
    }
    payload["regions"].update(material_regions)
    return transition


def patch_controller(html: str) -> str:
    html = html.replace(
        "function htmlTransitionTitle(item, mode) {\n"
        "  const startName = item.startName || \"0\";\n"
        "  const endName = item.endName || \"1\";\n"
        "  const suffix = mode.key === \"classic\" ? \"\" : ` (${mode.short_title || mode.title})`;\n"
        "  return `(1−λ)<i>H</i><sub>${startName}</sub>+λ<i>H</i><sub>${endName}</sub>${suffix}`;\n"
        "}",
        "function htmlTransitionTitle(item, mode) {\n"
        "  const startName = item.startName || \"0\";\n"
        "  const endName = item.endName || \"1\";\n"
        "  const suffix = mode.key === \"classic\" ? \"\" : ` (${mode.short_title || mode.title})`;\n"
        "  if (item.displayTitle) return `${item.displayTitle}${suffix}`;\n"
        "  return `(1−λ)<i>H</i><sub>${startName}</sub>+λ<i>H</i><sub>${endName}</sub>${suffix}`;\n"
        "}",
    )
    html = html.replace(
        "  const skeletonKind = region.skeleton.edge_count ? \"edge skeleton\" : \"vertex-only skeleton\";\n"
        "  const parameterLabel = \"E\";\n"
        "  const parameterValue = region.parameterValue ?? region.Gamma;\n"
        "  Plotly.react(\"geometry\", geometryTraces(region), layout, { responsive: true, displaylogo: false });\n"
        "  document.getElementById(\"status\").innerHTML = `\n"
        "    <strong>${region.title} - ${region.modeTitle}</strong>\n"
        "    <span>region ${region.regionId}, phase ${region.phaseId}; representative λ=${region.lambda}, ${parameterLabel}=${parameterValue}; ${region.source}, ${region.coreMode}; ${skeletonKind}; core V=${region.nodes}, E=${region.edges}, β=${region.cycleRank}; surface triangles=${region.surface.triangle_count}</span>\n"
        "    <span>${region.startName ? `endpoints: ${region.startName} → ${region.endName}` : \"\"}</span>\n"
        "    <span>${topologyStatusLine(region.topology)}</span>\n"
        "    <span>Classic regions: ${compactIdList(region.classicRegionIds)}; Classic phases: ${compactIdList(region.classicPhaseIds)}; cells=${region.cellCount}</span>\n"
        "    <span class=\"poly\">Yamada: ${yamadaDisplay(region.polynomial)}</span>\n"
        "  `;\n",
        "  const skeletonKind = region.skeleton.edge_count ? \"edge skeleton\" : \"vertex-only skeleton\";\n"
        "  const parameterLabel = \"E\";\n"
        "  const parameterValue = region.parameterValue ?? region.Gamma;\n"
        "  const materialLine = region.lambdaParameter ? `<span>${region.lambdaParameter}(λ)=${region.materialParameterValue}; band pair=${compactIdList(region.bandPair)}; N=${region.dimension}; ${region.spanLabel}</span>` : \"\";\n"
        "  const polynomialLabel = region.source === \"large-core\" ? \"Phase label\" : \"Yamada\";\n"
        "  Plotly.react(\"geometry\", geometryTraces(region), layout, { responsive: true, displaylogo: false });\n"
        "  document.getElementById(\"status\").innerHTML = `\n"
        "    <strong>${region.title} - ${region.modeTitle}</strong>\n"
        "    <span>region ${region.regionId}, phase ${region.phaseId}; representative λ=${region.lambda}, ${parameterLabel}=${parameterValue}; ${region.source}, ${region.coreMode}; ${skeletonKind}; core V=${region.nodes}, E=${region.edges}, β=${region.cycleRank}; surface triangles=${region.surface.triangle_count}</span>\n"
        "    <span>${region.startName ? `endpoints: ${region.startName} → ${region.endName}` : \"\"}</span>\n"
        "    ${materialLine}\n"
        "    <span>${topologyStatusLine(region.topology)}</span>\n"
        "    <span>Classic regions: ${compactIdList(region.classicRegionIds)}; Classic phases: ${compactIdList(region.classicPhaseIds)}; cells=${region.cellCount}</span>\n"
        "    <span class=\"poly\">${polynomialLabel}: ${yamadaDisplay(region.polynomial || region.phaseLabel || region.stableSignature)}</span>\n"
        "  `;\n",
    )
    if "const attachmentLine =" not in html:
        html = html.replace(
            "  const materialLine = region.lambdaParameter ? `<span>${region.lambdaParameter}(λ)=${region.materialParameterValue}; band pair=${compactIdList(region.bandPair)}; N=${region.dimension}; ${region.spanLabel}</span>` : \"\";\n"
            "  const polynomialLabel = region.source === \"large-core\" ? \"Phase label\" : \"Yamada\";\n",
            "  const materialLine = region.lambdaParameter ? `<span>${region.lambdaParameter}(λ)=${region.materialParameterValue}; band pair=${compactIdList(region.bandPair)}; N=${region.dimension}; ${region.spanLabel}</span>` : \"\";\n"
            "  const attachment = region.attachmentComponent;\n"
            "  const attachmentLine = attachment && attachment.full_component_count > 1 ? `<span>attached component: kept 1 of ${attachment.full_component_count} volume components; ignored ${attachment.discarded_components} disconnected components (${attachment.discarded_voxels} voxels)</span>` : \"\";\n"
            "  const polynomialLabel = region.source === \"large-core\" ? \"Phase label\" : \"Yamada\";\n",
        )
        html = html.replace(
            "    ${materialLine}\n"
            "    <span>${topologyStatusLine(region.topology)}</span>\n",
            "    ${materialLine}\n"
            "    ${attachmentLine}\n"
            "    <span>${topologyStatusLine(region.topology)}</span>\n",
        )
    return html


def main() -> None:
    started = time.perf_counter()
    html = TARGET_HTML.read_text(encoding="utf-8")
    payload, prefix, suffix = parse_payload(html)
    records = json.loads(RECORDS_JSON.read_text(encoding="utf-8"))
    summary = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
    lambdas = np.asarray(summary["lambdas"], dtype=float)

    material_keys = {str(item["key"]) for item in summary["families"]}
    all_known_material_keys = {
        family.key
        for dim in {80, *[int(item["dimension"]) for item in summary["families"]]}
        for family in material_families(int(dim))
    }
    payload["transitions"] = [
        item for item in payload["transitions"]
        if str(item.get("key")) not in all_known_material_keys
    ]
    payload["regions"] = {
        key: region for key, region in payload["regions"].items()
        if str(region.get("transition")) not in all_known_material_keys
    }
    for key in all_known_material_keys:
        payload.get("stable_min_component_cells", {}).pop(key, None)

    existing_keys = {item["key"] for item in payload["transitions"]}
    payload["version"] = "hamiltonian_yamada_plotly_region_geometry_v18_manual_phase_merges"
    payload.setdefault("material_parameter_scan", {
        "records": str(RECORDS_JSON),
        "summary": str(SUMMARY_JSON),
        "note": "Material Hamiltonian coefficient scans appended from the current material_parameter_phase_maps run.",
    })

    records_by_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        records_by_family[str(record["material"])].append(record)

    families_by_dimension = {
        int(dim): {family.key: family for family in material_families(int(dim))}
        for dim in sorted({int(item["dimension"]) for item in summary["families"]})
    }

    new_transitions = []
    for family_summary in summary["families"]:
        key = str(family_summary["key"])
        if key in existing_keys:
            continue
        family = families_by_dimension[int(family_summary["dimension"])][key]
        print(f"\n[{key}] adding {family.title}", flush=True)
        transition = build_material_transition(
            payload,
            family_summary,
            family,
            records_by_family[key],
            lambdas,
        )
        new_transitions.append(transition)

    payload["transitions"].extend(new_transitions)
    payload["stable_min_component_cells"].update(
        {item["key"]: MIN_STABLE_CELLS for item in new_transitions}
    )
    payload["material_parameter_scan"]["family_count"] = len(new_transitions)
    payload["material_parameter_scan"]["source_counts"] = dict(Counter(record["source"] for record in records))
    payload["material_parameter_scan"]["post_processing"] = {
        "min_island_cells": int(MIN_STABLE_CELLS),
        "island_merge_strategy": str(ISLAND_MERGE_STRATEGY),
        "primary_component_min_fraction": float(PRIMARY_COMPONENT_MIN_FRACTION),
        "island_rule": "small connected lambda-E islands are immersed into the spatially closest large-phase cell below in E, with vertical separation used only as a tie-breaker and nearest-large fallback only when no lower phase exists",
        "manual_phase_merges": {
            item["key"]: item["modes"]["classic"].get("manualPhaseMerges", [])
            for item in new_transitions
        },
        "geometry_rule": "surface/skeleton click attachments are restricted to the selected connected volume component",
    }

    merged = prefix + json.dumps(payload, separators=(",", ":"), ensure_ascii=False) + suffix
    merged = patch_controller(merged)
    OUT_HTML.write_text(merged, encoding="utf-8")
    print(f"\nwrote: {OUT_HTML}", flush=True)
    print(f"new transitions: {[item['key'] for item in new_transitions]}", flush=True)
    print(f"total transitions: {len(payload['transitions'])}", flush=True)
    print(f"total regions: {len(payload['regions'])}", flush=True)
    print(f"elapsed: {time.perf_counter() - started:.1f}s", flush=True)


if __name__ == "__main__":
    main()
