#!/usr/bin/env python3
"""Merge material-parameter phase maps into the existing region-geometry HTML."""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pyvista as pv
from skimage.measure import euler_number, label as label_volume, marching_cubes

from knotted_graph.applications.phase_map_examples._materials import (
    K_SYMBOLS,
    MaterialFermiSurface,
    PHASE_COLORS,
    apply_signature_phase_merges,
    attach_c6_invalid_components_to_lower_phase,
    boundary_filling_groups,
    enclosed_void_masks,
    resolve_mask_components,
    connected_components_for_label,
    material_families,
    remap_labels_contiguous,
    reset_gap_threshold,
    stable_labels,
)
from knotted_graph.extraction import skeletonize_volume  # noqa: E402


THIS_REPO = Path(__file__).resolve().parents[1]
TARGET_HTML = (
    Path(__file__).resolve().parents[5]
    / "doc"
    / "assets"
    / "demos"
    / "hamiltonian_yamada_phase_map.html"
)
MATERIAL_DIR = THIS_REPO / "data"

MAX_SURFACE_TRIANGLES = 1200
ROUND_DIGITS = 5
MIN_STABLE_CELLS = 8
ISLAND_MERGE_STRATEGY = "below"


C6_REPRESENTATIVE_HINTS = {
    (
        "tib2_d6_F",
        "yamada-set:outer[yamada:Integer(0)]inner[yamada:Integer(-1)]",
    ): (3.883729, 3.78),
}


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

    old_signature_to_id = {
        signature: idx + 1 for idx, signature in enumerate(signatures)
    }
    old_id_to_label = {
        idx + 1: next(
            str(record["phase_label"])
            for record in records
            if record["phase_signature"] == signature
        )
        for idx, signature in enumerate(signatures)
    }
    old_id_to_signature = {
        idx + 1: signature for idx, signature in enumerate(signatures)
    }

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
    stable_old, c6_symmetry_merges = attach_c6_invalid_components_to_lower_phase(
        stable_old, old_signature_to_id, str(records[0]["material"])
    )
    stable, stable_label_lookup = remap_labels_contiguous(stable_old, old_id_to_label)
    used_old = sorted(int(value) for value in set(stable_old.ravel()))
    old_to_new = {old: new for new, old in enumerate(used_old, start=1)}
    raw_remapped = np.vectorize(lambda value: old_to_new.get(int(value), 0))(
        raw_old
    ).astype(int)
    stable_signature_lookup = {
        new: old_id_to_signature[old] for old, new in old_to_new.items()
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
        "c6_symmetry_merges": c6_symmetry_merges,
        "raw_classes": len(signatures),
    }


def choose_representative(
    component: list[tuple[int, int]],
    phase_id: int,
    raw_remapped: np.ndarray,
    lookup: dict[tuple[float, float], dict[str, Any]],
    lambdas: np.ndarray,
    energies: list[float],
    preferred_parameter_energy: tuple[float, float] | None = None,
) -> dict[str, Any]:
    center = np.asarray(component, dtype=float).mean(axis=0)

    def score(cell: tuple[int, int]) -> tuple[int, int, float, float, int, int]:
        row, col = cell
        same_raw = int(raw_remapped[row, col]) == int(phase_id)
        record = lookup[
            (round(float(energies[row]), 12), round(float(lambdas[col]), 12))
        ]
        computed = bool(record.get("classification_computed", True))
        preferred_distance = 0.0
        if preferred_parameter_energy is not None:
            parameter_value, energy = preferred_parameter_energy
            preferred_distance = abs(
                float(record["parameter_value"]) - float(parameter_value)
            ) + abs(float(record["energy"]) - float(energy))
        distance = float(np.linalg.norm(np.asarray([row, col], dtype=float) - center))
        return (
            0 if same_raw else 1,
            0 if computed else 1,
            preferred_distance,
            distance,
            row,
            col,
        )

    row, col = min(component, key=score)
    return lookup[(round(float(energies[row]), 12), round(float(lambdas[col]), 12))]


def empty_surface() -> dict[str, Any]:
    return {
        "x": [],
        "y": [],
        "z": [],
        "i": [],
        "j": [],
        "k": [],
        "point_count": 0,
        "triangle_count": 0,
    }


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
            "void_components": 0,
            "boundary_components": 0,
            "boundary_faces": [],
            "touches_boundary": False,
            "closed_in_window": True,
        }
    _, component_count = label_volume(mask, connectivity=3, return_num=True)
    void_count = len(enclosed_void_masks(mask))
    euler = int(euler_number(mask, connectivity=3))
    faces = boundary_faces_for_mask(mask)
    return {
        "interior_voxels": int(mask.sum()),
        "components": int(component_count),
        "euler_characteristic": euler,
        "handle_rank": int(max(0, int(component_count) + int(void_count) - euler)),
        "void_components": int(void_count),
        "boundary_components": int(component_count) + int(void_count),
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
    raw_mask = np.asarray(obj._interior_mask, dtype=bool)
    total_voxels = int(raw_mask.sum())
    if total_voxels == 0:
        return {
            "mask": raw_mask,
            "component_id": 0,
            "component_count": 0,
            "kept_voxels": 0,
            "discarded_voxels": 0,
            "discarded_components": 0,
            "centroid": np.asarray(
                [(lo + hi) * 0.5 for lo, hi in obj.span], dtype=float
            ),
            "topology": topology_for_mask(raw_mask),
        }

    _, raw_component_count = label_volume(raw_mask, connectivity=3, return_num=True)
    component_mask, resolution = resolve_mask_components(raw_mask)
    kept_voxels = int(component_mask.sum())
    center_idx = np.argwhere(component_mask).mean(axis=0)
    centroid = np.asarray(obj._idx_to_coord(center_idx.reshape(1, 3))[0], dtype=float)
    return {
        "mask": component_mask,
        "component_id": 1,
        "component_count": int(raw_component_count),
        "kept_voxels": kept_voxels,
        "discarded_voxels": int(resolution["removed_component_voxels"]),
        "discarded_components": int(max(0, int(raw_component_count) - 1)),
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
        print(
            f"    closest surface component extraction failed: {type(exc).__name__}: {exc}",
            flush=True,
        )
    try:
        filtered = mesh.connectivity("largest")
        if filtered.n_points > 0 and filtered.n_cells > 0:
            return filtered
    except Exception as exc:
        print(
            f"    largest surface component extraction failed: {type(exc).__name__}: {exc}",
            flush=True,
        )
    return mesh


def surface_payload(
    obj: MaterialFermiSurface, component_info: dict[str, Any] | None = None
) -> dict[str, Any]:
    try:
        mask = np.asarray(
            component_info["mask"]
            if component_info is not None
            else obj._interior_mask,
            dtype=bool,
        )
        if not mask.any() or mask.all():
            return empty_surface()
        vertices, faces, _, _ = marching_cubes(mask.astype(np.uint8), level=0.5)
        points = np.asarray(obj._idx_to_coord(vertices), dtype=float)
        pv_faces = np.column_stack(
            [np.full(len(faces), 3, dtype=np.int64), np.asarray(faces, dtype=np.int64)]
        ).ravel()
        mesh = (
            pv.PolyData(points, pv_faces)
            .extract_surface(algorithm="dataset_surface")
            .triangulate()
            .clean()
        )
        if mesh.n_cells > MAX_SURFACE_TRIANGLES:
            reduction = 1.0 - (MAX_SURFACE_TRIANGLES / float(mesh.n_cells))
            reduction = min(0.96, max(0.0, reduction))
            try:
                mesh = (
                    mesh.decimate_pro(reduction, preserve_topology=True)
                    .triangulate()
                    .clean()
                )
            except Exception:
                mesh = (
                    mesh.decimate_pro(reduction, preserve_topology=False)
                    .triangulate()
                    .clean()
                )
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


def vertex_coordinate(
    obj: MaterialFermiSurface, component_info: dict[str, Any] | None = None
) -> np.ndarray:
    try:
        mask = np.asarray(
            component_info["mask"]
            if component_info is not None
            else obj._interior_mask,
            dtype=bool,
        )
        if mask.any():
            center_idx = np.argwhere(mask).mean(axis=0)
            return np.asarray(
                obj._idx_to_coord(center_idx.reshape(1, 3))[0], dtype=float
            )
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
            component_info["mask"]
            if component_info is not None
            else obj._interior_mask,
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
        print(
            f"    component vertex extraction failed: {type(exc).__name__}: {exc}",
            flush=True,
        )
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
        print(
            f"    component skeletonization failed: {type(exc).__name__}: {exc}",
            flush=True,
        )
    return None


def skeleton_payload_from_graph(
    obj: MaterialFermiSurface, graph: nx.MultiGraph
) -> dict[str, Any]:
    if graph.number_of_nodes() == 0:
        coord = vertex_coordinate(obj)
        return single_node_skeleton(coord)

    node_positions = np.asarray(
        [data["pos"] for _, data in graph.nodes(data=True)], dtype=float
    )
    node_coords = np.asarray(obj._idx_to_coord(node_positions), dtype=float)
    line_x: list[Any] = []
    line_y: list[Any] = []
    line_z: list[Any] = []
    for u, v, _, data in graph.edges(keys=True, data=True):
        pts = np.asarray(data.get("pts", []), dtype=float)
        if pts.ndim != 2 or pts.shape[0] < 2 or pts.shape[1] != 3:
            pts = np.asarray(
                [graph.nodes[u]["pos"], graph.nodes[v]["pos"]], dtype=float
            )
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
        "components": int(nx.number_connected_components(graph))
        if graph.number_of_nodes()
        else 0,
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
        core_mode = (
            "vertex-primary-component" if component_info else "vertex-after-prune"
        )
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
        print(
            f"    skeleton extraction failed: {type(exc).__name__}: {exc}", flush=True
        )
        coord = vertex_coordinate(obj, component_info)
        return (single_node_skeleton(coord), "skeleton-unavailable")


def geometry_for_record(
    family: Any,
    record: dict[str, Any],
    object_cache: dict[float, MaterialFermiSurface],
) -> tuple[dict[str, Any], dict[str, Any], str, dict[str, Any]]:
    if str(record["source"]) == "boundary-open":
        coord = np.asarray([(lo + hi) * 0.5 for lo, hi in family.span], dtype=float)
        return (
            empty_surface(),
            single_node_skeleton(coord),
            "boundary-open-not-skeletonized",
            {
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
            },
        )

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
    skeleton, core_mode = boundary_resolved_skeleton_payload(obj, component_info)
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
    colors = [
        PHASE_COLORS[(phase_id - 1) % len(PHASE_COLORS)] for phase_id in phase_ids
    ]

    base_region_keys = [
        ["" for _ in range(stable.shape[1])] for _ in range(stable.shape[0])
    ]
    material_regions: dict[str, dict[str, Any]] = {}
    object_cache: dict[float, MaterialFermiSurface] = {}
    region_id = 0

    for phase_id in phase_ids:
        components = connected_components_for_label(stable, phase_id)
        components.sort(
            key=lambda comp: (min(row for row, _ in comp), min(col for _, col in comp))
        )
        for component in components:
            region_id += 1
            stable_signature = grids["stable_signature_lookup"].get(int(phase_id), "")
            representative = choose_representative(
                component,
                phase_id,
                raw_remapped,
                grids["lookup"],
                lambdas,
                energies,
                C6_REPRESENTATIVE_HINTS.get((family.key, stable_signature)),
            )
            print(
                f"  representative {family.key} region {region_id}: phase={phase_id}, "
                f"lambda={float(representative['lam']):.4f}, E={float(representative['energy']):.6g}, "
                f"source={representative['source']}, cells={len(component)}",
                flush=True,
            )
            surface, skeleton, core_mode, component_info = geometry_for_record(
                family, representative, object_cache
            )
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
                "stableSignature": grids["stable_signature_lookup"].get(
                    int(phase_id), str(representative["phase_signature"])
                ),
                "phaseLabel": grids["stable_label_lookup"].get(
                    int(phase_id), str(representative["phase_label"])
                ),
                "lambda": rounded(float(representative["lam"]), 6),
                "Gamma": rounded(float(representative["energy"]), 6),
                "parameterName": "E",
                "parameterLabel": "E",
                "parameterValue": rounded(float(representative["energy"]), 6),
                "lambdaParameter": family.lambda_param,
                "materialParameterValue": rounded(
                    float(representative["parameter_value"]), 6
                ),
                "dimension": int(family.dimension),
                "spanLabel": span_label(family.span),
                "bandPair": list(family.band_pair),
                "cellCount": int(len(component)),
                "source": str(representative["source"]),
                "coreMode": core_mode,
                "nodes": int(skeleton.get("node_count", representative["nodes"])),
                "edges": int(skeleton.get("edge_count", representative["edges"])),
                "components": int(
                    skeleton.get("components", representative["components"])
                ),
                "cycleRank": int(
                    skeleton.get("cycle_rank", representative["cycle_rank"])
                ),
                "polynomial": str(
                    representative["polynomial"] or representative["phase_label"]
                ),
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
                    "euler_characteristic": int(
                        topology.get("euler_characteristic", 0)
                    ),
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
                region["modeTitle"] = (
                    mode_title
                    if mode_key == "classic"
                    else "Material scan in contraction view"
                )
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
            "c6SymmetryMerges": grids["c6_symmetry_merges"],
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
            "islandRule": "connected lambda-E components at or below the cutoff are relabeled to the closest large phase below in E, falling back to nearest large phase if no lower phase exists",
            "manualPhaseMerges": grids["manual_merges"],
            "c6SymmetryMerges": grids["c6_symmetry_merges"],
            "c6SymmetryRule": "for the exactly C6-covariant TiB2 F(lambda) family, audited non-C6 spine components are attached to the closest C6-valid phase strictly below in E",
            "geometryRule": "click attachments retain the resolved dominant volume component and draw one spine for each outer or nested boundary filling used in its Yamada invariant",
        },
        "modes": modes,
    }
    payload["regions"].update(material_regions)
    return transition


def patch_controller(html: str) -> str:
    html = html.replace(
        "function htmlTransitionTitle(item, mode) {\n"
        '  const startName = item.startName || "0";\n'
        '  const endName = item.endName || "1";\n'
        '  const suffix = mode.key === "classic" ? "" : ` (${mode.short_title || mode.title})`;\n'
        "  return `(1−λ)<i>H</i><sub>${startName}</sub>+λ<i>H</i><sub>${endName}</sub>${suffix}`;\n"
        "}",
        "function htmlTransitionTitle(item, mode) {\n"
        '  const startName = item.startName || "0";\n'
        '  const endName = item.endName || "1";\n'
        '  const suffix = mode.key === "classic" ? "" : ` (${mode.short_title || mode.title})`;\n'
        "  if (item.displayTitle) return `${item.displayTitle}${suffix}`;\n"
        "  return `(1−λ)<i>H</i><sub>${startName}</sub>+λ<i>H</i><sub>${endName}</sub>${suffix}`;\n"
        "}",
    )
    html = html.replace(
        '  const skeletonKind = region.skeleton.edge_count ? "edge skeleton" : "vertex-only skeleton";\n'
        '  const parameterLabel = "E";\n'
        "  const parameterValue = region.parameterValue ?? region.Gamma;\n"
        '  Plotly.react("geometry", geometryTraces(region), layout, { responsive: true, displaylogo: false });\n'
        '  document.getElementById("status").innerHTML = `\n'
        "    <strong>${region.title} - ${region.modeTitle}</strong>\n"
        "    <span>region ${region.regionId}, phase ${region.phaseId}; representative λ=${region.lambda}, ${parameterLabel}=${parameterValue}; ${region.source}, ${region.coreMode}; ${skeletonKind}; core V=${region.nodes}, E=${region.edges}, β=${region.cycleRank}; surface triangles=${region.surface.triangle_count}</span>\n"
        '    <span>${region.startName ? `endpoints: ${region.startName} → ${region.endName}` : ""}</span>\n'
        "    <span>${topologyStatusLine(region.topology)}</span>\n"
        "    <span>Classic regions: ${compactIdList(region.classicRegionIds)}; Classic phases: ${compactIdList(region.classicPhaseIds)}; cells=${region.cellCount}</span>\n"
        '    <span class="poly">Yamada: ${yamadaDisplay(region.polynomial)}</span>\n'
        "  `;\n",
        '  const skeletonKind = region.skeleton.edge_count ? "edge skeleton" : "vertex-only skeleton";\n'
        '  const parameterLabel = "E";\n'
        "  const parameterValue = region.parameterValue ?? region.Gamma;\n"
        '  const materialLine = region.lambdaParameter ? `<span>${region.lambdaParameter}(λ)=${region.materialParameterValue}; band pair=${compactIdList(region.bandPair)}; N=${region.dimension}; ${region.spanLabel}</span>` : "";\n'
        '  const polynomialLabel = region.source === "large-core" ? "Phase label" : "Yamada";\n'
        '  Plotly.react("geometry", geometryTraces(region), layout, { responsive: true, displaylogo: false });\n'
        '  document.getElementById("status").innerHTML = `\n'
        "    <strong>${region.title} - ${region.modeTitle}</strong>\n"
        "    <span>region ${region.regionId}, phase ${region.phaseId}; representative λ=${region.lambda}, ${parameterLabel}=${parameterValue}; ${region.source}, ${region.coreMode}; ${skeletonKind}; core V=${region.nodes}, E=${region.edges}, β=${region.cycleRank}; surface triangles=${region.surface.triangle_count}</span>\n"
        '    <span>${region.startName ? `endpoints: ${region.startName} → ${region.endName}` : ""}</span>\n'
        "    ${materialLine}\n"
        "    <span>${topologyStatusLine(region.topology)}</span>\n"
        "    <span>Classic regions: ${compactIdList(region.classicRegionIds)}; Classic phases: ${compactIdList(region.classicPhaseIds)}; cells=${region.cellCount}</span>\n"
        '    <span class="poly">${polynomialLabel}: ${yamadaDisplay(region.polynomial || region.phaseLabel || region.stableSignature)}</span>\n'
        "  `;\n",
    )
    if "const attachmentLine =" not in html:
        html = html.replace(
            '  const materialLine = region.lambdaParameter ? `<span>${region.lambdaParameter}(λ)=${region.materialParameterValue}; band pair=${compactIdList(region.bandPair)}; N=${region.dimension}; ${region.spanLabel}</span>` : "";\n'
            '  const polynomialLabel = region.source === "large-core" ? "Phase label" : "Yamada";\n',
            '  const materialLine = region.lambdaParameter ? `<span>${region.lambdaParameter}(λ)=${region.materialParameterValue}; band pair=${compactIdList(region.bandPair)}; N=${region.dimension}; ${region.spanLabel}</span>` : "";\n'
            "  const attachment = region.attachmentComponent;\n"
            '  const attachmentLine = attachment && attachment.full_component_count > 1 ? `<span>attached component: kept 1 of ${attachment.full_component_count} volume components; ignored ${attachment.discarded_components} disconnected components (${attachment.discarded_voxels} voxels)</span>` : "";\n'
            '  const polynomialLabel = region.source === "large-core" ? "Phase label" : "Yamada";\n',
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
    global MIN_STABLE_CELLS, ISLAND_MERGE_STRATEGY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-html", type=Path, default=TARGET_HTML)
    parser.add_argument("--data-dir", type=Path, default=MATERIAL_DIR)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("_build/new_phase_maps/material_region_geometry.html"),
    )
    parser.add_argument("--min-stable-cells", type=int, default=8)
    parser.add_argument(
        "--island-merge-strategy", choices=("below", "nearest"), default="below"
    )
    parser.add_argument(
        "--primary-component-min-fraction", type=float, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--historical-display-processing",
        action="store_true",
        required=True,
        help="Acknowledge smoothing, manual/C6 display merges and resolved dominant-component boundary geometry.",
    )
    args = parser.parse_args()
    MIN_STABLE_CELLS = args.min_stable_cells
    ISLAND_MERGE_STRATEGY = args.island_merge_strategy
    if args.primary_component_min_fraction is not None:
        parser.error(
            "--primary-component-min-fraction belongs to the older component-selection rule; omit it for the boundary-resolved dominant component."
        )
    if MIN_STABLE_CELLS < 1:
        parser.error("Use min-stable-cells >= 1.")
    records_json = args.data_dir / "material_parameter_phase_map_records.json"
    summary_json = args.data_dir / "material_parameter_phase_map_summary.json"
    started = time.perf_counter()
    html = args.input_html.read_text(encoding="utf-8")
    payload, prefix, suffix = parse_payload(html)
    records = json.loads(records_json.read_text(encoding="utf-8"))
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    lambdas = np.asarray(summary["lambdas"], dtype=float)

    material_keys = {str(item["key"]) for item in summary["families"]}
    payload["transitions"] = [
        item
        for item in payload["transitions"]
        if str(item.get("key")) not in material_keys
    ]
    payload["regions"] = {
        key: region
        for key, region in payload["regions"].items()
        if str(region.get("transition")) not in material_keys
    }
    for key in material_keys:
        payload.get("stable_min_component_cells", {}).pop(key, None)

    existing_keys = {item["key"] for item in payload["transitions"]}
    payload["version"] = (
        "hamiltonian_yamada_plotly_region_geometry_v19_boundary_resolved"
    )
    payload["material_parameter_scan"] = {
        "records": str(records_json),
        "summary": str(summary_json),
        "note": "Material Hamiltonian coefficient scans updated from boundary-resolved records; transitions absent from this run are retained from the input HTML.",
    }

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
    payload["material_parameter_scan"]["source_counts"] = dict(
        Counter(record["source"] for record in records)
    )
    payload["material_parameter_scan"]["post_processing"] = {
        "min_island_cells": int(MIN_STABLE_CELLS),
        "island_merge_strategy": str(ISLAND_MERGE_STRATEGY),
        "island_rule": "connected lambda-E islands at or below the cutoff are immersed into the spatially closest large-phase cell below in E, with vertical separation used only as a tie-breaker and nearest-large fallback only when no lower phase exists",
        "manual_phase_merges": {
            item["key"]: item["modes"]["classic"].get("manualPhaseMerges", [])
            for item in new_transitions
        },
        "c6_symmetry_merges": {
            item["key"]: item["modes"]["classic"].get("c6SymmetryMerges", [])
            for item in new_transitions
        },
        "geometry_rule": "surface attachments use the resolved dominant volume component; skeletons use the outer and nested-boundary fillings used by the classifier",
    }

    merged = (
        prefix + json.dumps(payload, separators=(",", ":"), ensure_ascii=False) + suffix
    )
    merged = patch_controller(merged)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(merged, encoding="utf-8")
    print(f"\nwrote: {args.output}", flush=True)
    print(f"new transitions: {[item['key'] for item in new_transitions]}", flush=True)
    print(f"total transitions: {len(payload['transitions'])}", flush=True)
    print(f"total regions: {len(payload['regions'])}", flush=True)
    print(f"elapsed: {time.perf_counter() - started:.1f}s", flush=True)


def combine_skeleton_payloads(parts: list[dict[str, Any]]) -> dict[str, Any]:
    combined = {
        "lines": {"x": [], "y": [], "z": []},
        "nodes": {"x": [], "y": [], "z": []},
        "node_count": 0,
        "edge_count": 0,
        "components": 0,
        "cycle_rank": 0,
    }
    for part in parts:
        for axis in ("x", "y", "z"):
            combined["lines"][axis].extend(part["lines"][axis])
            combined["nodes"][axis].extend(part["nodes"][axis])
        for key in ("node_count", "edge_count", "components", "cycle_rank"):
            combined[key] += int(part[key])
    return combined


def boundary_resolved_skeleton_payload(
    obj: MaterialFermiSurface,
    component_info: dict[str, Any],
) -> tuple[dict[str, Any], str]:
    """Render the same outer and nested-boundary fillings used by classification."""
    parts: list[dict[str, Any]] = []
    fillings = []
    for outer_filling, voids in boundary_filling_groups(component_info["mask"]):
        fillings.extend((outer_filling, *voids))

    for filling in fillings:
        topology = topology_for_mask(filling)
        if int(topology["handle_rank"]) == 0:
            center_idx = np.argwhere(filling).mean(axis=0)
            coord = np.asarray(
                obj._idx_to_coord(center_idx.reshape(1, 3))[0], dtype=float
            )
            parts.append(single_node_skeleton(coord))
            continue
        try:
            skeleton_image = np.asarray(skeletonize_volume(filling), dtype=bool)
            obj.skeleton_graph_cache = None
            obj.skeleton_graph_cache_args = None
            graph = obj.skeleton_graph(
                skeleton_image=skeleton_image,
                smooth_epsilon=0,
                simplify=True,
                force_small_edge_contraction=True,
                small_edge_limit=math.pi * 0.1,
                previous_n_edgepoint=20,
            )
            parts.append(skeleton_payload_from_graph(obj, graph))
        except Exception as exc:
            print(
                f"    boundary skeleton extraction failed: {type(exc).__name__}: {exc}",
                flush=True,
            )
            center_idx = np.argwhere(filling).mean(axis=0)
            coord = np.asarray(
                obj._idx_to_coord(center_idx.reshape(1, 3))[0], dtype=float
            )
            parts.append(single_node_skeleton(coord))

    if not parts:
        return single_node_skeleton(
            vertex_coordinate(obj, component_info)
        ), "boundary-empty"
    return combine_skeleton_payloads(parts), f"boundary-resolved-{len(parts)}"


if __name__ == "__main__":
    main()
