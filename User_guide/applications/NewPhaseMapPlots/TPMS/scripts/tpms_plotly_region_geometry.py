#!/usr/bin/env python3
"""Build a clickable TPMS phase-map/geometry HTML from TPMS scan outputs."""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import math
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
from skimage.measure import marching_cubes
from knotted_graph.applications.phase_map_examples._runtime import resolve_geometry_path


THIS_REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCAN_DIR = THIS_REPO / "data"
DEFAULT_TRANSITION_ORDER = (
    "schwarz_p_to_diamond",
    "gyroid_to_schwarz_p",
    "gyroid_to_diamond",
)

BASE_COLORS = (
    "#0f766e",
    "#2563eb",
    "#dc2626",
    "#7c3aed",
    "#ea580c",
    "#059669",
    "#9333ea",
    "#0ea5e9",
    "#be123c",
    "#65a30d",
    "#f59e0b",
    "#475569",
    "#0891b2",
    "#9f1239",
    "#166534",
    "#7f1d1d",
)


def load_tpms_module(scan_script: Path | None):
    if scan_script is None:
        from knotted_graph.applications.phase_map_examples import _tpms

        return _tpms
    spec = importlib.util.spec_from_file_location(scan_script.stem, scan_script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {scan_script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_path(path_text: str, scan_dir: Path = DEFAULT_SCAN_DIR) -> Path:
    # Legacy records contain the researcher's old tmp/ prefix. Relocate using
    # the explicitly selected dataset, not a machine-specific fallback repo.
    return resolve_geometry_path(path_text, scan_dir)


def load_geometry(record: dict[str, Any]) -> dict[str, Any]:
    path = Path(record["geometry_path"])
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def graph_from_geometry(geometry: dict[str, Any]) -> nx.MultiGraph:
    graph = nx.MultiGraph()
    for node in geometry["graph"].get("nodes", []):
        graph.add_node(str(node["id"]))
    for edge in geometry["graph"].get("edges", []):
        graph.add_edge(str(edge["u"]), str(edge["v"]))
    return graph


def graph_summary(graph: nx.MultiGraph) -> dict[str, Any]:
    nodes = graph.number_of_nodes()
    edges = graph.number_of_edges()
    components = nx.number_connected_components(graph) if nodes else 0
    cycle_rank = edges - nodes + components
    return {
        "nodes": int(nodes),
        "edges": int(edges),
        "components": int(components),
        "cycle_rank": int(cycle_rank),
        "degree_sequence": sorted((int(d) for _, d in graph.degree()), reverse=True),
    }


def multigraph_hash(graph: nx.MultiGraph) -> tuple[Any, ...]:
    simple = nx.Graph()
    loop_counts: Counter[Any] = Counter()
    multiplicities: Counter[tuple[Any, Any]] = Counter()
    for node in graph.nodes:
        simple.add_node(node)
    for u, v in graph.edges():
        if u == v:
            loop_counts[u] += 1
        else:
            key = tuple(sorted((u, v), key=str))
            multiplicities[key] += 1
    for node in simple.nodes:
        simple.nodes[node]["loops"] = str(loop_counts[node])
    for (u, v), count in multiplicities.items():
        simple.add_edge(u, v, m=str(count))
    wl_hash = nx.weisfeiler_lehman_graph_hash(
        simple,
        node_attr="loops",
        edge_attr="m",
        iterations=3,
    )
    summary = graph_summary(graph)
    return (
        summary["nodes"],
        summary["edges"],
        summary["components"],
        summary["cycle_rank"],
        tuple(summary["degree_sequence"]),
        wl_hash,
    )


def is_seen_graph(
    graph: nx.MultiGraph,
    buckets: dict[tuple[Any, ...], list[nx.MultiGraph]],
) -> bool:
    key = multigraph_hash(graph)
    bucket = buckets[key]
    for candidate in bucket:
        if nx.is_isomorphic(graph, candidate):
            return True
    bucket.append(graph.copy())
    return False


def contract_edge(
    graph: nx.MultiGraph,
    selected: tuple[Any, Any, Any],
) -> nx.MultiGraph:
    u, v, selected_key = selected
    if u == v:
        raise ValueError("cannot contract a loop edge")
    keep = u
    kill = v
    contracted = nx.MultiGraph()
    for node in graph.nodes:
        contracted.add_node(keep if node == kill else node)
    for a, b, key in graph.edges(keys=True):
        if a == u and b == v and key == selected_key:
            continue
        na = keep if a == kill else a
        nb = keep if b == kill else b
        contracted.add_edge(na, nb)
    return contracted


def can_contract_to(
    source: nx.MultiGraph,
    target: nx.MultiGraph,
    *,
    max_states: int,
    max_depth: int,
    max_edges: int,
) -> bool:
    source_summary = graph_summary(source)
    target_summary = graph_summary(target)
    if source_summary["components"] != target_summary["components"]:
        return False
    if source_summary["cycle_rank"] != target_summary["cycle_rank"]:
        return False
    if source_summary["nodes"] < target_summary["nodes"]:
        return False
    if source_summary["edges"] < target_summary["edges"]:
        return False
    if max(source_summary["edges"], target_summary["edges"]) > max_edges:
        return False
    node_drop = source_summary["nodes"] - target_summary["nodes"]
    edge_drop = source_summary["edges"] - target_summary["edges"]
    if node_drop != edge_drop:
        return False
    if node_drop > max_depth:
        return False
    if (
        source_summary["nodes"] == target_summary["nodes"]
        and source_summary["edges"] == target_summary["edges"]
    ):
        return nx.is_isomorphic(source, target)

    queue: deque[nx.MultiGraph] = deque([source.copy()])
    seen: dict[tuple[Any, ...], list[nx.MultiGraph]] = defaultdict(list)
    is_seen_graph(source, seen)
    states = 0
    while queue:
        current = queue.popleft()
        states += 1
        if states > max_states:
            return False
        current_summary = graph_summary(current)
        if current_summary["nodes"] < target_summary["nodes"]:
            continue
        if current_summary["edges"] < target_summary["edges"]:
            continue
        if (
            current_summary["nodes"] == target_summary["nodes"]
            and current_summary["edges"] == target_summary["edges"]
            and nx.is_isomorphic(current, target)
        ):
            return True
        for edge in list(current.edges(keys=True)):
            if edge[0] == edge[1]:
                continue
            contracted = contract_edge(current, edge)
            contracted_summary = graph_summary(contracted)
            if contracted_summary["nodes"] < target_summary["nodes"]:
                continue
            if contracted_summary["edges"] < target_summary["edges"]:
                continue
            if contracted_summary["cycle_rank"] != target_summary["cycle_rank"]:
                continue
            if is_seen_graph(contracted, seen):
                continue
            queue.append(contracted)
    return False


class DisjointSet:
    def __init__(self, values: list[int]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: int) -> int:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: int, right: int) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if root_right < root_left:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left


def connected_regions(label_grid: np.ndarray) -> list[dict[str, Any]]:
    rows, cols = label_grid.shape
    seen = np.zeros(label_grid.shape, dtype=bool)
    regions: list[dict[str, Any]] = []
    for row in range(rows):
        for col in range(cols):
            if seen[row, col]:
                continue
            label = int(label_grid[row, col])
            queue = deque([(row, col)])
            seen[row, col] = True
            cells: list[tuple[int, int]] = []
            while queue:
                current_row, current_col = queue.popleft()
                cells.append((current_row, current_col))
                for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    next_row = current_row + dr
                    next_col = current_col + dc
                    if not (0 <= next_row < rows and 0 <= next_col < cols):
                        continue
                    if seen[next_row, next_col]:
                        continue
                    if int(label_grid[next_row, next_col]) != label:
                        continue
                    seen[next_row, next_col] = True
                    queue.append((next_row, next_col))
            regions.append({"phase_id": label, "cells": cells})
    return regions


def phase_representatives(
    grid: np.ndarray,
    record_grid: list[list[dict[str, Any]]],
) -> dict[int, dict[str, Any]]:
    reps: dict[int, dict[str, Any]] = {}
    for phase_id in sorted(set(int(value) for value in grid.ravel())):
        cells = [
            (row, col)
            for row in range(grid.shape[0])
            for col in range(grid.shape[1])
            if int(grid[row, col]) == phase_id
        ]
        rep = choose_representative(cells, record_grid, phase_id=phase_id)
        reps[phase_id] = rep
    return reps


def choose_representative(
    cells: list[tuple[int, int]],
    record_grid: list[list[dict[str, Any]]],
    *,
    phase_id: int | None = None,
) -> dict[str, Any]:
    center_row = sum(row for row, _ in cells) / len(cells)
    center_col = sum(col for _, col in cells) / len(cells)
    candidate_cells = cells
    if phase_id is not None:
        filtered = [
            (row, col)
            for row, col in cells
            if int(record_grid[row][col].get("_classic_phase_id", -1)) == int(phase_id)
        ]
        if filtered:
            candidate_cells = filtered

    def score(cell: tuple[int, int]) -> tuple[float, float]:
        row, col = cell
        record = record_grid[row][col]
        source_rank = {
            "yamada": 0,
            "large-core": 1,
            "vertex": 2,
            "error": 3,
        }.get(record.get("source", ""), 4)
        distance = (row - center_row) ** 2 + (col - center_col) ** 2
        return float(source_rank), float(distance)

    row, col = min(candidate_cells, key=score)
    return record_grid[row][col]


def build_record_grid(
    records: list[dict[str, Any]],
    family_key: str,
    lambdas: list[float],
    thresholds: list[float],
) -> list[list[dict[str, Any]]]:
    lookup = {
        (
            round(float(record["lam"]), 12),
            round(float(record["threshold_c"]), 12),
        ): record
        for record in records
        if record["family"] == family_key
    }
    grid: list[list[dict[str, Any]]] = []
    for threshold in thresholds:
        row = []
        for lam in lambdas:
            row.append(lookup[(round(float(lam), 12), round(float(threshold), 12))])
        grid.append(row)
    return grid


def contraction_phase_map(
    phase_ids: list[int],
    phase_graphs: dict[int, nx.MultiGraph],
    *,
    max_states: int,
    max_depth: int,
    max_edges: int,
) -> dict[int, int]:
    dsu = DisjointSet(phase_ids)
    for index, left in enumerate(phase_ids):
        for right in phase_ids[index + 1 :]:
            g_left = phase_graphs[left]
            g_right = phase_graphs[right]
            if can_contract_to(
                g_left,
                g_right,
                max_states=max_states,
                max_depth=max_depth,
                max_edges=max_edges,
            ) or can_contract_to(
                g_right,
                g_left,
                max_states=max_states,
                max_depth=max_depth,
                max_edges=max_edges,
            ):
                dsu.union(left, right)
    roots = sorted({dsu.find(phase_id) for phase_id in phase_ids})
    root_to_class = {root: index + 1 for index, root in enumerate(roots)}
    return {phase_id: root_to_class[dsu.find(phase_id)] for phase_id in phase_ids}


def remapped_surface_mesh(
    values: np.ndarray,
    level: float,
    *,
    spacing: np.ndarray,
    origin: np.ndarray,
    max_faces: int,
    max_step: int = 6,
) -> dict[str, Any]:
    if not (float(np.nanmin(values)) <= level <= float(np.nanmax(values))):
        return {
            "x": [],
            "y": [],
            "z": [],
            "i": [],
            "j": [],
            "k": [],
            "triangle_count": 0,
            "marching_cubes_step_size": 1,
        }
    verts = np.empty((0, 3), dtype=float)
    faces = np.empty((0, 3), dtype=np.int64)
    step_used = 1
    for step in range(1, int(max_step) + 1):
        verts, faces, _, _ = marching_cubes(
            values,
            level=float(level),
            spacing=tuple(float(v) for v in spacing),
            step_size=step,
        )
        step_used = step
        if len(faces) <= max_faces:
            break
    verts = verts + origin[None, :]
    if not len(faces):
        return {
            "x": [],
            "y": [],
            "z": [],
            "i": [],
            "j": [],
            "k": [],
            "triangle_count": 0,
            "marching_cubes_step_size": int(step_used),
        }
    used = np.unique(faces.ravel())
    remap = {int(old): new for new, old in enumerate(used)}
    compact_faces = np.vectorize(lambda value: remap[int(value)])(faces)
    compact_verts = verts[used]
    compact_verts = np.round(compact_verts, 4)
    return {
        "x": compact_verts[:, 0].tolist(),
        "y": compact_verts[:, 1].tolist(),
        "z": compact_verts[:, 2].tolist(),
        "i": compact_faces[:, 0].astype(int).tolist(),
        "j": compact_faces[:, 1].astype(int).tolist(),
        "k": compact_faces[:, 2].astype(int).tolist(),
        "triangle_count": int(len(compact_faces)),
        "marching_cubes_step_size": int(step_used),
    }


def skeleton_payload(geometry: dict[str, Any]) -> dict[str, Any]:
    node_x: list[float] = []
    node_y: list[float] = []
    node_z: list[float] = []
    for node in geometry["graph"].get("nodes", []):
        coord = node.get("coord", [0.0, 0.0, 0.0])
        node_x.append(round(float(coord[0]), 4))
        node_y.append(round(float(coord[1]), 4))
        node_z.append(round(float(coord[2]), 4))

    line_x: list[float | None] = []
    line_y: list[float | None] = []
    line_z: list[float | None] = []
    for edge in geometry["graph"].get("edges", []):
        points = edge.get("points_coord", [])
        for point in points:
            line_x.append(round(float(point[0]), 4))
            line_y.append(round(float(point[1]), 4))
            line_z.append(round(float(point[2]), 4))
        line_x.append(None)
        line_y.append(None)
        line_z.append(None)
    return {
        "nodes": {"x": node_x, "y": node_y, "z": node_z},
        "lines": {"x": line_x, "y": line_y, "z": line_z},
        "edge_count": int(len(geometry["graph"].get("edges", []))),
    }


def compact_record_value(value: Any) -> Any:
    if isinstance(value, float):
        if math.isfinite(value):
            return round(value, 6)
        return None
    if isinstance(value, list):
        return [compact_record_value(item) for item in value]
    return value


def build_region_payload(
    *,
    mode_key: str,
    mode_title: str,
    family: Any,
    region_id: int,
    phase_id: int,
    cells: list[tuple[int, int]],
    record_grid: list[list[dict[str, Any]]],
    classic_region_grid: np.ndarray,
    tpms_module: Any,
    max_surface_faces: int,
    color: str,
    representative_phase_id: int | None,
) -> dict[str, Any]:
    record = choose_representative(cells, record_grid, phase_id=representative_phase_id)
    geometry = load_geometry(record)
    _, x, y, z = tpms_module.sample_grid(family.span, family.dimension)
    raw_values = family.field_at(float(record["lam"]), x, y, z)
    if hasattr(family, "domain") and hasattr(family, "compact_level_field"):
        domain_values = family.domain.sample(x, y, z)
        values = family.compact_level_field(
            raw_values,
            domain_values,
            float(record["threshold_c"]),
        )
        surface_level = 0.0
    else:
        values = raw_values
        surface_level = float(record["threshold_c"])
    region = tpms_module.ImplicitSolidRegion(
        values,
        threshold_c=float(surface_level),
        span=family.span,
    )
    surface = remapped_surface_mesh(
        values,
        float(surface_level),
        spacing=region.spacing * region.axis_scale,
        origin=region.origin,
        max_faces=max_surface_faces,
    )
    classic_region_ids = sorted(
        {int(classic_region_grid[row, col]) for row, col in cells}
    )
    classic_phase_ids = sorted(
        {int(record_grid[row][col]["_classic_phase_id"]) for row, col in cells}
    )
    topology = {
        "components": int(record["interior_components"]),
        "euler_characteristic": int(record["euler_characteristic"]),
        "handle_rank": int(record["handle_rank"]),
        "touches_boundary": bool(record["touches_boundary"]),
        "boundary_faces": list(record["boundary_faces"]),
        "volume_fraction": compact_record_value(float(record["interior_fraction"])),
        "interior_voxels": int(record["interior_voxels"]),
        "skeleton_voxels": int(record["skeleton_voxels"]),
        "surface_points": int(record["surface_points"]),
        "surface_cells": int(record["surface_cells"]),
        "surface_open_edges": int(record.get("surface_open_edges", 0)),
        "surface_nonmanifold_edges": int(record.get("surface_nonmanifold_edges", 0)),
        "surface_is_closed": bool(record.get("surface_is_closed", False)),
    }
    domain = getattr(family, "domain", None)
    return {
        "title": family.title,
        "mode": mode_key,
        "modeTitle": mode_title,
        "regionId": int(region_id),
        "phaseId": int(phase_id),
        "classicRegionIds": classic_region_ids,
        "classicPhaseIds": classic_phase_ids,
        "lambda": compact_record_value(float(record["lam"])),
        "parameterValue": compact_record_value(float(record["threshold_c"])),
        "parameterName": "c",
        "startName": family.start.title,
        "endName": family.end.title,
        "source": record["source"],
        "coreMode": "same-volume-skeleton-pipeline",
        "compactDomain": {
            "key": getattr(domain, "key", "none"),
            "title": getattr(domain, "title", "none"),
            "formula": getattr(domain, "formula", ""),
        },
        "nodes": int(record["nodes"]),
        "edges": int(record["edges"]),
        "components": int(record["components"]),
        "cycleRank": int(record["cycle_rank"]),
        "degreeSequence": list(record["degree_sequence"]),
        "totalEdgePoints": int(record["total_edge_points"]),
        "isTrivalent": bool(record["is_trivalent"]),
        "polynomial": record["polynomial"] or record["phase_label"],
        "phaseLabel": record["phase_label"],
        "phaseSignature": record["phase_signature"],
        "cellCount": int(len(cells)),
        "color": color,
        "surface": surface,
        "skeleton": skeleton_payload(geometry),
        "topology": topology,
    }


def colors_for_count(count: int) -> list[str]:
    return [BASE_COLORS[index % len(BASE_COLORS)] for index in range(max(1, count))]


def family_from_summary(
    tpms_module: Any, family_summary: dict[str, Any], scan_parameters: dict[str, Any]
):
    """Use recorded domain parameters, never silently regenerate a different solid."""
    allowed = {"domain_kind", "span_half_width", "domain_radius_fraction"}
    options = {key: value for key, value in scan_parameters.items() if key in allowed}
    family = next(
        f
        for f in tpms_module.tpms_families(
            dimension=int(family_summary["dimension"]),
            thresholds=tuple(float(value) for value in family_summary["thresholds"]),
            **options,
        )
        if f.key == family_summary["key"]
    )
    domain = family_summary["compact_domain"]
    if (
        not np.array_equal(family.span, family_summary["span"])
        or family.domain.key != domain["key"]
        or family.domain.formula != domain["formula"]
    ):
        raise ValueError(
            "Recorded TPMS domain differs from the reconstructed family. "
            "Provide exact scan_parameters in the summary; do not substitute default geometry."
        )
    return family


def build_payload(
    scan_dir: Path,
    *,
    scan_script: Path | None,
    max_surface_faces: int,
    max_contraction_states: int,
    max_contraction_depth: int,
    max_contraction_edges: int,
    stable_min_cells: int,
) -> dict[str, Any]:
    tpms_module = load_tpms_module(scan_script)
    records = load_json(scan_dir / "tpms_parameter_phase_map_records.json")
    for record in records:
        record["geometry_path"] = str(resolve_path(record["geometry_path"], scan_dir))
    source_data = load_json(scan_dir / "tpms_parameter_phase_map_source_data.json")
    summary = load_json(scan_dir / "tpms_parameter_phase_map_summary.json")
    lambdas = [float(value) for value in source_data["lambdas"]]

    transitions = []
    regions_payload: dict[str, Any] = {}
    global_region_rows = []

    order_index = {key: index for index, key in enumerate(DEFAULT_TRANSITION_ORDER)}
    family_summaries = sorted(
        summary["families"],
        key=lambda item: (order_index.get(item["key"], len(order_index)), item["key"]),
    )

    for family_summary in family_summaries:
        family_key = family_summary["key"]
        thresholds = [float(value) for value in family_summary["thresholds"]]
        family = family_from_summary(
            tpms_module, family_summary, summary.get("scan_parameters", {})
        )
        map_data = source_data["phase_maps"][family_key]
        classic_grid = np.asarray(map_data["stable_grid"], dtype=int)
        record_grid = build_record_grid(records, family_key, lambdas, thresholds)
        for row in range(classic_grid.shape[0]):
            for col in range(classic_grid.shape[1]):
                record_grid[row][col]["_classic_phase_id"] = int(classic_grid[row, col])

        phase_reps = phase_representatives(classic_grid, record_grid)
        phase_graphs = {
            phase_id: graph_from_geometry(load_geometry(record))
            for phase_id, record in phase_reps.items()
        }
        phase_ids = sorted(phase_graphs)
        contraction_class = contraction_phase_map(
            phase_ids,
            phase_graphs,
            max_states=max_contraction_states,
            max_depth=max_contraction_depth,
            max_edges=max_contraction_edges,
        )
        contraction_grid = np.vectorize(lambda value: contraction_class[int(value)])(
            classic_grid
        ).astype(int)

        transition_payload = {
            "key": family.key,
            "title": family.title,
            "lambdas": [round(value, 6) for value in lambdas],
            "thresholds": [round(value, 6) for value in thresholds],
            "parameters": [round(value, 6) for value in thresholds],
            "parameterName": "c",
            "parameterLabel": "c",
            "startName": family.start.title,
            "endName": family.end.title,
            "formula": f"(1-lambda){family.start.key} + lambda {family.end.key}",
            "modes": {},
        }

        classic_region_grid = np.zeros_like(classic_grid, dtype=int)
        classic_region_keys = [["" for _ in lambdas] for _ in thresholds]
        classic_colors = colors_for_count(int(classic_grid.max()))
        classic_regions = connected_regions(classic_grid)
        for region_index, region_info in enumerate(classic_regions, start=1):
            for row, col in region_info["cells"]:
                classic_region_grid[row, col] = region_index
                classic_region_keys[row][col] = f"classic:{family.key}:{region_index}"
            color = classic_colors[int(region_info["phase_id"]) - 1]
            region_payload = build_region_payload(
                mode_key="classic",
                mode_title="Classic",
                family=family,
                region_id=region_index,
                phase_id=int(region_info["phase_id"]),
                cells=region_info["cells"],
                record_grid=record_grid,
                classic_region_grid=classic_region_grid,
                tpms_module=tpms_module,
                max_surface_faces=max_surface_faces,
                color=color,
                representative_phase_id=int(region_info["phase_id"]),
            )
            regions_payload[f"classic:{family.key}:{region_index}"] = region_payload
            global_region_rows.append(
                {
                    "transition": family.key,
                    "mode": "classic",
                    "region_id": region_index,
                    "phase_id": int(region_info["phase_id"]),
                    "cell_count": len(region_info["cells"]),
                    "polynomial": region_payload["polynomial"],
                    "nodes": region_payload["nodes"],
                    "edges": region_payload["edges"],
                    "cycle_rank": region_payload["cycleRank"],
                }
            )

        contraction_region_grid = np.zeros_like(contraction_grid, dtype=int)
        contraction_region_keys = [["" for _ in lambdas] for _ in thresholds]
        contraction_colors = colors_for_count(int(contraction_grid.max()))
        contraction_regions = connected_regions(contraction_grid)
        for region_index, region_info in enumerate(contraction_regions, start=1):
            for row, col in region_info["cells"]:
                contraction_region_grid[row, col] = region_index
                contraction_region_keys[row][col] = (
                    f"contraction:{family.key}:{region_index}"
                )
            color = contraction_colors[int(region_info["phase_id"]) - 1]
            region_payload = build_region_payload(
                mode_key="contraction",
                mode_title="Up to contraction moves",
                family=family,
                region_id=region_index,
                phase_id=int(region_info["phase_id"]),
                cells=region_info["cells"],
                record_grid=record_grid,
                classic_region_grid=classic_region_grid,
                tpms_module=tpms_module,
                max_surface_faces=max_surface_faces,
                color=color,
                representative_phase_id=None,
            )
            regions_payload[f"contraction:{family.key}:{region_index}"] = region_payload
            global_region_rows.append(
                {
                    "transition": family.key,
                    "mode": "contraction",
                    "region_id": region_index,
                    "phase_id": int(region_info["phase_id"]),
                    "classic_region_ids": region_payload["classicRegionIds"],
                    "classic_phase_ids": region_payload["classicPhaseIds"],
                    "cell_count": len(region_info["cells"]),
                    "polynomial": region_payload["polynomial"],
                    "nodes": region_payload["nodes"],
                    "edges": region_payload["edges"],
                    "cycle_rank": region_payload["cycleRank"],
                }
            )

        transition_payload["modes"]["classic"] = {
            "z": classic_grid.astype(int).tolist(),
            "regionKeys": classic_region_keys,
            "colors": classic_colors,
            "classes": int(len(set(int(value) for value in classic_grid.ravel()))),
            "components": int(len(classic_regions)),
            "rawClasses": int(len(set(int(value) for value in classic_grid.ravel()))),
        }
        if stable_min_cells > 1:
            stable_grid, stable_changed_cells = tpms_module.stable_labels(
                classic_grid,
                min_cells=int(stable_min_cells),
            )
            stable_region_grid = np.zeros_like(stable_grid, dtype=int)
            stable_region_keys = [["" for _ in lambdas] for _ in thresholds]
            stable_colors = colors_for_count(
                int(max(classic_grid.max(), stable_grid.max()))
            )
            stable_regions = connected_regions(stable_grid)
            for region_index, region_info in enumerate(stable_regions, start=1):
                for row, col in region_info["cells"]:
                    stable_region_grid[row, col] = region_index
                    stable_region_keys[row][col] = f"stable:{family.key}:{region_index}"
                color = stable_colors[int(region_info["phase_id"]) - 1]
                region_payload = build_region_payload(
                    mode_key="stable",
                    mode_title=f"Stable display (min {stable_min_cells} cells)",
                    family=family,
                    region_id=region_index,
                    phase_id=int(region_info["phase_id"]),
                    cells=region_info["cells"],
                    record_grid=record_grid,
                    classic_region_grid=classic_region_grid,
                    tpms_module=tpms_module,
                    max_surface_faces=max_surface_faces,
                    color=color,
                    representative_phase_id=int(region_info["phase_id"]),
                )
                regions_payload[f"stable:{family.key}:{region_index}"] = region_payload
                global_region_rows.append(
                    {
                        "transition": family.key,
                        "mode": "stable",
                        "region_id": region_index,
                        "phase_id": int(region_info["phase_id"]),
                        "classic_region_ids": region_payload["classicRegionIds"],
                        "classic_phase_ids": region_payload["classicPhaseIds"],
                        "cell_count": len(region_info["cells"]),
                        "polynomial": region_payload["polynomial"],
                        "nodes": region_payload["nodes"],
                        "edges": region_payload["edges"],
                        "cycle_rank": region_payload["cycleRank"],
                    }
                )
            transition_payload["modes"]["stable"] = {
                "z": stable_grid.astype(int).tolist(),
                "regionKeys": stable_region_keys,
                "colors": stable_colors,
                "classes": int(len(set(int(value) for value in stable_grid.ravel()))),
                "components": int(len(stable_regions)),
                "rawClasses": int(
                    len(set(int(value) for value in classic_grid.ravel()))
                ),
                "stableMinCells": int(stable_min_cells),
                "stableReassignedCells": int(stable_changed_cells),
            }
        transition_payload["modes"]["contraction"] = {
            "z": contraction_grid.astype(int).tolist(),
            "regionKeys": contraction_region_keys,
            "colors": contraction_colors,
            "classes": int(len(set(int(value) for value in contraction_grid.ravel()))),
            "components": int(len(contraction_regions)),
            "contractionMergedClassicRegions": int(
                max(0, len(classic_regions) - len(contraction_regions))
            ),
            "phaseToContractionClass": {
                str(phase_id): int(contraction_class[phase_id])
                for phase_id in phase_ids
            },
        }
        if stable_min_cells > 1:
            stable_contraction_grid, stable_contraction_changed_cells = (
                tpms_module.stable_labels(
                    contraction_grid,
                    min_cells=int(stable_min_cells),
                )
            )
            stable_contraction_region_grid = np.zeros_like(
                stable_contraction_grid, dtype=int
            )
            stable_contraction_region_keys = [["" for _ in lambdas] for _ in thresholds]
            stable_contraction_colors = colors_for_count(
                int(max(contraction_grid.max(), stable_contraction_grid.max()))
            )
            stable_contraction_regions = connected_regions(stable_contraction_grid)
            for region_index, region_info in enumerate(
                stable_contraction_regions, start=1
            ):
                for row, col in region_info["cells"]:
                    stable_contraction_region_grid[row, col] = region_index
                    stable_contraction_region_keys[row][col] = (
                        f"stable_contraction:{family.key}:{region_index}"
                    )
                color = stable_contraction_colors[int(region_info["phase_id"]) - 1]
                region_payload = build_region_payload(
                    mode_key="stable_contraction",
                    mode_title=f"Stable up to contraction (min {stable_min_cells} cells)",
                    family=family,
                    region_id=region_index,
                    phase_id=int(region_info["phase_id"]),
                    cells=region_info["cells"],
                    record_grid=record_grid,
                    classic_region_grid=classic_region_grid,
                    tpms_module=tpms_module,
                    max_surface_faces=max_surface_faces,
                    color=color,
                    representative_phase_id=None,
                )
                regions_payload[f"stable_contraction:{family.key}:{region_index}"] = (
                    region_payload
                )
                global_region_rows.append(
                    {
                        "transition": family.key,
                        "mode": "stable_contraction",
                        "region_id": region_index,
                        "phase_id": int(region_info["phase_id"]),
                        "classic_region_ids": region_payload["classicRegionIds"],
                        "classic_phase_ids": region_payload["classicPhaseIds"],
                        "cell_count": len(region_info["cells"]),
                        "polynomial": region_payload["polynomial"],
                        "nodes": region_payload["nodes"],
                        "edges": region_payload["edges"],
                        "cycle_rank": region_payload["cycleRank"],
                    }
                )
            transition_payload["modes"]["stable_contraction"] = {
                "z": stable_contraction_grid.astype(int).tolist(),
                "regionKeys": stable_contraction_region_keys,
                "colors": stable_contraction_colors,
                "classes": int(
                    len(set(int(value) for value in stable_contraction_grid.ravel()))
                ),
                "components": int(len(stable_contraction_regions)),
                "rawClasses": int(
                    len(set(int(value) for value in contraction_grid.ravel()))
                ),
                "stableMinCells": int(stable_min_cells),
                "stableReassignedCells": int(stable_contraction_changed_cells),
                "sourceMode": "contraction",
            }
        transitions.append(transition_payload)

    modes = []
    if stable_min_cells > 1:
        modes.append(
            {
                "key": "stable_contraction",
                "title": f"Stable up to contraction (min {stable_min_cells})",
                "short_title": "Stable contraction",
                "description": (
                    "Contraction phase map with connected components smaller "
                    f"than {stable_min_cells} cells reassigned to neighboring regions. "
                    "Raw cell signatures and representative audit data remain in the CSV/JSON."
                ),
            }
        )
        modes.append(
            {
                "key": "stable",
                "title": f"Stable raw audit (min {stable_min_cells})",
                "short_title": "Stable raw",
                "description": (
                    "Display-only raw Yamada/audit phase map with connected components "
                    f"smaller than {stable_min_cells} cells reassigned to neighboring regions."
                ),
            }
        )
    modes.extend(
        [
            {
                "key": "classic",
                "title": "Raw audit",
                "short_title": "Raw",
                "description": "Connected-region partition of Yamada/audit signatures.",
            },
            {
                "key": "contraction",
                "title": "Up to contraction moves",
                "short_title": "Contraction",
                "description": (
                    "Classic regions are merged when representative core graphs "
                    "are related by non-loop edge contractions within the state/depth/edge limits."
                ),
            },
        ]
    )
    return {
        "version": "tpms_yamada_plotly_region_geometry_v1",
        "scanDir": str(scan_dir),
        "scanScript": str(scan_script)
        if scan_script
        else "knotted_graph.applications.phase_map_examples._tpms",
        "defaultMode": "stable_contraction" if stable_min_cells > 1 else "contraction",
        "lambdas": [round(value, 6) for value in lambdas],
        "modes": modes,
        "transitions": transitions,
        "regions": regions_payload,
        "regionRows": global_region_rows,
        "summary": {
            "reconstructedVolumes": int(
                summary["map_info"]["timing"]["reconstructed_volumes"]
            ),
            "elapsedSeconds": float(summary["map_info"]["timing"]["elapsed_seconds"]),
            "sourceCounts": summary["source_counts"],
            "errors": len(summary["errors"]),
            "maxSurfaceFacesPerRepresentative": int(max_surface_faces),
            "maxContractionStatesPerTest": int(max_contraction_states),
            "maxContractionDepthPerTest": int(max_contraction_depth),
            "maxContractionEdgesPerTest": int(max_contraction_edges),
            "stableDisplayMinCells": int(stable_min_cells),
            "scopeNote": summary["scientific_scope_note"],
        },
    }


def html_document(payload: dict[str, Any]) -> str:
    payload_json = json.dumps(
        payload, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>TPMS Yamada geometry QA</title>
  <script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
  <style>
    :root {{ color-scheme: light; --ink:#111111; --muted:#555555; --line:#111111; --paper:#ffffff; --soft:#f4f4f4; }}
    body {{ margin:0; background:var(--paper); color:var(--ink); font-family:"Times New Roman", Times, "DejaVu Serif", serif; }}
    #app {{ max-width:1320px; margin:0 auto; padding:14px 16px 18px; }}
    .topbar {{ display:flex; flex-wrap:wrap; align-items:center; gap:10px 18px; margin-bottom:10px; }}
    .buttonGroup {{ display:flex; flex-wrap:wrap; align-items:center; gap:8px; }}
    .buttonGroup::before {{ content:attr(data-label); font-size:12px; font-weight:700; color:var(--muted); text-transform:uppercase; letter-spacing:0.04em; margin-right:2px; font-family:Arial, Helvetica, sans-serif; }}
    .topbar button {{ border:1px solid var(--line); background:var(--paper); color:var(--ink); padding:6px 10px; font-size:13px; cursor:pointer; font-family:Arial, Helvetica, sans-serif; }}
    .topbar button[aria-pressed="true"] {{ background:var(--ink); color:var(--paper); }}
    .layout {{ display:grid; grid-template-columns:minmax(380px, 0.92fr) minmax(440px, 1.08fr); gap:14px; align-items:start; }}
    #phaseMap, #geometry {{ width:100%; height:640px; border:1.5px solid var(--line); box-sizing:border-box; }}
    .status {{ min-height:112px; margin-top:10px; display:grid; grid-template-columns:1fr; gap:3px; font-size:13px; color:var(--muted); font-family:Arial, Helvetica, sans-serif; }}
    .status strong {{ color:var(--ink); font-size:15px; }}
    .poly {{ font-family:Menlo, Consolas, monospace; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; max-width:100%; }}
    @media (max-width: 900px) {{ .layout {{ grid-template-columns:1fr; }} #phaseMap, #geometry {{ height:520px; }} }}
  </style>
</head>
<body>
<div id="app">
  <div class="topbar">
    <div class="buttonGroup" id="modeButtons" data-label="Mode" aria-label="Phase modes"></div>
    <div class="buttonGroup" id="transitionButtons" data-label="Transition" aria-label="Transitions"></div>
  </div>
  <div class="layout">
    <div>
      <div id="phaseMap" aria-label="Clickable TPMS phase map"></div>
    </div>
    <div>
      <div id="geometry" aria-label="Representative TPMS surface and skeleton"></div>
      <div class="status" id="status" aria-live="polite"></div>
    </div>
  </div>
</div>
<script>
const payload = {payload_json};
const transitions = payload.transitions;
const regions = payload.regions;
const modes = payload.modes;
let activeTransition = transitions[0].key;
let activeMode = payload.defaultMode || modes[0].key;

function discreteColorscale(colors) {{
  const scale = [];
  const n = colors.length;
  if (n === 1) return [[0, colors[0]], [1, colors[0]]];
  colors.forEach((color, index) => {{
    scale.push([index / n, color]);
    scale.push([(index + 1) / n, color]);
  }});
  return scale;
}}

function transitionByKey(key) {{
  return transitions.find(item => item.key === key) || transitions[0];
}}

function modeByKey(key) {{
  return modes.find(item => item.key === key) || modes[0];
}}

function htmlTransitionTitle(item, mode) {{
  const suffix = mode.key === "classic" ? "" : ` (${{mode.short_title || mode.title}})`;
  return `(1−λ)<i>F</i><sub>${{item.startName}}</sub>+λ<i>F</i><sub>${{item.endName}}</sub>${{suffix}}`;
}}

function edgeCoordinates(values) {{
  if (values.length === 1) return [values[0] - 0.5, values[0] + 0.5];
  const edges = [];
  edges.push(values[0] - 0.5 * (values[1] - values[0]));
  for (let index = 1; index < values.length; index += 1) {{
    edges.push(0.5 * (values[index - 1] + values[index]));
  }}
  edges.push(values[values.length - 1] + 0.5 * (values[values.length - 1] - values[values.length - 2]));
  return edges;
}}

function phaseBoundaryTrace(xValues, yValues, labels) {{
  const xEdges = edgeCoordinates(xValues);
  const yEdges = edgeCoordinates(yValues);
  const xs = [];
  const ys = [];
  for (let row = 0; row < labels.length; row += 1) {{
    for (let col = 0; col < labels[row].length - 1; col += 1) {{
      if (labels[row][col] === labels[row][col + 1]) continue;
      xs.push(xEdges[col + 1], xEdges[col + 1], null);
      ys.push(yEdges[row], yEdges[row + 1], null);
    }}
  }}
  for (let row = 0; row < labels.length - 1; row += 1) {{
    for (let col = 0; col < labels[row].length; col += 1) {{
      if (labels[row][col] === labels[row + 1][col]) continue;
      xs.push(xEdges[col], xEdges[col + 1], null);
      ys.push(yEdges[row + 1], yEdges[row + 1], null);
    }}
  }}
  return {{
    type: "scatter",
    mode: "lines",
    x: xs,
    y: ys,
    line: {{ color: "#111111", width: 1.2 }},
    hoverinfo: "skip",
    showlegend: false,
    name: "phase boundary"
  }};
}}

function plotPhaseMap(transitionKey = activeTransition, modeKey = activeMode) {{
  const item = transitionByKey(transitionKey);
  const mode = modeByKey(modeKey);
  const data = item.modes[mode.key];
  const isContractionMode = mode.key === "contraction" || mode.key === "stable_contraction";
  activeTransition = item.key;
  activeMode = mode.key;
  const trace = {{
    type: "heatmap",
    x: item.lambdas,
    y: item.parameters,
    z: data.z,
    customdata: data.regionKeys,
    zmin: 0.5,
    zmax: data.colors.length + 0.5,
    colorscale: discreteColorscale(data.colors),
    showscale: true,
    colorbar: {{
      title: {{ text: isContractionMode ? "Υ/contract" : "Υ", side: "top", font: {{ family: "Times New Roman, Times, serif", size: 20, color: "#111111" }} }},
      tickmode: "array",
      tickvals: [],
      ticks: "",
      showticklabels: false,
      outlinecolor: "#111111",
      outlinewidth: 1.2,
      thickness: 18,
      len: 0.86
    }},
    hovertemplate: `λ=%{{x:.3f}}<br>c=%{{y:.3f}}<br>region=%{{customdata}}<br>phase=%{{z}}<extra></extra>`
  }};
  const boundaryTrace = phaseBoundaryTrace(item.lambdas, item.parameters, data.z);
  const layout = {{
    title: {{
      text: htmlTransitionTitle(item, mode),
      x: 0.5,
      xanchor: "center",
      font: {{ family: "Times New Roman, Times, serif", size: 24, color: "#111111" }}
    }},
    margin: {{ l: 74, r: 72, t: 68, b: 72 }},
    paper_bgcolor: "#ffffff",
    plot_bgcolor: "#ffffff",
    font: {{ family: "Times New Roman, Times, serif", color: "#111111" }},
    xaxis: {{ title: {{ text: "λ", font: {{ size: 28, family: "Times New Roman, Times, serif" }} }}, range: [item.lambdas[0], item.lambdas[item.lambdas.length - 1]], mirror: true, ticks: "outside", tickfont: {{ size: 18 }}, linewidth: 2, linecolor: "#111111", showgrid: false, zeroline: false }},
    yaxis: {{ title: {{ text: "c", font: {{ size: 28, family: "Times New Roman, Times, serif" }} }}, range: [item.parameters[0], item.parameters[item.parameters.length - 1]], mirror: true, ticks: "outside", tickfont: {{ size: 18 }}, linewidth: 2, linecolor: "#111111", showgrid: false, zeroline: false }}
  }};
  Plotly.react("phaseMap", [trace, boundaryTrace], layout, {{ responsive: true, displaylogo: false }});
  const phaseMap = document.getElementById("phaseMap");
  phaseMap.removeAllListeners?.("plotly_click");
  phaseMap.on("plotly_click", event => {{
    const regionKey = event.points?.[0]?.customdata;
    if (regionKey) selectRegion(regionKey);
  }});
  const firstRegion = Object.keys(regions).find(regionKey => regionKey.startsWith(`${{mode.key}}:${{item.key}}:`));
  if (firstRegion) selectRegion(firstRegion);
}}

function geometryTraces(region) {{
  const traces = [];
  if (region.surface && region.surface.triangle_count) {{
    traces.push({{
      type: "mesh3d",
      name: "level surface",
      x: region.surface.x,
      y: region.surface.y,
      z: region.surface.z,
      i: region.surface.i,
      j: region.surface.j,
      k: region.surface.k,
      color: region.color,
      opacity: 0.76,
      flatshading: false,
      lighting: {{ ambient: 0.48, diffuse: 0.72, specular: 0.18, roughness: 0.62, fresnel: 0.12 }},
      lightposition: {{ x: 120, y: 160, z: 220 }},
      hoverinfo: "skip"
    }});
  }}
  if (region.skeleton.lines.x.length) {{
    traces.push({{
      type: "scatter3d",
      name: "skeleton",
      x: region.skeleton.lines.x,
      y: region.skeleton.lines.y,
      z: region.skeleton.lines.z,
      mode: "lines",
      line: {{ color: "#111111", width: 8 }},
      hoverinfo: "skip"
    }});
  }}
  traces.push({{
    type: "scatter3d",
    name: "vertices",
    x: region.skeleton.nodes.x,
    y: region.skeleton.nodes.y,
    z: region.skeleton.nodes.z,
    mode: "markers",
    marker: {{ size: 7, color: "#d7191c", line: {{ color: "#ffffff", width: 1.2 }} }},
    hovertemplate: "core vertex<extra></extra>"
  }});
  return traces;
}}

function compactIdList(values, limit = 8) {{
  if (!values || values.length === 0) return "none";
  if (values.length <= limit) return values.join(", ");
  return `${{values.slice(0, limit).join(", ")}} +${{values.length - limit}} more`;
}}

function topologyStatusLine(topology) {{
  if (!topology) return "Solid topology: unavailable";
  const faces = topology.boundary_faces && topology.boundary_faces.length ? topology.boundary_faces.join(", ") : "none";
  const boxText = topology.touches_boundary ? "touches sampling box" : "inside sampling box";
  const surfaceText = topology.surface_is_closed ? "watertight surface" : `open surface, ${{topology.surface_open_edges}} open edges`;
  return `Solid topology: ${{surfaceText}}, ${{boxText}}; voxel handle rank=${{topology.handle_rank}}, components=${{topology.components}}, Euler=${{topology.euler_characteristic}}; box faces ${{faces}}`;
}}

function yamadaDisplay(polynomial) {{
  return String(polynomial || "");
}}

function selectRegion(regionKey) {{
  const region = regions[regionKey];
  if (!region) return;
  const layout = {{
    title: {{
      text: `${{region.title}}: ${{region.modeTitle}}, region ${{region.regionId}}, phase ${{region.phaseId}}`,
      x: 0,
      xanchor: "left",
      font: {{ family: "Times New Roman, serif", size: 21, color: "#111111" }}
    }},
    margin: {{ l: 0, r: 0, t: 50, b: 0 }},
    paper_bgcolor: "#ffffff",
    plot_bgcolor: "#ffffff",
    showlegend: true,
    legend: {{ x: 0.02, y: 0.98, bgcolor: "rgba(255,255,255,0.75)" }},
    scene: {{
      xaxis: {{ title: "x", showbackground: false, gridcolor: "#dddddd", zerolinecolor: "#cccccc" }},
      yaxis: {{ title: "y", showbackground: false, gridcolor: "#dddddd", zerolinecolor: "#cccccc" }},
      zaxis: {{ title: "z", showbackground: false, gridcolor: "#dddddd", zerolinecolor: "#cccccc" }},
      aspectmode: "data",
      camera: {{ eye: {{ x: 1.45, y: 1.55, z: 1.15 }} }}
    }}
  }};
  const skeletonKind = region.skeleton.edge_count ? "edge skeleton" : "vertex-only skeleton";
  Plotly.react("geometry", geometryTraces(region), layout, {{ responsive: true, displaylogo: false }});
  document.getElementById("status").innerHTML = `
    <strong>${{region.title}} - ${{region.modeTitle}}</strong>
    <span>region ${{region.regionId}}, phase ${{region.phaseId}}; representative λ=${{region.lambda}}, c=${{region.parameterValue}}; ${{region.source}}, ${{region.coreMode}}; ${{skeletonKind}}; core V=${{region.nodes}}, E=${{region.edges}}, β=${{region.cycleRank}}; surface triangles=${{region.surface.triangle_count}}, MC step=${{region.surface.marching_cubes_step_size || 1}}</span>
    <span>endpoints: ${{region.startName}} → ${{region.endName}}</span>
    <span>finite domain: ${{region.compactDomain.title}}; ${{region.compactDomain.formula}}</span>
    <span>${{topologyStatusLine(region.topology)}}; volume fraction=${{region.topology.volume_fraction}}; skeleton voxels=${{region.topology.skeleton_voxels}}</span>
    <span>Classic regions: ${{compactIdList(region.classicRegionIds)}}; Classic phases: ${{compactIdList(region.classicPhaseIds)}}; cells=${{region.cellCount}}</span>
    <span class="poly">Yamada/signature: ${{yamadaDisplay(region.polynomial)}}</span>
  `;
}}

function syncButtons() {{
  document.querySelectorAll("#modeButtons button").forEach(button => {{
    button.setAttribute("aria-pressed", button.dataset.key === activeMode ? "true" : "false");
  }});
  document.querySelectorAll("#transitionButtons button").forEach(button => {{
    button.setAttribute("aria-pressed", button.dataset.key === activeTransition ? "true" : "false");
  }});
}}

function buildButtons() {{
  const modeHolder = document.getElementById("modeButtons");
  modeHolder.innerHTML = "";
  modes.forEach(mode => {{
    const button = document.createElement("button");
    button.type = "button";
    button.dataset.key = mode.key;
    button.textContent = mode.title;
    button.setAttribute("aria-pressed", mode.key === activeMode ? "true" : "false");
    button.addEventListener("click", () => {{
      activeMode = mode.key;
      syncButtons();
      plotPhaseMap(activeTransition, activeMode);
    }});
    modeHolder.appendChild(button);
  }});

  const transitionHolder = document.getElementById("transitionButtons");
  transitionHolder.innerHTML = "";
  transitions.forEach(item => {{
    const button = document.createElement("button");
    button.type = "button";
    button.dataset.key = item.key;
    button.textContent = item.title;
    button.setAttribute("aria-pressed", item.key === activeTransition ? "true" : "false");
    button.addEventListener("click", () => {{
      activeTransition = item.key;
      syncButtons();
      plotPhaseMap(activeTransition, activeMode);
    }});
    transitionHolder.appendChild(button);
  }});
}}

buildButtons();
plotPhaseMap(activeTransition, activeMode);
window.addEventListener("resize", () => {{
  Plotly.Plots.resize(document.getElementById("phaseMap"));
  Plotly.Plots.resize(document.getElementById("geometry"));
}});
</script>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", type=Path, default=DEFAULT_SCAN_DIR)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("_build/new_phase_maps/tpms_region_geometry.html"),
    )
    parser.add_argument(
        "--scan-script",
        type=Path,
        default=None,
        help=(
            "optional custom Python implementation; by default use the installed compact-TPMS engine"
        ),
    )
    parser.add_argument("--max-surface-faces", type=int, default=2600)
    parser.add_argument("--max-contraction-states", type=int, default=8000)
    parser.add_argument(
        "--max-contraction-depth",
        type=int,
        default=4,
        help="largest node/edge drop allowed in one representative contraction test",
    )
    parser.add_argument(
        "--max-contraction-edges",
        type=int,
        default=24,
        help="largest representative core edge count eligible for contraction tests",
    )
    parser.add_argument(
        "--stable-min-cells",
        type=int,
        default=1,
        help="add a display-only stable mode by reassigning smaller connected regions",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = build_payload(
        args.scan_dir,
        scan_script=args.scan_script,
        max_surface_faces=int(args.max_surface_faces),
        max_contraction_states=int(args.max_contraction_states),
        max_contraction_depth=int(args.max_contraction_depth),
        max_contraction_edges=int(args.max_contraction_edges),
        stable_min_cells=int(args.stable_min_cells),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html_document(payload), encoding="utf-8")
    print(f"wrote: {args.output}")
    print(f"transitions: {len(payload['transitions'])}")
    print(f"regions: {len(payload['regions'])}")
    for transition in payload["transitions"]:
        classic = transition["modes"]["classic"]
        contraction = transition["modes"]["contraction"]
        print(
            f"{transition['key']}: classic {classic['classes']} classes/"
            f"{classic['components']} regions; contraction {contraction['classes']} "
            f"classes/{contraction['components']} regions"
        )


if __name__ == "__main__":
    main()
