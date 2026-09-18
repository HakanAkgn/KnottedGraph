#!/usr/bin/env python3
"""Validate compact TPMS phase-map outputs and clickable HTML payload."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

import numpy as np
from knotted_graph.applications.phase_map_examples._runtime import (
    geometry_directory,
    resolve_geometry_path,
)


REFERENCE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_SCAN_DIR = REFERENCE_DIR / "data"
DEFAULT_HTML = (
    REFERENCE_DIR
    / "html"
    / "tpms_compact_c0_03_dense_stable_yamada_plotly_region_geometry.html"
)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_html_payload(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"const payload = (.*?);\nconst transitions", text)
    if match is None:
        raise ValueError(f"could not find embedded payload in {path}")
    return json.loads(match.group(1))


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
                for delta_row, delta_col in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    next_row = current_row + delta_row
                    next_col = current_col + delta_col
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


def stable_labels(
    labels: np.ndarray,
    *,
    min_cells: int,
) -> tuple[np.ndarray, int]:
    if min_cells <= 1:
        return labels.copy(), 0
    stable = labels.copy()
    changed_cells = 0
    for _ in range(3):
        changed = False
        for phase_id in sorted(set(int(value) for value in stable.ravel())):
            components = []
            seen = np.zeros(stable.shape, dtype=bool)
            for start_row in range(stable.shape[0]):
                for start_col in range(stable.shape[1]):
                    if (
                        seen[start_row, start_col]
                        or int(stable[start_row, start_col]) != phase_id
                    ):
                        continue
                    queue = deque([(start_row, start_col)])
                    seen[start_row, start_col] = True
                    component: list[tuple[int, int]] = []
                    while queue:
                        row, col = queue.popleft()
                        component.append((row, col))
                        for delta_row, delta_col in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                            next_row = row + delta_row
                            next_col = col + delta_col
                            if not (
                                0 <= next_row < stable.shape[0]
                                and 0 <= next_col < stable.shape[1]
                            ):
                                continue
                            if seen[next_row, next_col]:
                                continue
                            if int(stable[next_row, next_col]) != phase_id:
                                continue
                            seen[next_row, next_col] = True
                            queue.append((next_row, next_col))
                    components.append(component)
            for component in components:
                if len(component) >= min_cells:
                    continue
                neighbor_counts: Counter[int] = Counter()
                for row, col in component:
                    for delta_row, delta_col in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        next_row = row + delta_row
                        next_col = col + delta_col
                        if (
                            0 <= next_row < stable.shape[0]
                            and 0 <= next_col < stable.shape[1]
                        ):
                            neighbor = int(stable[next_row, next_col])
                            if neighbor != phase_id:
                                neighbor_counts[neighbor] += 1
                if not neighbor_counts:
                    continue
                replacement = neighbor_counts.most_common(1)[0][0]
                for row, col in component:
                    stable[row, col] = replacement
                    changed_cells += 1
                changed = True
        if not changed:
            break
    return stable, changed_cells


def records_by_family(records: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["family"])].append(record)
    return dict(grouped)


def recompute_raw_grid(
    family_records: list[dict[str, Any]],
    lambdas: list[float],
    thresholds: list[float],
) -> tuple[np.ndarray, dict[str, int]]:
    signatures: list[str] = []
    seen: set[str] = set()
    for record in family_records:
        signature = str(record["phase_signature"])
        if signature in seen:
            continue
        signatures.append(signature)
        seen.add(signature)
    signature_to_id = {
        signature: index + 1 for index, signature in enumerate(signatures)
    }
    lookup = {
        (
            round(float(record["threshold_c"]), 12),
            round(float(record["lam"]), 12),
        ): record
        for record in family_records
    }
    grid = np.zeros((len(thresholds), len(lambdas)), dtype=int)
    for row, threshold in enumerate(thresholds):
        for col, lam in enumerate(lambdas):
            record = lookup[(round(float(threshold), 12), round(float(lam), 12))]
            grid[row, col] = signature_to_id[str(record["phase_signature"])]
    return grid, signature_to_id


def phase_signature_stats(records: list[dict[str, Any]]) -> dict[str, Any]:
    exact_groups: dict[str, set[str]] = defaultdict(set)
    large_core_without_wl = 0
    gross_groups: dict[tuple[Any, ...], set[str]] = defaultdict(set)
    for record in records:
        source = str(record["source"])
        signature = str(record["phase_signature"])
        if source == "yamada":
            exact_groups[signature].add(str(record["polynomial"]))
        if source == "large-core" and "wl=" not in signature:
            large_core_without_wl += 1
        gross_key = (
            record["family"],
            source,
            record["nodes"],
            record["edges"],
            record["cycle_rank"],
            record["interior_components"],
            record["handle_rank"],
            tuple(record["degree_sequence"]),
        )
        gross_groups[gross_key].add(signature)
    exact_polynomial_conflicts = {
        key: sorted(values) for key, values in exact_groups.items() if len(values) > 1
    }
    graph_similar_but_signature_distinct = [
        {
            "family": key[0],
            "source": key[1],
            "nodes": key[2],
            "edges": key[3],
            "cycle_rank": key[4],
            "interior_components": key[5],
            "handle_rank": key[6],
            "degree_sequence": list(key[7]),
            "distinct_signatures": len(values),
        }
        for key, values in gross_groups.items()
        if len(values) > 1
    ]
    graph_similar_but_signature_distinct.sort(
        key=lambda item: (
            str(item["family"]),
            -int(item["distinct_signatures"]),
            int(item["edges"]),
        )
    )
    return {
        "exact_yamada_polynomial_conflicts": exact_polynomial_conflicts,
        "large_core_records_missing_wl_hash": large_core_without_wl,
        "grossly_similar_graph_groups_split_by_signature": (
            graph_similar_but_signature_distinct[:20]
        ),
    }


def validate(args: argparse.Namespace) -> dict[str, Any]:
    scan_dir = Path(args.scan_dir)
    html_path = Path(args.html)
    records = load_json(scan_dir / "tpms_parameter_phase_map_records.json")
    source_data = load_json(scan_dir / "tpms_parameter_phase_map_source_data.json")
    summary = load_json(scan_dir / "tpms_parameter_phase_map_summary.json")
    payload = load_html_payload(html_path)

    failures: list[str] = []
    warnings: list[str] = []
    family_records = records_by_family(records)
    lambdas = [float(value) for value in source_data["lambdas"]]

    geometry_files = list(geometry_directory(scan_dir).glob("*.json.gz"))
    if len(geometry_files) != len(records):
        failures.append(
            f"geometry file count {len(geometry_files)} does not match record count {len(records)}"
        )
    for record in records:
        try:
            resolve_geometry_path(str(record["geometry_path"]), scan_dir)
        except FileNotFoundError as exc:
            failures.append(str(exc))
            break

    compactness = {
        "record_count": len(records),
        "box_boundary_touch_count": sum(
            int(record["touches_boundary"]) for record in records
        ),
        "closed_surface_count": sum(
            int(record["surface_is_closed"]) for record in records
        ),
        "open_surface_count": sum(
            int(not record["surface_is_closed"]) for record in records
        ),
        "surface_open_edge_max": max(
            int(record["surface_open_edges"]) for record in records
        ),
        "surface_nonmanifold_edge_max": max(
            int(record["surface_nonmanifold_edges"]) for record in records
        ),
        "empty_solid_count": sum(
            int(record["interior_voxels"] == 0) for record in records
        ),
    }
    if compactness["box_boundary_touch_count"]:
        failures.append("one or more compact solids touch the sampling box boundary")
    if compactness["open_surface_count"]:
        failures.append("one or more marching-cubes surfaces are not closed")
    if (
        compactness["surface_open_edge_max"]
        or compactness["surface_nonmanifold_edge_max"]
    ):
        failures.append("one or more surfaces have open or nonmanifold edges")
    if compactness["empty_solid_count"]:
        failures.append("one or more parameter cells generated an empty solid")

    endpoint_failures = []
    for endpoint in summary["endpoint_verification"]:
        expected = endpoint["expected_field"]
        correlation = float(endpoint["correlations"][expected])
        if float(endpoint["self_max_abs_error"]) > 1e-10:
            endpoint_failures.append(
                f"{endpoint['family']} lambda={endpoint['lambda']} field mismatch"
            )
        if correlation < 1.0 - 1e-10:
            endpoint_failures.append(
                f"{endpoint['family']} lambda={endpoint['lambda']} low correlation"
            )
        if endpoint["interior"]["touches_boundary"]:
            endpoint_failures.append(
                f"{endpoint['family']} lambda={endpoint['lambda']} touches box"
            )
        if not endpoint["surface"]["surface_is_closed"]:
            endpoint_failures.append(
                f"{endpoint['family']} lambda={endpoint['lambda']} open surface"
            )
    failures.extend(endpoint_failures)

    source_counts = dict(Counter(record["source"] for record in records))
    if source_counts != summary["source_counts"]:
        failures.append(
            f"summary source_counts {summary['source_counts']} != records {source_counts}"
        )

    family_checks: dict[str, Any] = {}
    html_transition_lookup = {
        transition["key"]: transition for transition in payload["transitions"]
    }
    for family_key, map_data in source_data["phase_maps"].items():
        thresholds = [float(value) for value in map_data["thresholds"]]
        computed_raw, signature_to_id = recompute_raw_grid(
            family_records[family_key],
            lambdas,
            thresholds,
        )
        raw_grid = np.asarray(map_data["raw_grid"], dtype=int)
        stable_grid = np.asarray(map_data["stable_grid"], dtype=int)
        if not np.array_equal(computed_raw, raw_grid):
            failures.append(f"{family_key}: raw grid does not match records")
        if not np.array_equal(raw_grid, stable_grid):
            warnings.append(f"{family_key}: stable grid differs from raw grid")

        transition = html_transition_lookup.get(family_key)
        if transition is None:
            failures.append(f"{family_key}: missing HTML transition")
            continue
        classic = transition["modes"]["classic"]
        html_grid = np.asarray(classic["z"], dtype=int)
        if not np.array_equal(stable_grid, html_grid):
            failures.append(
                f"{family_key}: HTML classic grid does not match source stable grid"
            )
        classic_regions = connected_regions(html_grid)
        if int(classic["components"]) != len(classic_regions):
            failures.append(f"{family_key}: HTML classic component count mismatch")
        region_key_grid = classic["regionKeys"]
        for region in classic_regions:
            keys = {region_key_grid[row][col] for row, col in region["cells"]}
            if len(keys) != 1:
                failures.append(
                    f"{family_key}: a connected classic region has multiple keys"
                )
                break
            key = next(iter(keys))
            payload_region = payload["regions"].get(key)
            if payload_region is None:
                failures.append(f"{family_key}: region key {key} missing from payload")
                break
            if int(payload_region["cellCount"]) != len(region["cells"]):
                failures.append(f"{family_key}: region {key} cell count mismatch")
                break

        stable_display_info = None
        if "stable" in transition["modes"]:
            stable_mode = transition["modes"]["stable"]
            stable_min_cells = int(stable_mode.get("stableMinCells", 1))
            expected_stable, changed_cells = stable_labels(
                raw_grid,
                min_cells=stable_min_cells,
            )
            html_stable = np.asarray(stable_mode["z"], dtype=int)
            if not np.array_equal(expected_stable, html_stable):
                failures.append(
                    f"{family_key}: HTML stable display grid does not match min-cells rule"
                )
            if int(stable_mode.get("stableReassignedCells", -1)) != int(changed_cells):
                failures.append(
                    f"{family_key}: HTML stable reassigned-cell count mismatch"
                )
            stable_regions = connected_regions(html_stable)
            if int(stable_mode["components"]) != len(stable_regions):
                failures.append(f"{family_key}: HTML stable component count mismatch")
            stable_region_key_grid = stable_mode["regionKeys"]
            for region in stable_regions:
                keys = {
                    stable_region_key_grid[row][col] for row, col in region["cells"]
                }
                if len(keys) != 1:
                    failures.append(
                        f"{family_key}: a connected stable region has multiple keys"
                    )
                    break
                key = next(iter(keys))
                payload_region = payload["regions"].get(key)
                if payload_region is None:
                    failures.append(
                        f"{family_key}: stable region key {key} missing from payload"
                    )
                    break
                if int(payload_region["cellCount"]) != len(region["cells"]):
                    failures.append(
                        f"{family_key}: stable region {key} cell count mismatch"
                    )
                    break
            stable_display_info = {
                "min_cells": stable_min_cells,
                "phase_count": int(
                    len(set(int(value) for value in html_stable.ravel()))
                ),
                "connected_regions": len(stable_regions),
                "reassigned_cells": int(changed_cells),
                "single_cell_regions": sum(
                    1 for region in stable_regions if len(region["cells"]) == 1
                ),
                "small_regions_le_3_cells": sum(
                    1 for region in stable_regions if len(region["cells"]) <= 3
                ),
            }

        family_checks[family_key] = {
            "records": len(family_records[family_key]),
            "raw_phase_count": int(raw_grid.max()),
            "stable_phase_count": int(
                len(set(int(value) for value in stable_grid.ravel()))
            ),
            "classic_connected_regions": len(classic_regions),
            "stable_display": stable_display_info,
            "html_contraction_classes": int(
                transition["modes"]["contraction"]["classes"]
            ),
            "html_contraction_regions": int(
                transition["modes"]["contraction"]["components"]
            ),
            "signature_count": len(signature_to_id),
            "source_counts": dict(
                Counter(record["source"] for record in family_records[family_key])
            ),
        }

    html_region_failures = []
    for key, region in payload["regions"].items():
        topology = region["topology"]
        if not topology["surface_is_closed"]:
            html_region_failures.append(f"{key}: representative surface not closed")
        if topology["touches_boundary"]:
            html_region_failures.append(f"{key}: representative touches box boundary")
        if int(region["surface"]["triangle_count"]) <= 0:
            html_region_failures.append(
                f"{key}: representative surface has no triangles"
            )
        if region["compactDomain"]["key"] == "none":
            html_region_failures.append(f"{key}: missing compact domain metadata")
    failures.extend(html_region_failures)

    report = {
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "warnings": warnings,
        "compactness": compactness,
        "endpoint_failure_count": len(endpoint_failures),
        "source_counts": source_counts,
        "summary_timing": summary["map_info"]["timing"],
        "html": {
            "path": str(html_path),
            "transition_count": len(payload["transitions"]),
            "region_count": len(payload["regions"]),
            "modes": [mode["key"] for mode in payload["modes"]],
            "summary": payload["summary"],
        },
        "families": family_checks,
        "signature_stats": phase_signature_stats(records),
    }
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", type=Path, default=DEFAULT_SCAN_DIR)
    parser.add_argument("--html", type=Path, default=DEFAULT_HTML)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("_build/new_phase_maps/tpms_validation.json"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = validate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"status: {report['status']}")
    print(f"wrote: {args.output}")
    print(f"failures: {len(report['failures'])}")
    print(f"warnings: {len(report['warnings'])}")
    print(f"compactness: {report['compactness']}")
    for family_key, family in report["families"].items():
        print(
            f"{family_key}: phases={family['raw_phase_count']}, "
            f"classic_regions={family['classic_connected_regions']}, "
            f"contraction_regions={family['html_contraction_regions']}"
        )
    if report["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
