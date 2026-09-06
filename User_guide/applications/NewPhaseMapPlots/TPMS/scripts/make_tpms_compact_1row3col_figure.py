#!/usr/bin/env python3
"""Render the compact TPMS phase maps as a publication PDF."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
from collections import Counter, deque
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/codex-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/codex-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap, to_hex
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes


DEFAULT_SCAN_DIR = Path("tmp/tpms_compact_scaffold_phase_maps_c0_03_dense")
DEFAULT_HTML = (
    DEFAULT_SCAN_DIR
    / "tpms_compact_c0_03_dense_stable_yamada_plotly_region_geometry.html"
)
DEFAULT_OUTPUT = (
    DEFAULT_SCAN_DIR
    / "tpms_compact_c0_03_stable_up_to_contraction_2row_cba.pdf"
)
DEFAULT_TRANSITION_ORDER = (
    "schwarz_p_to_diamond",
    "gyroid_to_schwarz_p",
    "gyroid_to_diamond",
)
PANEL_LABELS = ("(a)", "(b)", "(c)")
PANEL_C_TRANSITION_KEY = "gyroid_to_diamond"
PANEL_C_SLICE_VIEW_ANGLES = (
    (26.0, -64.0),
    (20.0, -30.0),
    (28.0, -48.0),
    (18.0, -76.0),
    (30.0, -24.0),
)
FIELD_KEYS_BY_TRANSITION = {
    "gyroid_to_diamond": ("gyroid", "diamond"),
    "gyroid_to_schwarz_p": ("gyroid", "schwarz_p"),
    "schwarz_p_to_diamond": ("schwarz_p", "diamond"),
}

PALETTE = (
    "#fdae1b",
    "#df2027",
    "#2c7fb8",
    "#9467bd",
    "#2ca02c",
    "#8c564b",
    "#17becf",
    "#bcbd22",
    "#e377c2",
    "#7f7f7f",
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
    "#666666",
    "#1f78b4",
    "#b2df8a",
    "#fb9a99",
    "#cab2d6",
    "#ffff99",
    "#6a3d9a",
    "#b15928",
    "#a6cee3",
    "#33a02c",
    "#fb8072",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
)


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
            for region in connected_regions(stable):
                if int(region["phase_id"]) != phase_id:
                    continue
                component = region["cells"]
                if len(component) >= min_cells:
                    continue
                neighbor_counts: Counter[int] = Counter()
                for row, col in component:
                    for delta_row, delta_col in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        next_row = row + delta_row
                        next_col = col + delta_col
                        if 0 <= next_row < stable.shape[0] and 0 <= next_col < stable.shape[1]:
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


def cell_edges(values: list[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if len(arr) == 1:
        return np.asarray([arr[0] - 0.5, arr[0] + 0.5], dtype=float)
    edges = np.empty(len(arr) + 1, dtype=float)
    edges[0] = arr[0] - 0.5 * (arr[1] - arr[0])
    edges[1:-1] = 0.5 * (arr[:-1] + arr[1:])
    edges[-1] = arr[-1] + 0.5 * (arr[-1] - arr[-2])
    return edges


def phase_boundaries(
    ax: plt.Axes,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    labels: np.ndarray,
) -> None:
    rows, cols = labels.shape
    for row in range(rows):
        for col in range(cols - 1):
            if int(labels[row, col]) == int(labels[row, col + 1]):
                continue
            x = x_edges[col + 1]
            ax.plot(
                [x, x],
                [y_edges[row], y_edges[row + 1]],
                color="black",
                lw=1.35,
                solid_capstyle="butt",
                zorder=4,
            )
    for row in range(rows - 1):
        for col in range(cols):
            if int(labels[row, col]) == int(labels[row + 1, col]):
                continue
            y = y_edges[row + 1]
            ax.plot(
                [x_edges[col], x_edges[col + 1]],
                [y, y],
                color="black",
                lw=1.35,
                solid_capstyle="butt",
                zorder=4,
            )


def distinct_color_list(count: int) -> list[str]:
    colors: list[str] = []
    seen: set[str] = set()
    for color in PALETTE:
        key = color.lower()
        if key in seen:
            continue
        seen.add(key)
        colors.append(color)
    for name in ("tab20", "tab20b", "tab20c", "Set3", "Dark2", "Paired", "Accent"):
        if len(colors) >= count:
            break
        cmap = plt.get_cmap(name)
        for index in range(cmap.N):
            color = to_hex(cmap(index), keep_alpha=False)
            key = color.lower()
            if key in seen:
                continue
            seen.add(key)
            colors.append(color)
            if len(colors) >= count:
                break
    if len(colors) < count:
        hsv = plt.get_cmap("hsv")
        for index in range(count - len(colors)):
            color = to_hex(hsv(index / max(1, count - len(colors))), keep_alpha=False)
            key = color.lower()
            if key in seen:
                continue
            seen.add(key)
            colors.append(color)
    return colors[:count]


def title_for_transition(key: str) -> str:
    titles = {
        "gyroid_to_diamond": r"$(1-\lambda)F_{\mathrm{G}}+\lambda F_{\mathrm{D}}$",
        "gyroid_to_schwarz_p": r"$(1-\lambda)F_{\mathrm{G}}+\lambda F_{\mathrm{P}}$",
        "schwarz_p_to_diamond": r"$(1-\lambda)F_{\mathrm{P}}+\lambda F_{\mathrm{D}}$",
    }
    return titles.get(key, key.replace("_", r"\_"))


def field_values(
    key: str,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    if key == "gyroid":
        return (
            np.sin(x) * np.cos(y)
            + np.sin(y) * np.cos(z)
            + np.sin(z) * np.cos(x)
        )
    if key == "schwarz_p":
        return np.cos(x) + np.cos(y) + np.cos(z)
    if key == "diamond":
        return np.cos(x) * np.cos(y) * np.cos(z) - np.sin(x) * np.sin(y) * np.sin(z)
    raise ValueError(f"unsupported TPMS field key: {key}")


def compact_values_for_geometry(metadata: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    span = np.asarray(metadata["span"], dtype=float)
    dimension = int(metadata["dimension"])
    axes = [
        np.linspace(float(lo), float(hi), dimension, dtype=np.float64)
        for lo, hi in span
    ]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    lam = float(metadata["lambda"])
    threshold_c = float(metadata["threshold_c"])
    start = field_values(str(metadata["start_field"]), x, y, z)
    end = field_values(str(metadata["end_field"]), x, y, z)
    raw = (1.0 - lam) * start + lam * end
    half_width = min(abs(float(lo)) for lo, _hi in span)
    radius = 0.72 * half_width
    domain = np.sqrt(x * x + y * y + z * z) / radius - 1.0
    values = np.maximum(raw - threshold_c, domain)
    spacing = np.diff(span, axis=1).squeeze() / (dimension - 1)
    origin = span[:, 0]
    return values, spacing, origin


def representative_geometry_path(
    geometry_dir: Path,
    transition_key: str,
    lambda_index: int,
    c_index: int,
) -> Path:
    return geometry_dir / f"{transition_key}_lambda{lambda_index:03d}_c{c_index:03d}.json.gz"


def load_geometry_metadata(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def slice_representatives(
    transition: dict[str, Any],
    stable_original_labels: np.ndarray,
    display_labels: np.ndarray,
    *,
    target_c: float,
    geometry_dir: Path,
) -> dict[str, Any]:
    parameters = [float(value) for value in transition["parameters"]]
    lambdas = [float(value) for value in transition["lambdas"]]
    row_index = min(range(len(parameters)), key=lambda index: abs(parameters[index] - target_c))
    row_values = [int(value) for value in stable_original_labels[row_index]]
    row_display_values = [int(value) for value in display_labels[row_index]]

    runs: list[dict[str, Any]] = []
    start_index = 0
    current_value = row_values[0]
    for col_index, value in enumerate(row_values[1:], start=1):
        if value == current_value:
            continue
        end_index = col_index - 1
        rep_index = (start_index + end_index) // 2
        runs.append(
            {
                "original_phase_id": int(current_value),
                "display_label": int(row_display_values[rep_index]),
                "lambda_start": float(lambdas[start_index]),
                "lambda_end": float(lambdas[end_index]),
                "lambda_index": int(rep_index),
                "lambda": float(lambdas[rep_index]),
                "c_index": int(row_index),
                "c": float(parameters[row_index]),
                "geometry_path": str(
                    representative_geometry_path(
                        geometry_dir,
                        transition["key"],
                        rep_index,
                        row_index,
                    )
                ),
            }
        )
        start_index = col_index
        current_value = value

    end_index = len(row_values) - 1
    rep_index = (start_index + end_index) // 2
    runs.append(
        {
            "original_phase_id": int(current_value),
            "display_label": int(row_display_values[rep_index]),
            "lambda_start": float(lambdas[start_index]),
            "lambda_end": float(lambdas[end_index]),
            "lambda_index": int(rep_index),
            "lambda": float(lambdas[rep_index]),
            "c_index": int(row_index),
            "c": float(parameters[row_index]),
            "geometry_path": str(
                representative_geometry_path(
                    geometry_dir,
                    transition["key"],
                    rep_index,
                    row_index,
                )
            ),
        }
    )
    return {
        "target_c": float(target_c),
        "displayed_c": float(parameters[row_index]),
        "row_index": int(row_index),
        "runs": runs,
    }


def plot_handlebody_and_skeleton(
    ax: plt.Axes,
    metadata: dict[str, Any],
    *,
    surface_color: str,
    max_faces: int,
    view_angle: tuple[float, float] = (22.0, -38.0),
) -> None:
    values, spacing, origin = compact_values_for_geometry(metadata)
    verts, faces, _normals, _sampled_values = marching_cubes(
        values,
        level=0.0,
        spacing=tuple(float(value) for value in spacing),
    )
    verts = verts + origin[None, :]
    if len(faces) > max_faces:
        face_indices = np.linspace(0, len(faces) - 1, int(max_faces), dtype=int)
        faces = faces[face_indices]
    mesh = Poly3DCollection(
        verts[faces],
        facecolor=surface_color,
        edgecolor=(1.0, 1.0, 1.0, 0.0),
        linewidths=0.0,
        alpha=0.26,
    )
    mesh.set_rasterized(True)
    ax.add_collection3d(mesh)

    graph = metadata.get("graph", {})
    for edge in graph.get("edges", []):
        points = np.asarray(edge.get("points_coord", []), dtype=float)
        if points.ndim != 2 or points.shape[0] < 2:
            continue
        ax.plot3D(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            color="#111827",
            linewidth=2.2,
            solid_capstyle="round",
            zorder=10,
        )
    node_points = [
        node.get("coord")
        for node in graph.get("nodes", [])
        if isinstance(node.get("coord"), list) and len(node.get("coord")) == 3
    ]
    if node_points:
        nodes = np.asarray(node_points, dtype=float)
        ax.scatter(
            nodes[:, 0],
            nodes[:, 1],
            nodes[:, 2],
            s=56,
            color="#d7191c",
            edgecolors="#ffffff",
            linewidths=0.55,
            depthshade=False,
            zorder=11,
        )

    radius = 5.25
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_zlim(-radius, radius)
    try:
        ax.set_box_aspect((1, 1, 1), zoom=1.65)
    except TypeError:
        ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=float(view_angle[0]), azim=float(view_angle[1]))
    ax.set_axis_off()


def render_figure(
    payload: dict[str, Any],
    *,
    input_html: Path,
    output_pdf: Path,
    stable_min_cells: int,
    transition_order: tuple[str, ...],
    include_panel_c_slice: bool,
    geometry_dir: Path,
    slice_c: float,
    surface_max_faces: int,
) -> dict[str, Any]:
    transitions_by_key = {transition["key"]: transition for transition in payload["transitions"]}
    missing = [key for key in transition_order if key not in transitions_by_key]
    if missing:
        raise ValueError(f"transition keys not found in payload: {missing}")
    transitions = [transitions_by_key[key] for key in transition_order]
    grids: dict[str, np.ndarray] = {}
    diagnostics: dict[str, Any] = {}

    for transition in transitions:
        contraction = np.asarray(transition["modes"]["contraction"]["z"], dtype=int)
        stable, changed_cells = stable_labels(
            contraction,
            min_cells=int(stable_min_cells),
        )
        regions = connected_regions(stable)
        grids[transition["key"]] = stable
        diagnostics[transition["key"]] = {
            "raw_contraction_classes": int(transition["modes"]["contraction"]["classes"]),
            "raw_contraction_regions": int(transition["modes"]["contraction"]["components"]),
            "stable_up_to_contraction_classes": int(len(set(int(v) for v in stable.ravel()))),
            "stable_up_to_contraction_regions": int(len(regions)),
            "stable_min_cells": int(stable_min_cells),
            "reassigned_cells": int(changed_cells),
            "single_cell_regions": sum(1 for region in regions if len(region["cells"]) == 1),
            "small_regions_le_3_cells": sum(1 for region in regions if len(region["cells"]) <= 3),
        }

    used_phase_ids = sorted(
        {int(value) for transition in transitions for value in grids[transition["key"]].ravel()}
    )
    display_label_by_phase_id = {
        phase_id: index + 1 for index, phase_id in enumerate(used_phase_ids)
    }
    display_grids = {
        key: np.vectorize(display_label_by_phase_id.__getitem__)(value).astype(int)
        for key, value in grids.items()
    }
    max_classes = len(used_phase_ids)
    original_phase_colors = distinct_color_list(max(used_phase_ids))
    colors = [original_phase_colors[phase_id - 1] for phase_id in used_phase_ids]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(0.5, max_classes + 1.5), cmap.N)

    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "text.usetex": False,
            "figure.dpi": 220,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.linewidth": 1.5,
            "axes.edgecolor": "black",
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(
        figsize=(15.3, 8.5 if include_panel_c_slice else 7.0),
        constrained_layout=True,
        facecolor="white",
    )
    if include_panel_c_slice:
        grid = fig.add_gridspec(
            2,
            6,
            height_ratios=[0.88, 1.35],
            hspace=0.0,
        )
        slots = [
            grid[0, 0:2],
            grid[0, 2:4],
            grid[0, 4:6],
        ]
    else:
        grid = fig.add_gridspec(2, 4)
        slots = [
            grid[0, 0:2],
            grid[0, 2:4],
            grid[1, 1:3],
        ]
    mappable = None
    top_right_ax = None
    panel_c_ax = None
    panel_c_transition = None
    panel_c_representatives = None
    top_axes: list[plt.Axes] = []
    geo_axes: list[plt.Axes] = []
    arrow_axes: list[plt.Axes] = []
    geo_interval_labels: list[str] = []
    box_ax: plt.Axes | None = None

    for slot, transition, panel_label in zip(slots, transitions, PANEL_LABELS):
        ax = fig.add_subplot(slot)
        top_axes.append(ax)
        key = transition["key"]
        labels = display_grids[key]
        lambdas = [float(value) for value in transition["lambdas"]]
        parameters = [float(value) for value in transition["parameters"]]
        x_edges = cell_edges(lambdas)
        y_edges = cell_edges(parameters)

        mappable = ax.pcolormesh(
            x_edges,
            y_edges,
            labels,
            cmap=cmap,
            norm=norm,
            shading="flat",
            edgecolors=(1.0, 1.0, 1.0, 0.16),
            linewidth=0.28,
            antialiased=False,
        )
        phase_boundaries(ax, x_edges, y_edges, labels)
        if include_panel_c_slice and key == PANEL_C_TRANSITION_KEY:
            panel_c_representatives = slice_representatives(
                transition,
                grids[key],
                labels,
                target_c=float(slice_c),
                geometry_dir=geometry_dir,
            )
            displayed_c = float(panel_c_representatives["displayed_c"])
            ax.axhline(
                displayed_c,
                color="#111827",
                linewidth=1.1,
                linestyle=(0, (3.2, 2.0)),
                zorder=6,
            )
            marker_lambdas = [
                float(run["lambda"]) for run in panel_c_representatives["runs"]
            ]
            ax.scatter(
                marker_lambdas,
                [displayed_c] * len(marker_lambdas),
                s=18,
                color="#111827",
                edgecolor="white",
                linewidth=0.45,
                zorder=7,
            )
        ax.set_xlim(float(lambdas[0]), float(lambdas[-1]))
        ax.set_ylim(float(parameters[0]), float(parameters[-1]))
        title_size = 19 if include_panel_c_slice else 16
        label_size = 25 if include_panel_c_slice else 22
        tick_size = 18 if include_panel_c_slice else 16
        ax.set_title(title_for_transition(key), fontsize=title_size, fontweight="semibold", pad=10)
        ax.set_xlabel(r"$\lambda$", fontsize=label_size, labelpad=5)
        ax.set_ylabel(r"$c$", fontsize=label_size, labelpad=6)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.set_yticks([0.0, 0.1, 0.2, 0.3])
        ax.tick_params(axis="both", which="major", labelsize=tick_size, width=1.5, length=5)
        ax.tick_params(axis="both", which="minor", width=0.9, length=3)
        ax.minorticks_on()
        for spine in ax.spines.values():
            spine.set_linewidth(1.75)
            spine.set_color("black")
        ax.text(
            -0.22,
            1.05,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=23 if include_panel_c_slice else 20,
            fontweight="bold",
            family="serif",
            color="#111827",
            clip_on=False,
        )
        if panel_label == "(b)":
            top_right_ax = ax
        if key == PANEL_C_TRANSITION_KEY:
            panel_c_ax = ax
            panel_c_transition = transition

    if mappable is not None and top_right_ax is not None and not include_panel_c_slice:
        cax = top_right_ax.inset_axes([1.045, 0.0, 0.04, 1.0], transform=top_right_ax.transAxes)
        cbar = fig.colorbar(
            mappable,
            cax=cax,
            boundaries=np.arange(0.5, max_classes + 1.5, 1.0),
            ticks=[],
            spacing="proportional",
        )
        cbar.ax.set_title(r"$\Upsilon$", fontsize=18, fontweight="semibold", pad=7)
        cbar.set_ticks([])
        cbar.ax.tick_params(width=0, length=0)
        cbar.outline.set_linewidth(1.2)

    if include_panel_c_slice:
        if panel_c_ax is None or panel_c_transition is None or panel_c_representatives is None:
            raise ValueError(f"could not identify panel-c transition {PANEL_C_TRANSITION_KEY}")
        box_ax = fig.add_subplot(grid[1, :])
        box_ax.set_in_layout(False)
        box_ax.set_zorder(30)
        box_ax.patch.set_alpha(0.0)
        box_ax.set_xticks([])
        box_ax.set_yticks([])
        box_ax.set_xlim(0.0, 1.0)
        box_ax.set_ylim(0.0, 1.0)
        for spine in box_ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.45)
            spine.set_color("#111827")
        box_ax.text(
            0.0,
            0.965,
            "(d)",
            ha="left",
            va="top",
            transform=box_ax.transAxes,
            fontsize=26,
            fontweight="bold",
            family="serif",
            color="#111827",
            zorder=31,
        )
        geo_grid = grid[1, :].subgridspec(
            1,
            11,
            width_ratios=[
                0.62,
                1.46,
                0.055,
                1.46,
                0.055,
                1.46,
                0.055,
                1.46,
                0.055,
                1.46,
                0.62,
            ],
            wspace=0.005,
        )
        runs = list(panel_c_representatives["runs"])
        if len(runs) != 5:
            raise ValueError(f"expected five panel-c representatives, found {len(runs)}")
        for index, run in enumerate(runs):
            geo_col = 1 + index * 2
            ax_geo = fig.add_subplot(geo_grid[0, geo_col], projection="3d")
            ax_geo.set_in_layout(False)
            ax_geo.set_zorder(10)
            ax_geo.patch.set_alpha(0.0)
            geo_axes.append(ax_geo)
            geometry_path = Path(str(run["geometry_path"]))
            if not geometry_path.exists():
                raise FileNotFoundError(geometry_path)
            metadata = load_geometry_metadata(geometry_path)
            color = colors[int(run["display_label"]) - 1]
            plot_handlebody_and_skeleton(
                ax_geo,
                metadata,
                surface_color=color,
                max_faces=int(surface_max_faces),
                view_angle=PANEL_C_SLICE_VIEW_ANGLES[index],
            )
            if abs(float(run["lambda_start"]) - float(run["lambda_end"])) < 1e-12:
                interval = rf"$\lambda={run['lambda']:.2f}$"
            else:
                interval = (
                    rf"${run['lambda_start']:.2f}\leq\lambda\leq"
                    rf"{run['lambda_end']:.2f}$"
                )
            geo_interval_labels.append(interval)
            if index < len(runs) - 1:
                ax_arrow = fig.add_subplot(geo_grid[0, geo_col + 1])
                ax_arrow.set_in_layout(False)
                ax_arrow.set_zorder(10)
                ax_arrow.patch.set_alpha(0.0)
                arrow_axes.append(ax_arrow)
                ax_arrow.set_axis_off()
                ax_arrow.text(
                    0.5,
                    0.54,
                    r"$\rightarrow$",
                    ha="center",
                    va="center",
                    fontsize=34,
                    color="#111827",
                    transform=ax_arrow.transAxes,
                )

        if box_ax is not None:
            fig.canvas.draw()
            fig.set_constrained_layout(False)
            renderer = fig.canvas.get_renderer()
            top_extent = Bbox.union([ax.get_tightbbox(renderer) for ax in top_axes])
            top_extent = top_extent.transformed(fig.transFigure.inverted())
            current_box = box_ax.get_position()
            box_x0 = top_extent.x0
            box_x1 = top_extent.x1
            box_y0 = current_box.y0 + 0.004
            box_y1 = min(current_box.y1 + 0.012, top_extent.y0 - 0.006)
            if box_y1 <= box_y0:
                box_y1 = current_box.y1
            box_width = box_x1 - box_x0
            box_height = box_y1 - box_y0
            box_ax.set_position([box_x0, box_y0, box_width, box_height])

            inner_x0 = box_x0 + 0.045 * box_width
            inner_x1 = box_x1 - 0.035 * box_width
            object_width = 0.142 * box_width
            centers = np.linspace(
                inner_x0 + 0.5 * object_width,
                inner_x1 - 0.5 * object_width,
                len(geo_axes),
            )
            object_y0 = box_y0 + 0.055 * box_height
            object_height = 0.760 * box_height
            for ax_geo, center in zip(geo_axes, centers):
                ax_geo.set_position(
                    [
                        float(center - 0.5 * object_width),
                        float(object_y0),
                        float(object_width),
                        float(object_height),
                    ]
                )
            for label, center in zip(geo_interval_labels, centers):
                box_ax.text(
                    float((center - box_x0) / box_width),
                    0.925,
                    label,
                    ha="center",
                    va="top",
                    transform=box_ax.transAxes,
                    fontsize=19,
                    color="#111827",
                    zorder=32,
                )
            arrow_width = 0.030 * box_width
            arrow_height = 0.26 * box_height
            arrow_y0 = box_y0 + 0.38 * box_height
            for index, ax_arrow in enumerate(arrow_axes):
                arrow_center = 0.5 * (centers[index] + centers[index + 1])
                ax_arrow.set_position(
                    [
                        float(arrow_center - 0.5 * arrow_width),
                        float(arrow_y0),
                        float(arrow_width),
                        float(arrow_height),
                    ]
                )

    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output_png = output_pdf.with_suffix(".png")
    source_json = output_pdf.with_name(output_pdf.stem + "_source_data.json")
    fig.savefig(output_pdf, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(output_png, dpi=450, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)

    source = {
        "figure": str(output_pdf),
        "mode": "stable up to contraction",
        "layout": "2row_phase_maps_plus_panel_c_slice" if include_panel_c_slice else "2row_2_1",
        "transition_order": list(transition_order),
        "panel_labels_by_position": list(PANEL_LABELS),
        "panel_c_slice": panel_c_representatives,
        "stable_min_cells": int(stable_min_cells),
        "surface_max_faces": int(surface_max_faces),
        "input_html": str(input_html),
        "scan_dir": payload.get("scanDir"),
        "scan_summary": payload.get("summary"),
        "transitions": diagnostics,
        "original_phase_id_to_display_label": {
            str(phase_id): display_label
            for phase_id, display_label in display_label_by_phase_id.items()
        },
        "display_colors": colors,
        "grids_original_phase_ids": {
            key: value.astype(int).tolist() for key, value in grids.items()
        },
        "grids_display_labels": {
            key: value.astype(int).tolist() for key, value in display_grids.items()
        },
        "style_reference": (
            "07_hamiltonian_yamada_phase_maps_5transition_overview_"
            "up_to_contraction_3row_2_2_1.pdf"
        ),
    }
    source_json.write_text(json.dumps(source, indent=2), encoding="utf-8")
    return {
        "pdf": str(output_pdf),
        "png": str(output_png),
        "source_json": str(source_json),
        "diagnostics": diagnostics,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, default=DEFAULT_HTML)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stable-min-cells", type=int, default=4)
    parser.add_argument("--include-panel-c-slice", action="store_true")
    parser.add_argument("--geometry-dir", type=Path, default=DEFAULT_SCAN_DIR / "geometry")
    parser.add_argument("--slice-c", type=float, default=0.1)
    parser.add_argument("--surface-max-faces", type=int, default=100000)
    parser.add_argument(
        "--transition-order",
        nargs="+",
        default=list(DEFAULT_TRANSITION_ORDER),
        help="Transition keys to plot in order.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = load_html_payload(args.html)
    result = render_figure(
        payload,
        input_html=args.html,
        output_pdf=args.output,
        stable_min_cells=int(args.stable_min_cells),
        transition_order=tuple(args.transition_order),
        include_panel_c_slice=bool(args.include_panel_c_slice),
        geometry_dir=args.geometry_dir,
        slice_c=float(args.slice_c),
        surface_max_faces=int(args.surface_max_faces),
    )
    print("wrote:", result["pdf"])
    print("wrote:", result["png"])
    print("wrote:", result["source_json"])
    for key, value in result["diagnostics"].items():
        print(
            f"{key}: {value['raw_contraction_regions']} contraction regions -> "
            f"{value['stable_up_to_contraction_regions']} stable regions; "
            f"classes={value['stable_up_to_contraction_classes']}, "
            f"reassigned={value['reassigned_cells']}"
        )


if __name__ == "__main__":
    main()
