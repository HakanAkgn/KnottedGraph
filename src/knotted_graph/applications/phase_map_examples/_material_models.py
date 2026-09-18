#!/usr/bin/env python3
"""Reproduce the legacy Multiband.ipynb material surfaces with current code.

The old notebook is used only as a reference for Hamiltonian formulas,
parameters, k-space spans, grid dimensions, band pairs, and energy cuts.
Surface sampling and contour extraction go through the current
installed MaterialFermiSurface workflow.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import plotly.graph_objects as go
import pyvista as pv
import sympy as sp
from plotly.subplots import make_subplots


from ._runtime import package_location

from knotted_graph.applications.material_surface import MaterialFermiSurface  # noqa: E402


kx, ky, kz = sp.symbols("k_x k_y k_z", real=True)
K_SYMBOLS = (kx, ky, kz)


def H_D6_sympy(params: dict[str, float]) -> sp.Matrix:
    kplus = kx + sp.I * ky
    kminus = kx - sp.I * ky
    k2xy = kx**2 + ky**2

    eps1 = params["E1"] + params["A1"] * k2xy + params["B1"] * kz**2
    eps2 = (
        params["E2"]
        + params["A2"] * k2xy
        + params["B2"] * kz**2
        + params["L"] * k2xy**2
        + params["M"] * (kplus**6 + kminus**6)
    )

    h12 = params["C"] * kminus**2 + params["F"] * kplus**4
    h13 = params["D"] * kminus * kz
    h23 = params["D"] * kplus * kz

    return sp.Matrix(
        [
            [eps1, sp.expand(h12), sp.expand(h13)],
            [sp.conjugate(h12), eps1, sp.expand(h23)],
            [sp.conjugate(h13), sp.conjugate(h23), eps2],
        ]
    )


D6_BASE = {
    "E1": 1.787,
    "A1": 2.6,
    "B1": -3.8,
    "E2": -2.12,
    "A2": 1.63,
    "B2": 5.1,
    "L": 1.3,
    "C": 3.55,
}
D6_PANELS = {
    "a": {"M": 0.65, "F": 1.83, "D": 5.1},
    "b": {"M": 0.0, "F": 1.18, "D": 2.53},
    "c": {"M": 0.65, "F": 7.6, "D": 5.1},
}


def H_Ti3Al_sympy(params: dict[str, float] | None = None) -> sp.Matrix:
    p = {
        "A1": -9.66,
        "A2": 11.37,
        "B1": 36.22,
        "B2": -25.71,
        "M1": 0.12,
        "M2": -0.52,
        "C": 22.34,
    }
    if params:
        p.update(params)

    k2xy = kx**2 + ky**2
    h1 = p["A1"] * k2xy + p["B1"] * kz**2 + p["M1"]
    h2 = p["A2"] * k2xy + p["B2"] * kz**2 + p["M2"]
    h = 2 * p["C"] * kz
    eps = sp.Rational(1, 2) * (h1 + h2 + sp.sqrt((h1 - h2) ** 2 + h**2))
    return sp.Matrix([[eps, 0], [0, -eps]])


def H_YH3_sympy(params: dict[str, float] | None = None) -> sp.Matrix:
    p = {
        "m1": 2.99,
        "a1": 2.0,
        "r1": 1.032,
        "s1": 1.032,
        "t1": 1.032,
        "n1": 3,
        "m2": 2.96,
        "m3": 2.96,
        "a2": 4.0,
        "a3": 4.0,
    }
    if params:
        p.update(params)

    a = 3.14
    g1, g2, g3 = sp.sin(a * kz), sp.sin(a * kx), sp.sin(a * ky)
    h1 = p["a1"] * (
        p["r1"] * sp.cos(a * kx) ** p["n1"]
        + p["s1"] * sp.cos(a * ky) ** p["n1"]
        + p["t1"] * sp.cos(a * kz) ** p["n1"]
        - p["m1"]
    )
    h2 = p["a2"] * (sp.cos(a * kx) + sp.cos(a * ky) + sp.cos(a * kz) - p["m2"])
    h3 = p["a3"] * (sp.cos(a * kx) + sp.cos(a * ky) + sp.cos(a * kz) - p["m3"])
    eps = sp.sqrt((g1**2 + h1**2) * (g2**2 + h2**2) * (g3**2 + h3**2))
    return sp.Matrix([[eps, 0], [0, -eps]])


def H_Co2MnGa_TB6_sympy(params: dict[str, float] | None = None) -> sp.Matrix:
    p = {
        "t1": -0.31,
        "t2": -0.018,
        "t3": -0.01,
        "t4": 0.2,
        "t5": -0.02,
        "t6": 0.04,
        "t7": 0.28,
        "t8": -0.34,
        "eps_d": -0.6,
        "eps_p": 0.6,
    }
    if params:
        p.update(params)

    cx2, cy2, cz2 = sp.cos(kx / 2), sp.cos(ky / 2), sp.cos(kz / 2)
    sx2, sy2, sz2 = sp.sin(kx / 2), sp.sin(ky / 2), sp.sin(kz / 2)

    xi_d1 = (
        4 * p["t1"] * cx2 * cz2
        + 2 * p["t2"] * (sp.cos(kx) + sp.cos(kz))
        + 2 * p["t3"] * sp.cos(ky)
        + p["eps_d"]
    )
    xi_d2 = (
        4 * p["t1"] * cy2 * cz2
        + 2 * p["t2"] * (sp.cos(ky) + sp.cos(kz))
        + 2 * p["t3"] * sp.cos(kx)
        + p["eps_d"]
    )
    xi_d3 = (
        4 * p["t1"] * cx2 * cy2
        + 2 * p["t2"] * (sp.cos(kx) + sp.cos(ky))
        + 2 * p["t3"] * sp.cos(kz)
        + p["eps_d"]
    )

    xi_p1 = (
        4 * p["t4"] * cy2 * cz2
        + 2 * p["t5"] * (sp.cos(ky) + sp.cos(kz))
        + 2 * p["t6"] * sp.cos(kx)
        + p["eps_p"]
    )
    xi_p2 = (
        4 * p["t4"] * cx2 * cz2
        + 2 * p["t5"] * (sp.cos(kx) + sp.cos(kz))
        + 2 * p["t6"] * sp.cos(ky)
        + p["eps_p"]
    )
    xi_p3 = (
        4 * p["t4"] * cx2 * cy2
        + 2 * p["t5"] * (sp.cos(kx) + sp.cos(ky))
        + 2 * p["t6"] * sp.cos(kz)
        + p["eps_p"]
    )

    xi_p12 = -4 * p["t7"] * sx2 * sy2
    xi_p23 = -4 * p["t7"] * sy2 * sz2
    xi_p31 = -4 * p["t7"] * sx2 * sz2

    xi_dp11 = 2 * p["t8"] * sz2
    xi_dp22 = 2 * p["t8"] * sz2
    xi_dp13 = 2 * p["t8"] * sx2
    xi_dp32 = 2 * p["t8"] * sx2
    xi_dp23 = 2 * p["t8"] * sy2
    xi_dp31 = 2 * p["t8"] * sy2

    return sp.Matrix(
        [
            [xi_d1, 0, 0, xi_dp11, 0, xi_dp13],
            [0, xi_d2, 0, 0, xi_dp22, xi_dp23],
            [0, 0, xi_d3, xi_dp31, xi_dp32, 0],
            [xi_dp11, 0, xi_dp31, xi_p1, xi_p12, xi_p31],
            [0, xi_dp22, xi_dp32, xi_p12, xi_p2, xi_p23],
            [xi_dp13, xi_dp23, 0, xi_p31, xi_p23, xi_p3],
        ]
    )


@dataclass(frozen=True)
class SurfaceSpec:
    key: str
    title: str
    gap_tol: float
    span: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    dimension: int
    band_pair: tuple[int, int] = (0, 1)
    color: str = "#168fd3"
    opacity: float = 0.82
    old_notebook_cell: int | None = None


@dataclass(frozen=True)
class FigureSpec:
    key: str
    title: str
    hamiltonian: sp.Matrix
    params: dict[str, Any]
    surfaces: tuple[SurfaceSpec, ...]
    legacy_png: str


def cube_span(
    half_width: float,
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    value = float(half_width)
    return ((-value, value), (-value, value), (-value, value))


def positive_cube_span(
    width: float,
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    value = float(width)
    return ((0.0, value), (0.0, value), (0.0, value))


def legacy_specs() -> tuple[FigureSpec, ...]:
    panel_a_params = {**D6_BASE, **D6_PANELS["a"]}
    return (
        FigureSpec(
            key="tib2_d6_panel_a",
            title="TiB2 panel a",
            hamiltonian=H_D6_sympy(panel_a_params),
            params=panel_a_params,
            legacy_png="material_tib2_transition.png",
            surfaces=(
                SurfaceSpec(
                    "tib2_e_lt_2p8",
                    "E < 2.8 eV",
                    0.35,
                    cube_span(1.5),
                    200,
                    old_notebook_cell=9,
                ),
                SurfaceSpec(
                    "tib2_e_approx_2p8",
                    "E ~= 2.8 eV",
                    2.8,
                    cube_span(1.5),
                    200,
                    opacity=0.60,
                    old_notebook_cell=11,
                ),
                SurfaceSpec(
                    "tib2_mid",
                    "2.8 eV < E < 3.9 eV",
                    3.2,
                    cube_span(1.5),
                    200,
                    opacity=0.60,
                    old_notebook_cell=13,
                ),
                SurfaceSpec(
                    "tib2_high",
                    "3.9 eV < E",
                    4.0,
                    cube_span(1.5),
                    200,
                    opacity=0.60,
                    old_notebook_cell=16,
                ),
            ),
        ),
        FigureSpec(
            key="ti3al",
            title="Ti3Al",
            hamiltonian=H_Ti3Al_sympy(),
            params={
                "A1": -9.66,
                "A2": 11.37,
                "B1": 36.22,
                "B2": -25.71,
                "M1": 0.12,
                "M2": -0.52,
                "C": 22.34,
            },
            legacy_png="material_ti3al_transition.png",
            surfaces=(
                SurfaceSpec(
                    "ti3al_low",
                    "E < 0.25 eV",
                    0.3,
                    cube_span(math.pi),
                    250,
                    old_notebook_cell=45,
                ),
                SurfaceSpec(
                    "ti3al_high",
                    "0.25 eV < E",
                    0.6,
                    cube_span(math.pi),
                    250,
                    old_notebook_cell=47,
                ),
            ),
        ),
        FigureSpec(
            key="yh3",
            title="YH3",
            hamiltonian=H_YH3_sympy(),
            params={
                "m1": 2.99,
                "a1": 2.0,
                "r1": 1.032,
                "s1": 1.032,
                "t1": 1.032,
                "n1": 3,
                "m2": 2.96,
                "m3": 2.96,
                "a2": 4.0,
                "a3": 4.0,
                "a_lattice": 3.14,
            },
            legacy_png="material_yh3_transition.png",
            surfaces=(
                SurfaceSpec(
                    "yh3_low",
                    "E < 0.006 eV",
                    0.002,
                    cube_span(1.0),
                    200,
                    old_notebook_cell=52,
                ),
                SurfaceSpec(
                    "yh3_mid",
                    "0.006 eV < E < 0.2 eV",
                    0.007,
                    cube_span(1.0),
                    200,
                    opacity=0.60,
                    old_notebook_cell=54,
                ),
                SurfaceSpec(
                    "yh3_high",
                    "0.2 eV < E",
                    0.025,
                    cube_span(1.0),
                    200,
                    opacity=0.60,
                    old_notebook_cell=56,
                ),
            ),
        ),
        FigureSpec(
            key="co2mnga",
            title="Co2MnGa panel c",
            hamiltonian=H_Co2MnGa_TB6_sympy(),
            params={
                "t1": -0.31,
                "t2": -0.018,
                "t3": -0.01,
                "t4": 0.2,
                "t5": -0.02,
                "t6": 0.04,
                "t7": 0.28,
                "t8": -0.34,
                "eps_d": -0.6,
                "eps_p": 0.6,
            },
            legacy_png="material_co2mnga_transition.png",
            surfaces=(
                SurfaceSpec(
                    "co2mnga_low",
                    "E < 0.25 eV",
                    0.05,
                    cube_span(2.0 * math.pi),
                    200,
                    old_notebook_cell=63,
                ),
                SurfaceSpec(
                    "co2mnga_mid",
                    "0.25 eV < E < 1.15 eV",
                    0.5,
                    cube_span(2.0 * math.pi),
                    200,
                    old_notebook_cell=65,
                ),
                SurfaceSpec(
                    "co2mnga_high",
                    "1.15 eV < E",
                    1.2,
                    cube_span(2.05 * math.pi),
                    200,
                    opacity=0.60,
                    old_notebook_cell=67,
                ),
            ),
        ),
    )


def make_material_surface(
    hamiltonian: sp.Matrix,
    surface: SurfaceSpec,
    *,
    chunk_size: int,
) -> MaterialFermiSurface:
    return MaterialFermiSurface(
        hamiltonian,
        k_symbols=K_SYMBOLS,
        span=surface.span,
        dimension=surface.dimension,
        band_pair=surface.band_pair,
        gap_tol=surface.gap_tol,
        check_hermitian=False,
        check_pt_symmetry=False,
        chunk_size=chunk_size,
        force_small_edge_contraction=True,
        small_edge_limit=math.pi * 0.1,
    )


def contour_from_cached_gap(base: MaterialFermiSurface, gap_tol: float) -> pv.PolyData:
    gap = np.asarray(base.band_gap, dtype=np.float64)
    vol = pv.ImageData(
        dimensions=gap.shape,
        spacing=base.spacing * base.axis_scale,
        origin=base.origin,
    )
    vol.point_data["ES_helper"] = (gap - float(gap_tol)).ravel(order="F")
    return vol.contour(isosurfaces=[0.0], scalars="ES_helper")


def triangulate_for_plot(mesh: pv.PolyData, max_triangles: int) -> pv.PolyData:
    if mesh.n_points == 0 or mesh.n_cells == 0:
        return mesh
    tri = mesh.extract_geometry().triangulate().clean()
    if max_triangles > 0 and tri.n_cells > max_triangles:
        reduction = 1.0 - max_triangles / float(tri.n_cells)
        try:
            tri = tri.decimate_pro(reduction, preserve_topology=True).clean()
        except TypeError:
            tri = tri.decimate_pro(reduction).clean()
    return tri


def mesh_bounds(
    mesh: pv.PolyData, fallback_span: tuple[tuple[float, float], ...]
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    if mesh.n_points == 0:
        return fallback_span
    xmin, xmax, ymin, ymax, zmin, zmax = [float(v) for v in mesh.bounds]
    ranges = [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
    padded = []
    for (lo, hi), (flo, fhi) in zip(ranges, fallback_span):
        width = hi - lo
        if width <= 0:
            center = 0.5 * (lo + hi)
            width = max(abs(fhi - flo) * 0.05, 1e-3)
            lo, hi = center - 0.5 * width, center + 0.5 * width
        pad = max(0.06 * width, 1e-3)
        padded.append((lo - pad, hi + pad))
    return tuple(padded)  # type: ignore[return-value]


def camera_from_old_view(
    elevation_deg: float, azimuth_deg: float, distance: float = 1.75
) -> dict[str, dict[str, float]]:
    elevation = math.radians(elevation_deg)
    azimuth = math.radians(azimuth_deg)
    return {
        "eye": {
            "x": distance * math.cos(elevation) * math.cos(azimuth),
            "y": distance * math.cos(elevation) * math.sin(azimuth),
            "z": distance * math.sin(elevation),
        }
    }


def mesh_trace(
    mesh: pv.PolyData, surface: SurfaceSpec, fig_spec: FigureSpec
) -> go.BaseTraceType:
    tri = mesh.extract_geometry().triangulate()
    if tri.n_points == 0 or tri.n_cells == 0:
        return go.Scatter3d(
            x=[0],
            y=[0],
            z=[0],
            mode="text",
            text=[f"{surface.title}<br>empty contour"],
            showlegend=False,
        )

    faces = tri.faces.reshape((-1, 4))
    triangles = faces[faces[:, 0] == 3][:, 1:4]
    points = np.asarray(tri.points, dtype=float)

    hover = (
        f"{fig_spec.title}<br>"
        f"{surface.title}<br>"
        f"gap_tol={surface.gap_tol:g} eV<br>"
        f"N={surface.dimension}<br>"
        f"band_pair={surface.band_pair}<br>"
        f"old_cell={surface.old_notebook_cell}"
    )
    return go.Mesh3d(
        x=points[:, 0],
        y=points[:, 1],
        z=points[:, 2],
        i=triangles[:, 0],
        j=triangles[:, 1],
        k=triangles[:, 2],
        color=surface.color,
        opacity=surface.opacity,
        flatshading=False,
        lighting={
            "ambient": 0.35,
            "diffuse": 0.78,
            "specular": 0.12,
            "roughness": 0.48,
        },
        lightposition={"x": 120, "y": 180, "z": 240},
        hovertemplate=hover + "<extra></extra>",
        showscale=False,
        name=f"{fig_spec.key}:{surface.key}",
    )


def compute_surfaces(
    specs: tuple[FigureSpec, ...],
    *,
    only: set[str] | None,
    chunk_size: int,
    max_plot_triangles: int,
) -> tuple[list[dict[str, Any]], dict[str, pv.PolyData]]:
    records: list[dict[str, Any]] = []
    plot_meshes: dict[str, pv.PolyData] = {}

    for fig_spec in specs:
        if only and fig_spec.key not in only:
            continue
        print(f"\n[{fig_spec.key}] {fig_spec.title}", flush=True)
        groups: dict[tuple[Any, ...], list[SurfaceSpec]] = {}
        for surface in fig_spec.surfaces:
            group_key = (
                tuple(tuple(float(v) for v in pair) for pair in surface.span),
                int(surface.dimension),
                tuple(surface.band_pair),
            )
            groups.setdefault(group_key, []).append(surface)

        for _, group_surfaces in groups.items():
            first = group_surfaces[0]
            group_started = time.perf_counter()
            base = make_material_surface(
                fig_spec.hamiltonian, first, chunk_size=chunk_size
            )
            gap = base.band_gap
            group_elapsed = time.perf_counter() - group_started
            print(
                f"  sampled N={first.dimension}, span={first.span}, band={first.band_pair} "
                f"in {group_elapsed:.1f}s; gap range [{float(np.nanmin(gap)):.5g}, {float(np.nanmax(gap)):.5g}]",
                flush=True,
            )

            for surface in group_surfaces:
                contour_started = time.perf_counter()
                mesh = contour_from_cached_gap(base, surface.gap_tol)
                contour_elapsed = time.perf_counter() - contour_started
                plot_mesh = triangulate_for_plot(mesh, max_plot_triangles)
                plot_key = f"{fig_spec.key}:{surface.key}"
                plot_meshes[plot_key] = plot_mesh
                bounds = mesh_bounds(mesh, surface.span)
                record = {
                    "figure": fig_spec.key,
                    "figure_title": fig_spec.title,
                    "surface_key": surface.key,
                    "surface_title": surface.title,
                    "gap_tol_eV": surface.gap_tol,
                    "span": surface.span,
                    "dimension_N": surface.dimension,
                    "grid_points": int(surface.dimension**3),
                    "band_pair": surface.band_pair,
                    "old_notebook_cell": surface.old_notebook_cell,
                    "legacy_png": fig_spec.legacy_png,
                    "gap_min": float(np.nanmin(gap)),
                    "gap_max": float(np.nanmax(gap)),
                    "full_mesh_points": int(mesh.n_points),
                    "full_mesh_cells": int(mesh.n_cells),
                    "plot_mesh_points": int(plot_mesh.n_points),
                    "plot_mesh_cells": int(plot_mesh.n_cells),
                    "mesh_bounds": bounds,
                    "sample_elapsed_s": round(group_elapsed, 3),
                    "contour_elapsed_s": round(contour_elapsed, 3),
                    "material_api": "knotted_graph.applications.material_surface.MaterialFermiSurface",
                    "knotted_graph_root": package_location(),
                }
                records.append(record)
                print(
                    f"    {surface.title}: gap_tol={surface.gap_tol:g}, "
                    f"mesh={mesh.n_points} pts/{mesh.n_cells} cells, "
                    f"plot={plot_mesh.n_points} pts/{plot_mesh.n_cells} cells, "
                    f"contour {contour_elapsed:.1f}s",
                    flush=True,
                )

    return records, plot_meshes


def build_plotly_gallery(
    specs: tuple[FigureSpec, ...],
    records: list[dict[str, Any]],
    plot_meshes: dict[str, pv.PolyData],
    output_html: Path,
) -> None:
    rows = [
        spec
        for spec in specs
        if any(f"{spec.key}:{s.key}" in plot_meshes for s in spec.surfaces)
    ]
    cols = max(len(spec.surfaces) for spec in rows)
    titles = []
    for spec in rows:
        for surface in spec.surfaces:
            if f"{spec.key}:{surface.key}" in plot_meshes:
                titles.append(
                    f"{spec.title}<br>{surface.title}<br>cut={surface.gap_tol:g} eV"
                )
            else:
                titles.append("")
        titles.extend([""] * (cols - len(spec.surfaces)))

    fig = make_subplots(
        rows=len(rows),
        cols=cols,
        specs=[[{"type": "scene"} for _ in range(cols)] for _ in rows],
        horizontal_spacing=0.01,
        vertical_spacing=0.045,
        subplot_titles=titles,
    )

    record_by_key = {f"{r['figure']}:{r['surface_key']}": r for r in records}
    scene_layouts: dict[str, Any] = {}
    default_camera = camera_from_old_view(-10, 100)
    yh3_camera = camera_from_old_view(-25, 20)
    co_camera = camera_from_old_view(-5, 100)

    for row, fig_spec in enumerate(rows, start=1):
        for col, surface in enumerate(fig_spec.surfaces, start=1):
            key = f"{fig_spec.key}:{surface.key}"
            mesh = plot_meshes.get(key)
            if mesh is None:
                continue
            fig.add_trace(mesh_trace(mesh, surface, fig_spec), row=row, col=col)

            scene_index = (row - 1) * cols + col
            scene_name = "scene" if scene_index == 1 else f"scene{scene_index}"
            bounds = record_by_key[key]["mesh_bounds"]
            camera = default_camera
            if fig_spec.key == "yh3":
                camera = yh3_camera
            elif fig_spec.key == "co2mnga":
                camera = co_camera

            scene_layouts[scene_name] = {
                "xaxis": {
                    "title": "kx (1/A)",
                    "range": list(bounds[0]),
                    "showbackground": False,
                    "gridcolor": "#d9d9d9",
                    "zeroline": False,
                },
                "yaxis": {
                    "title": "ky (1/A)",
                    "range": list(bounds[1]),
                    "showbackground": False,
                    "gridcolor": "#d9d9d9",
                    "zeroline": False,
                },
                "zaxis": {
                    "title": "kz (1/A)",
                    "range": list(bounds[2]),
                    "showbackground": False,
                    "gridcolor": "#d9d9d9",
                    "zeroline": False,
                },
                "aspectmode": "cube",
                "camera": camera,
            }

    fig.update_layout(
        **scene_layouts,
        title={
            "text": "Reference material Hamiltonian surfaces reproduced with KnottedGraph",
            "x": 0.5,
            "xanchor": "center",
        },
        margin={"l": 12, "r": 12, "t": 115, "b": 24},
        height=max(520, 410 * len(rows)),
        width=1660,
        showlegend=False,
        paper_bgcolor="white",
        plot_bgcolor="white",
        font={"family": "Times New Roman, Times, serif", "size": 14, "color": "black"},
    )

    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(output_html),
        include_plotlyjs=True,
        full_html=True,
        config={"displayModeBar": True, "responsive": True},
    )


def write_summary(records: list[dict[str, Any]], output_json: Path) -> None:
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(records, indent=2), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("_build/new_phase_maps/material_surfaces"),
        help="Directory for the HTML gallery and JSON provenance.",
    )
    parser.add_argument(
        "--only",
        nargs="*",
        choices=[spec.key for spec in legacy_specs()],
        default=None,
        help="Optional subset of figure keys to compute.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100_000,
        help="Point chunk size for multiband eigensolves.",
    )
    parser.add_argument(
        "--max-plot-triangles",
        type=int,
        default=90_000,
        help="Decimate each displayed Plotly mesh to at most this many triangles. Set <=0 to disable.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    specs = legacy_specs()
    started = time.perf_counter()
    records, plot_meshes = compute_surfaces(
        specs,
        only=set(args.only) if args.only else None,
        chunk_size=args.chunk_size,
        max_plot_triangles=args.max_plot_triangles,
    )
    html_path = args.output_dir / "multiband_material_surface_check.html"
    json_path = args.output_dir / "multiband_material_surface_check_summary.json"
    build_plotly_gallery(specs, records, plot_meshes, html_path)
    write_summary(records, json_path)
    elapsed = time.perf_counter() - started
    print(f"\nWrote {html_path}")
    print(f"Wrote {json_path}")
    print(f"Total elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
