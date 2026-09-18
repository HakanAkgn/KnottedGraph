#!/usr/bin/env python3
"""Compact finite TPMS scaffold scans through the KnottedGraph pipeline.

This is a non-Hamiltonian companion to ``material_parameter_phase_maps.py``.
It uses the installed library's skeletonization, spatial-graph, and Yamada
tools. Use the package CLI for a small, explicitly bounded introductory scan.

Workflow per grid point:

    (lambda, c)
    -> sample F_lambda(x,y,z)
    -> compact scaffold mask Omega = {F_lambda <= c} cap {Phi_domain <= 0}
    -> Lee volume skeletonization
    -> sparse embedded spatial graph
    -> exact Yamada polynomial when tractable, otherwise graph/audit signature
    -> phase map of Yamada-distinguished signatures

The explicit finite design domain is the literature-backed correction compared
with a clipped periodic cell: every reconstructed handlebody is required to sit
strictly inside the sampling box.

The scalar-field side is intentionally generic: add another ``ImplicitField``
and another ``TPMSFamily`` to pass a new field through the same interface.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import gzip
import json
import math
import time
from collections import Counter, deque
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Callable, Iterable

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pyvista as pv
import sympy as sp
from matplotlib.colors import BoundaryNorm, ListedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import euler_number, label as label_volume, marching_cubes


from ._runtime import package_location

from knotted_graph.applications.phase_maps import (  # noqa: E402
    _closed_genus_zero_interior,
    _compute_yamada_audited,
    _graph_summary,
    _one_vertex_graph,
    volume_topology,
)
from knotted_graph.core import (  # noqa: E402
    is_trivalent,
    remove_leaf_nodes,
    simplify_edges,
    smooth_edges,
)
from knotted_graph.extraction import (  # noqa: E402
    skeleton_image_to_graph,
    skeletonize_volume,
)


Y = sp.Symbol("Y")
TEXT_COLOR = "#111827"
PHASE_COLORS = (
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

ScalarField = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]
Span3D = tuple[tuple[float, float], tuple[float, float], tuple[float, float]]


@dataclass(frozen=True)
class ImplicitField:
    """Analytic scalar field F(x,y,z) for a TPMS level-set family."""

    key: str
    title: str
    formula: str
    function: ScalarField
    reference_note: str

    def sample(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> np.ndarray:
        values = np.asarray(self.function(x, y, z), dtype=np.float64)
        if values.shape != x.shape:
            values = np.broadcast_to(values, x.shape).astype(np.float64, copy=True)
        if not np.isfinite(values).all():
            raise ValueError(f"{self.key} produced NaN or infinite samples")
        return values


@dataclass(frozen=True)
class SolidDomain:
    """Finite design domain Phi(x,y,z) <= 0 used to compactify TPMS scaffolds."""

    key: str
    title: str
    formula: str
    function: ScalarField
    reference_note: str

    def sample(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> np.ndarray:
        values = np.asarray(self.function(x, y, z), dtype=np.float64)
        if values.shape != x.shape:
            values = np.broadcast_to(values, x.shape).astype(np.float64, copy=True)
        if not np.isfinite(values).all():
            raise ValueError(f"{self.key} produced NaN or infinite samples")
        return values


@dataclass(frozen=True)
class TPMSFamily:
    """One interpolation family F_lambda = (1-lambda)F0 + lambda F1."""

    key: str
    title: str
    start: ImplicitField
    end: ImplicitField
    thresholds: tuple[float, ...]
    span: Span3D
    dimension: int
    domain: SolidDomain
    convention: str = "Omega(lambda,c) = {r : F_lambda(r) <= c and Phi_domain(r) <= 0}"
    rationale: str = ""

    def field_at(
        self,
        lam: float,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> np.ndarray:
        lam = float(lam)
        if not 0.0 <= lam <= 1.0:
            raise ValueError("lambda must lie in [0, 1]")
        left = self.start.sample(x, y, z)
        right = self.end.sample(x, y, z)
        return (1.0 - lam) * left + lam * right

    def compact_level_field(
        self,
        raw_values: np.ndarray,
        domain_values: np.ndarray,
        threshold_c: float,
    ) -> np.ndarray:
        """Return G=max(F_lambda-c, Phi_domain), whose G<=0 set is compact."""
        return np.maximum(
            np.asarray(raw_values, dtype=np.float64) - float(threshold_c),
            np.asarray(domain_values, dtype=np.float64),
        )


@dataclass
class PhaseCell:
    family: str
    title: str
    start_field: str
    end_field: str
    lam: float
    threshold_c: float
    dimension: int
    span: Span3D
    convention: str
    phase_signature: str
    phase_label: str
    source: str
    polynomial: str
    nodes: int
    edges: int
    components: int
    cycle_rank: int
    degree_sequence: tuple[int, ...]
    total_edge_points: int
    is_trivalent: bool
    interior_voxels: int
    interior_fraction: float
    interior_components: int
    euler_characteristic: int
    handle_rank: int
    touches_boundary: bool
    boundary_faces: tuple[str, ...]
    skeleton_voxels: int
    surface_points: int
    surface_cells: int
    surface_open_edges: int
    surface_nonmanifold_edges: int
    surface_is_closed: bool
    field_min: float
    field_max: float
    field_mean: float
    field_std: float
    compact_level_min: float
    compact_level_max: float
    compact_level_mean: float
    compact_level_std: float
    exact_yamada_attempted: bool
    geometry_path: str
    sample_seconds: float
    skeleton_seconds: float
    graph_seconds: float
    yamada_seconds: float
    surface_seconds: float
    total_seconds: float
    error: str | None = None
    evaluation_kind: str = "unrecorded"
    normalization: str | None = None
    projection: dict[str, Any] | None = None
    is_subcubic: bool | None = None
    classification_computed: bool = True
    component_count_matches: bool | None = None
    void_components: int = 0


class ImplicitSolidRegion:
    """Adapter exposing the same region-to-graph surface used by the pipeline."""

    def __init__(
        self,
        values: np.ndarray,
        *,
        threshold_c: float,
        span: Span3D,
        axis_scale: tuple[float, float, float] = (1.0, 1.0, 1.0),
    ) -> None:
        arr = np.asarray(values, dtype=np.float64)
        if arr.ndim != 3:
            raise ValueError("values must be a three-dimensional scalar field")
        if len(set(arr.shape)) != 1:
            raise ValueError("this adapter expects a cubic scalar grid")
        if not np.isfinite(arr).all():
            raise ValueError("values contain NaN or infinite samples")

        self.values = arr
        self.threshold_c = float(threshold_c)
        self.span = np.asarray(span, dtype=float)
        self.dimension = int(arr.shape[0])
        self.spacing = np.diff(self.span, axis=1).squeeze() / (self.dimension - 1)
        self.axis_scale = np.asarray(axis_scale, dtype=float)
        self.origin = self.span[:, 0]
        self.skeleton_graph_cache = None
        self.skeleton_graph_cache_args = None

    @property
    def _interior_mask(self) -> np.ndarray:
        return self.values <= self.threshold_c

    @cached_property
    def _skeleton_image(self) -> np.ndarray:
        return skeletonize_volume(self._interior_mask)

    @cached_property
    def fields_pv(self) -> pv.ImageData:
        vol = pv.ImageData(
            dimensions=self.values.shape,
            spacing=self.spacing * self.axis_scale,
            origin=self.origin,
        )
        vol.point_data["field"] = self.values.ravel(order="F")
        vol.point_data["level_helper"] = (self.values - self.threshold_c).ravel(
            order="F"
        )
        vol.point_data["solid"] = self._interior_mask.astype(np.uint8).ravel(order="F")
        return vol

    @cached_property
    def exceptional_surface_pv(self) -> pv.PolyData:
        return self.fields_pv.contour(isosurfaces=[0.0], scalars="level_helper")

    def skeleton_graph(
        self,
        simplify: bool = True,
        smooth_epsilon: int = 0,
        *,
        skeleton_image: np.ndarray | nx.Graph | nx.MultiGraph | None = None,
        force_small_edge_contraction: bool = False,
        small_edge_limit: float = 0.0,
        previous_n_edgepoint: int = 20,
    ) -> nx.MultiGraph:
        """Convert the solid-region skeleton to the standard embedded graph."""
        del previous_n_edgepoint
        args = (
            smooth_epsilon,
            simplify,
            id(skeleton_image),
            bool(force_small_edge_contraction),
            float(small_edge_limit),
        )
        if (
            self.skeleton_graph_cache is not None
            and self.skeleton_graph_cache_args == args
        ):
            return self.skeleton_graph_cache

        if skeleton_image is None:
            topology = volume_topology(self._interior_mask)
            if topology.enclosed_voids:
                raise ValueError("volume has enclosed voids; a single graph is not a spine")
            graph = skeleton_image_to_graph(
                self._skeleton_image,
                expected_cycle_rank=topology.handle_rank,
                expected_components=topology.connected_components,
            )
        elif isinstance(skeleton_image, (nx.Graph, nx.MultiGraph)):
            graph = (
                skeleton_image
                if isinstance(skeleton_image, nx.MultiGraph)
                else nx.MultiGraph(skeleton_image)
            )
        else:
            graph = skeleton_image_to_graph(np.asarray(skeleton_image, dtype=bool))

        if simplify:
            graph = remove_leaf_nodes(graph)
            if graph.number_of_edges():
                graph = simplify_edges(graph)

        if force_small_edge_contraction and small_edge_limit > 0:
            from knotted_graph.core import contract_short_edges

            graph = contract_short_edges(
                graph,
                min_length=float(small_edge_limit),
                copy=False,
            )

        if graph.number_of_edges():
            graph = smooth_edges(graph, epsilon=float(smooth_epsilon), copy=False)
        graph.graph["is_trivalent"] = is_trivalent(graph)
        self.skeleton_graph_cache = graph
        self.skeleton_graph_cache_args = args
        return graph

    def idx_to_world(self, points: Any) -> np.ndarray:
        arr = np.asarray(points, dtype=float)
        return arr * (self.spacing * self.axis_scale) + self.origin


def gyroid(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def schwarz_p(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    return np.cos(x) + np.cos(y) + np.cos(z)


def diamond(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    return np.cos(x) * np.cos(y) * np.cos(z) - np.sin(x) * np.sin(y) * np.sin(z)


def spherical_domain(radius: float) -> ScalarField:
    radius = float(radius)

    def field(
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> np.ndarray:
        return np.sqrt(x * x + y * y + z * z) / radius - 1.0

    return field


def cylindrical_domain(radius: float, half_height: float) -> ScalarField:
    radius = float(radius)
    half_height = float(half_height)

    def field(
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> np.ndarray:
        radial = np.sqrt(x * x + y * y) / radius - 1.0
        axial = np.abs(z) / half_height - 1.0
        return np.maximum(radial, axial)

    return field


def compact_domain(
    *,
    kind: str,
    span: Span3D,
    radius_fraction: float,
) -> SolidDomain:
    half_widths = [min(abs(lo), abs(hi)) for lo, hi in span]
    min_half_width = float(min(half_widths))
    radius = float(radius_fraction) * min_half_width
    if not 0.0 < radius < min_half_width:
        raise ValueError("domain radius must be strictly inside the sampled span")
    if kind == "sphere":
        return SolidDomain(
            key="sphere",
            title="Finite spherical scaffold coupon",
            formula=f"sqrt(x^2+y^2+z^2)/{radius:.6g} - 1 <= 0",
            function=spherical_domain(radius),
            reference_note=(
                "Finite TPMS scaffold formed by intersecting the periodic "
                "level-set phase with an explicit outer design domain."
            ),
        )
    if kind == "cylinder":
        half_height = radius
        return SolidDomain(
            key="cylinder",
            title="Finite cylindrical scaffold coupon",
            formula=(
                f"max(sqrt(x^2+y^2)/{radius:.6g} - 1, |z|/{half_height:.6g} - 1) <= 0"
            ),
            function=cylindrical_domain(radius, half_height),
            reference_note=(
                "Cylindrical TPMS scaffold coupon, matching the common finite "
                "specimen geometry used in scaffold/mechanics studies."
            ),
        )
    raise ValueError(f"unsupported compact domain kind: {kind}")


def implicit_fields() -> dict[str, ImplicitField]:
    """Return the endpoint fields used in the TPMS scans."""
    return {
        "gyroid": ImplicitField(
            key="gyroid",
            title="Gyroid",
            formula=("sin(x) cos(y) + sin(y) cos(z) + sin(z) cos(x)"),
            function=gyroid,
            reference_note=(
                "Standard gyroid TPMS nodal field; no package-specific field "
                "was found in this checkout."
            ),
        ),
        "schwarz_p": ImplicitField(
            key="schwarz_p",
            title="Schwarz-P",
            formula="cos(x) + cos(y) + cos(z)",
            function=schwarz_p,
            reference_note="Standard Schwarz primitive TPMS nodal field.",
        ),
        "diamond": ImplicitField(
            key="diamond",
            title="Diamond",
            formula="cos(x) cos(y) cos(z) - sin(x) sin(y) sin(z)",
            function=diamond,
            reference_note=(
                "Standard Schwarz-D / diamond nodal approximation used in "
                "TPMS scaffold literature."
            ),
        ),
    }


def cube_span(half_width: float = math.pi) -> Span3D:
    return ((-float(half_width), float(half_width)),) * 3


def tpms_families(
    *,
    dimension: int,
    thresholds: tuple[float, ...],
    domain_kind: str = "sphere",
    span_half_width: float = 2.25 * math.pi,
    domain_radius_fraction: float = 0.72,
) -> tuple[TPMSFamily, ...]:
    fields = implicit_fields()
    span = cube_span(float(span_half_width))
    domain = compact_domain(
        kind=domain_kind,
        span=span,
        radius_fraction=float(domain_radius_fraction),
    )
    return (
        TPMSFamily(
            key="gyroid_to_diamond",
            title="Gyroid -> Diamond",
            start=fields["gyroid"],
            end=fields["diamond"],
            thresholds=thresholds,
            span=span,
            dimension=dimension,
            domain=domain,
            rationale=(
                "Interpolates a chiral gyroid scaffold into the standard "
                "Schwarz-D/diamond TPMS scaffold inside a finite design domain."
            ),
        ),
        TPMSFamily(
            key="gyroid_to_schwarz_p",
            title="Gyroid -> Schwarz-P",
            start=fields["gyroid"],
            end=fields["schwarz_p"],
            thresholds=thresholds,
            span=span,
            dimension=dimension,
            domain=domain,
            rationale=(
                "Interpolates the gyroid field into the primitive Schwarz-P "
                "porous-material scaffold inside a finite design domain."
            ),
        ),
        TPMSFamily(
            key="schwarz_p_to_diamond",
            title="Schwarz-P -> Diamond",
            start=fields["schwarz_p"],
            end=fields["diamond"],
            thresholds=thresholds,
            span=span,
            dimension=dimension,
            domain=domain,
            rationale=(
                "Interpolates between two standard cubic TPMS scaffold fields "
                "inside a finite design domain."
            ),
        ),
    )


def sample_grid(
    span: Span3D,
    dimension: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    axes = [
        np.linspace(float(lo), float(hi), int(dimension), dtype=np.float64)
        for lo, hi in span
    ]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    return axes[0], x, y, z


def boundary_faces(mask: np.ndarray) -> tuple[str, ...]:
    mask = np.asarray(mask, dtype=bool)
    faces: list[str] = []
    for label, touched in (
        ("x_min", mask[0, :, :].any()),
        ("x_max", mask[-1, :, :].any()),
        ("y_min", mask[:, 0, :].any()),
        ("y_max", mask[:, -1, :].any()),
        ("z_min", mask[:, :, 0].any()),
        ("z_max", mask[:, :, -1].any()),
    ):
        if bool(touched):
            faces.append(label)
    return tuple(faces)


def interior_summary(mask: np.ndarray) -> dict[str, Any]:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return {
            "void_components": 0,
            "interior_voxels": 0,
            "interior_fraction": 0.0,
            "interior_components": 0,
            "euler_characteristic": 0,
            "handle_rank": 0,
            "touches_boundary": False,
            "boundary_faces": tuple(),
        }
    _, component_count = label_volume(mask, connectivity=3, return_num=True)
    euler = int(euler_number(mask, connectivity=3))
    topology = volume_topology(mask)
    handle_rank = topology.handle_rank
    faces = boundary_faces(mask)
    return {
        "void_components": int(topology.enclosed_voids),
        "interior_voxels": int(mask.sum()),
        "interior_fraction": float(mask.mean()),
        "interior_components": int(component_count),
        "euler_characteristic": euler,
        "handle_rank": int(handle_rank),
        "touches_boundary": bool(faces),
        "boundary_faces": faces,
    }


def multigraph_wl_hash(graph: nx.MultiGraph) -> str:
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
    return nx.weisfeiler_lehman_graph_hash(
        simple,
        node_attr="loops",
        edge_attr="m",
        iterations=3,
    )


def graph_signature(
    summary: dict[str, Any],
    interior: dict[str, Any],
    graph: nx.MultiGraph,
) -> str:
    degree_hist = tuple(sorted(Counter(summary["degree_sequence"]).items()))
    return (
        "core:"
        f"ic={interior['interior_components']};"
        f"h={interior['handle_rank']};"
        f"b={int(interior['touches_boundary'])};"
        f"gcomp={summary['components']};"
        f"beta={summary['cycle_rank']};"
        f"deg={degree_hist};"
        f"wl={multigraph_wl_hash(graph)}"
    )


def short_signature_label(signature: str, source: str, polynomial: str) -> str:
    if source == "vertex":
        return f"isolated vertices; Yamada {polynomial}"
    if source == "yamada":
        return f"Yamada {polynomial}"
    if source == "diagram-yamada":
        return f"fixed-diagram Yamada {polynomial}"
    if source == "large-core":
        chunks = dict(
            part.split("=", 1)
            for part in signature.removeprefix("core:").split(";")
            if "=" in part
        )
        return (
            f"core h={chunks.get('h', '?')}, beta={chunks.get('beta', '?')}, "
            f"deg={chunks.get('deg', '?')}, WL={chunks.get('wl', '?')[:10]}"
        )
    if source == "error":
        return "extraction error"
    return source


def canonical_yamada_string(expr: sp.Expr) -> tuple[str, str]:
    canonical = sp.factor(sp.together(sp.expand(expr)))
    return str(canonical), "yamada:" + sp.srepr(canonical)


def surface_stats_and_sample(
    values: np.ndarray,
    threshold_c: float,
    *,
    spacing: np.ndarray,
    origin: np.ndarray,
    max_sample_points: int,
) -> tuple[dict[str, Any], list[list[float]]]:
    if not (float(np.nanmin(values)) <= threshold_c <= float(np.nanmax(values))):
        return {
            "surface_points": 0,
            "surface_cells": 0,
            "surface_open_edges": 0,
            "surface_nonmanifold_edges": 0,
            "surface_is_closed": False,
        }, []
    try:
        verts, faces, _, _ = marching_cubes(
            values,
            level=float(threshold_c),
            spacing=tuple(float(v) for v in spacing),
        )
    except Exception:
        return {
            "surface_points": 0,
            "surface_cells": 0,
            "surface_open_edges": 0,
            "surface_nonmanifold_edges": 0,
            "surface_is_closed": False,
        }, []
    verts = verts + origin[None, :]
    edge_counts: Counter[tuple[int, int]] = Counter()
    for face in faces:
        a, b, c = (int(v) for v in face)
        for u, v in ((a, b), (b, c), (c, a)):
            if v < u:
                u, v = v, u
            edge_counts[(u, v)] += 1
    open_edges = sum(1 for count in edge_counts.values() if count == 1)
    nonmanifold_edges = sum(1 for count in edge_counts.values() if count != 2)
    if len(verts) and max_sample_points > 0:
        count = min(int(max_sample_points), len(verts))
        indices = np.linspace(0, len(verts) - 1, count, dtype=int)
        sample = verts[indices].round(6).tolist()
    else:
        sample = []
    return {
        "surface_points": int(len(verts)),
        "surface_cells": int(len(faces)),
        "surface_open_edges": int(open_edges),
        "surface_nonmanifold_edges": int(nonmanifold_edges),
        "surface_is_closed": bool(edge_counts and nonmanifold_edges == 0),
    }, sample


def downsample_polyline(points: np.ndarray, max_points: int) -> list[list[float]]:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3 or len(arr) == 0:
        return []
    if len(arr) > max_points:
        indices = np.linspace(0, len(arr) - 1, int(max_points), dtype=int)
        arr = arr[indices]
    return arr.round(6).tolist()


def graph_geometry(
    region: ImplicitSolidRegion,
    graph: nx.MultiGraph,
    *,
    max_edge_points: int,
) -> dict[str, Any]:
    nodes = []
    for node, data in graph.nodes(data=True):
        pos = np.asarray(data.get("pos", (0.0, 0.0, 0.0)), dtype=float)
        nodes.append(
            {
                "id": str(node),
                "index_pos": pos.tolist(),
                "coord": region.idx_to_world(pos).round(6).tolist(),
            }
        )

    edges = []
    for u, v, key, data in graph.edges(keys=True, data=True):
        pts = np.asarray(data.get("pts", []), dtype=float)
        if pts.ndim != 2 or pts.shape[1] != 3 or len(pts) == 0:
            start = np.asarray(graph.nodes[u]["pos"], dtype=float)
            end = np.asarray(graph.nodes[v]["pos"], dtype=float)
            pts = np.vstack([start, end])
        edges.append(
            {
                "u": str(u),
                "v": str(v),
                "key": str(key),
                "points_index": pts.tolist(),
                "display_downsampled": bool(len(pts) > max_edge_points),
                "points_coord": downsample_polyline(
                    region.idx_to_world(pts),
                    max_edge_points,
                ),
                "original_point_count": int(len(pts)),
            }
        )

    return {
        "geometry_role": "full-index-polylines-with-display-world-coordinates",
        "nodes": nodes,
        "edges": edges,
        "node_count": int(graph.number_of_nodes()),
        "edge_count": int(graph.number_of_edges()),
    }


def write_geometry_record(
    path: Path,
    *,
    family: TPMSFamily,
    lam: float,
    threshold_c: float,
    raw_values: np.ndarray,
    compact_values: np.ndarray,
    region: ImplicitSolidRegion,
    graph: nx.MultiGraph,
    interior: dict[str, Any],
    surface_stats: dict[str, Any],
    surface_sample: list[list[float]],
    max_edge_points: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "family": family.key,
        "title": family.title,
        "start_field": family.start.key,
        "end_field": family.end.key,
        "lambda": float(lam),
        "threshold_c": float(threshold_c),
        "convention": family.convention,
        "compact_domain": {
            "key": family.domain.key,
            "title": family.domain.title,
            "formula": family.domain.formula,
            "reference_note": family.domain.reference_note,
        },
        "compact_level_convention": (
            "G_lambda,c(r)=max(F_lambda(r)-c, Phi_domain(r)); "
            "the finite scaffold is {G_lambda,c <= 0}."
        ),
        "dimension": int(family.dimension),
        "span": family.span,
        "raw_field_stats": {
            "min": float(np.min(raw_values)),
            "max": float(np.max(raw_values)),
            "mean": float(np.mean(raw_values)),
            "std": float(np.std(raw_values)),
        },
        "compact_level_stats": {
            "min": float(np.min(compact_values)),
            "max": float(np.max(compact_values)),
            "mean": float(np.mean(compact_values)),
            "std": float(np.std(compact_values)),
        },
        "interior": {
            **{k: v for k, v in interior.items() if k != "boundary_faces"},
            "boundary_faces": list(interior["boundary_faces"]),
        },
        "surface": surface_stats,
        "skeleton_voxels": int(np.sum(region._skeleton_image))
        if "_skeleton_image" in region.__dict__
        else None,
        "graph": graph_geometry(region, graph, max_edge_points=max_edge_points),
        "surface_point_sample": surface_sample,
        "note": (
            "This geometry is a finite compact per-cell audit sample, not a complete "
            "surface mesh. Full topology classification is recorded by the "
            "spatial graph and Yamada/audit signature."
        ),
    }
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def evaluate_cell(
    family: TPMSFamily,
    raw_values: np.ndarray,
    domain_values: np.ndarray,
    lam: float,
    threshold_c: float,
    *,
    output_dir: Path,
    lambda_index: int,
    threshold_index: int,
    sample_seconds: float,
    max_exact_yamada_edges: int,
    max_surface_sample_points: int,
    max_edge_points: int,
) -> PhaseCell:
    t_cell = time.perf_counter()
    compact_values = family.compact_level_field(
        raw_values,
        domain_values,
        float(threshold_c),
    )
    region = ImplicitSolidRegion(
        compact_values,
        threshold_c=0.0,
        span=family.span,
    )
    mask = np.asarray(region._interior_mask, dtype=bool)
    interior = interior_summary(mask)

    skeleton_seconds = 0.0
    graph_seconds = 0.0
    yamada_seconds = 0.0
    surface_seconds = 0.0
    surface_stats = {
        "surface_points": 0,
        "surface_cells": 0,
        "surface_open_edges": 0,
        "surface_nonmanifold_edges": 0,
        "surface_is_closed": False,
    }
    surface_sample: list[list[float]] = []
    graph = nx.MultiGraph()
    summary = _graph_summary(graph)
    polynomial = ""
    signature = ""
    source = ""
    attempted = False
    error = None

    audit: dict[str, Any] = {}
    skeleton_voxels = 0
    component_count_matches = None
    try:
        t0 = time.perf_counter()
        if _closed_genus_zero_interior(mask):
            graph = _one_vertex_graph()
            graph.nodes[0]["pos"] = np.argwhere(mask).mean(axis=0)
            skeleton_seconds = time.perf_counter() - t0
        else:
            _ = region._skeleton_image
            skeleton_voxels = int(np.sum(region._skeleton_image))
            skeleton_seconds = time.perf_counter() - t0
            t0 = time.perf_counter()
            graph = region.skeleton_graph(smooth_epsilon=0, simplify=True)
            graph_seconds = time.perf_counter() - t0
        summary = _graph_summary(graph)
        component_count_matches = summary["components"] == interior["interior_components"]
        if not component_count_matches:
            raise ValueError(
                "reconstruction component mismatch: "
                f"volume={interior['interior_components']}, graph={summary['components']}"
            )
        if graph.number_of_nodes() == 0:
            raise ValueError("reconstruction produced an empty graph")
        if interior["void_components"]:
            raise ValueError("volume has enclosed voids; a single graph is not a spine")
        if summary["cycle_rank"] != interior["handle_rank"]:
            raise ValueError(
                "reconstruction cycle-rank mismatch: "
                f"volume={interior['handle_rank']}, graph={summary['cycle_rank']}"
            )
        if graph.number_of_edges() <= max_exact_yamada_edges:
            attempted = True
            t0 = time.perf_counter()
            polynomial_expr, audit = _compute_yamada_audited(graph, Y, {"normalize": True})
            polynomial, signature = canonical_yamada_string(polynomial_expr)
            yamada_seconds = time.perf_counter() - t0
            if not audit["is_subcubic"]:
                signature = signature.replace("yamada:", "diagram-yamada:", 1)
                source = "diagram-yamada"
            else:
                source = "vertex" if graph.number_of_edges() == 0 else "yamada"
        else:
            signature = graph_signature(summary, interior, graph)
            source = "large-core"
    except Exception as exc:
        # Preserve the actual graph and failure, never invent a genus-zero
        # vertex or substitute a crossing-free polynomial after an error.
        summary = _graph_summary(graph)
        polynomial = ""
        signature = "error:" + type(exc).__name__ + ":" + str(exc)[:120]
        source = "error"
        error = type(exc).__name__ + ": " + str(exc)[:240]

    t0 = time.perf_counter()
    try:
        surface_stats, surface_sample = surface_stats_and_sample(
            compact_values, 0.0,
            spacing=region.spacing * region.axis_scale,
            origin=region.origin,
            max_sample_points=max_surface_sample_points,
        )
    except Exception as exc:
        error = error or "surface:" + type(exc).__name__ + ": " + str(exc)[:200]
    surface_seconds = time.perf_counter() - t0

    geometry_path = (
        output_dir
        / "geometry"
        / f"{family.key}_lambda{lambda_index:03d}_c{threshold_index:03d}.json.gz"
    )
    try:
        write_geometry_record(
            geometry_path,
            family=family,
            lam=float(lam),
            threshold_c=float(threshold_c),
            raw_values=raw_values,
            compact_values=compact_values,
            region=region,
            graph=graph,
            interior=interior,
            surface_stats=surface_stats,
            surface_sample=surface_sample,
            max_edge_points=max_edge_points,
        )
    except Exception as exc:
        if error is None:
            error = "geometry-write:" + type(exc).__name__ + ": " + str(exc)[:200]

    total_seconds = time.perf_counter() - t_cell
    label = short_signature_label(signature, source, polynomial)
    return PhaseCell(
        family=family.key,
        title=family.title,
        start_field=family.start.key,
        end_field=family.end.key,
        lam=float(lam),
        threshold_c=float(threshold_c),
        dimension=int(family.dimension),
        span=family.span,
        convention=family.convention,
        phase_signature=signature,
        phase_label=label,
        source=source,
        polynomial=polynomial,
        nodes=int(summary["nodes"]),
        edges=int(summary["edges"]),
        components=int(summary["components"]),
        cycle_rank=int(summary["cycle_rank"]),
        degree_sequence=tuple(int(value) for value in summary["degree_sequence"]),
        total_edge_points=int(summary["total_edge_points"]),
        is_trivalent=bool(graph.graph.get("is_trivalent", is_trivalent(graph))),
        skeleton_voxels=int(skeleton_voxels),
        surface_points=int(surface_stats["surface_points"]),
        surface_cells=int(surface_stats["surface_cells"]),
        surface_open_edges=int(surface_stats["surface_open_edges"]),
        surface_nonmanifold_edges=int(surface_stats["surface_nonmanifold_edges"]),
        surface_is_closed=bool(surface_stats["surface_is_closed"]),
        field_min=float(np.min(raw_values)),
        field_max=float(np.max(raw_values)),
        field_mean=float(np.mean(raw_values)),
        field_std=float(np.std(raw_values)),
        compact_level_min=float(np.min(compact_values)),
        compact_level_max=float(np.max(compact_values)),
        compact_level_mean=float(np.mean(compact_values)),
        compact_level_std=float(np.std(compact_values)),
        exact_yamada_attempted=bool(attempted),
        geometry_path=str(geometry_path),
        sample_seconds=float(sample_seconds),
        skeleton_seconds=float(skeleton_seconds),
        graph_seconds=float(graph_seconds),
        yamada_seconds=float(yamada_seconds),
        surface_seconds=float(surface_seconds),
        total_seconds=float(total_seconds),
        error=error,
        evaluation_kind=audit.get("evaluation_kind", "failed" if source == "error" else "not-evaluated"),
        normalization=audit.get("normalization"),
        projection=audit.get("projection"),
        is_subcubic=max(dict(graph.degree()).values(), default=0) <= 3,
        classification_computed=True,
        component_count_matches=component_count_matches,
        **interior,
    )


def compute_family(
    family: TPMSFamily,
    lambdas: np.ndarray,
    *,
    output_dir: Path,
    max_exact_yamada_edges: int,
    max_surface_sample_points: int,
    max_edge_points: int,
) -> list[PhaseCell]:
    records: list[PhaseCell] = []
    _, x, y, z = sample_grid(family.span, family.dimension)
    domain_values = family.domain.sample(x, y, z)
    for lambda_index, lam in enumerate(lambdas):
        t0 = time.perf_counter()
        raw_values = family.field_at(float(lam), x, y, z)
        sample_seconds = time.perf_counter() - t0
        print(
            f"  {family.key} lambda {lambda_index + 1:02d}/{len(lambdas)} "
            f"lambda={float(lam):.4f}: field "
            f"[{float(raw_values.min()):.5g}, {float(raw_values.max()):.5g}] "
            f"in {sample_seconds:.3f}s",
            flush=True,
        )
        for threshold_index, threshold_c in enumerate(family.thresholds):
            records.append(
                evaluate_cell(
                    family,
                    raw_values,
                    domain_values,
                    float(lam),
                    float(threshold_c),
                    output_dir=output_dir,
                    lambda_index=lambda_index,
                    threshold_index=threshold_index,
                    sample_seconds=sample_seconds,
                    max_exact_yamada_edges=max_exact_yamada_edges,
                    max_surface_sample_points=max_surface_sample_points,
                    max_edge_points=max_edge_points,
                )
            )
    return records


def connected_components_for_label(
    labels: np.ndarray,
    phase_id: int,
) -> list[list[tuple[int, int]]]:
    h, w = labels.shape
    seen = np.zeros(labels.shape, dtype=bool)
    components: list[list[tuple[int, int]]] = []
    for i in range(h):
        for j in range(w):
            if seen[i, j] or int(labels[i, j]) != int(phase_id):
                continue
            queue = deque([(i, j)])
            seen[i, j] = True
            component: list[tuple[int, int]] = []
            while queue:
                ci, cj = queue.popleft()
                component.append((ci, cj))
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ni, nj = ci + di, cj + dj
                    if (
                        0 <= ni < h
                        and 0 <= nj < w
                        and not seen[ni, nj]
                        and int(labels[ni, nj]) == int(phase_id)
                    ):
                        seen[ni, nj] = True
                        queue.append((ni, nj))
            components.append(component)
    return components


def stable_labels(
    labels: np.ndarray,
    *,
    min_cells: int,
    protected_labels: Iterable[int] = (),
) -> tuple[np.ndarray, int]:
    """Optional display-only filtering; never evidence of topological equivalence.

    In each of at most three synchronous passes, a four-connected component
    smaller than ``min_cells`` can join an adjacent component already at least
    that size. Boundary contacts supply votes; a tied maximum leaves the small
    component unchanged. Protected labels are neither changed nor recipients.
    This rule is invariant under renaming label IDs. Return distinct changed
    cells relative to the input, not the number of intermediate assignments.
    """
    stable = labels.copy()
    if min_cells <= 1:
        return stable, 0
    protected = {int(value) for value in protected_labels}
    for _ in range(3):
        large = np.zeros(stable.shape, dtype=bool)
        components = []
        for phase_id in sorted(set(int(v) for v in stable.ravel())):
            if phase_id in protected:
                continue
            for component in connected_components_for_label(stable, phase_id):
                if len(component) >= min_cells:
                    for i, j in component:
                        large[i, j] = True
                else:
                    components.append((phase_id, component))
        updated = stable.copy()
        for phase_id, component in components:
            neighbor_counts: Counter[int] = Counter()
            for i, j in component:
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ni, nj = i + di, j + dj
                    if 0 <= ni < stable.shape[0] and 0 <= nj < stable.shape[1]:
                        neighbor = int(stable[ni, nj])
                        if large[ni, nj] and neighbor != phase_id:
                            neighbor_counts[neighbor] += 1
            if not neighbor_counts:
                continue
            maximum = max(neighbor_counts.values())
            winners = [k for k, count in neighbor_counts.items() if count == maximum]
            if len(winners) != 1:
                continue
            for i, j in component:
                updated[i, j] = winners[0]
        if np.array_equal(updated, stable):
            break
        stable = updated
    return stable, int(np.count_nonzero(stable != labels))


def phase_grid(
    records: list[PhaseCell],
    lambdas: np.ndarray,
    thresholds: tuple[float, ...],
) -> tuple[np.ndarray, dict[int, str], dict[str, int]]:
    signatures: list[str] = []
    seen: set[str] = set()
    for record in records:
        if record.phase_signature not in seen:
            signatures.append(record.phase_signature)
            seen.add(record.phase_signature)
    signature_to_id = {signature: idx + 1 for idx, signature in enumerate(signatures)}
    label_lookup = {
        idx + 1: next(r.phase_label for r in records if r.phase_signature == signature)
        for idx, signature in enumerate(signatures)
    }
    lookup = {(round(r.threshold_c, 12), round(r.lam, 12)): r for r in records}
    grid = np.zeros((len(thresholds), len(lambdas)), dtype=int)
    for row, threshold_c in enumerate(thresholds):
        for col, lam in enumerate(lambdas):
            record = lookup[(round(float(threshold_c), 12), round(float(lam), 12))]
            grid[row, col] = signature_to_id[record.phase_signature]
    return grid, label_lookup, signature_to_id


def record_dict(record: PhaseCell) -> dict[str, Any]:
    row = dataclasses.asdict(record)
    row["degree_sequence"] = list(record.degree_sequence)
    row["boundary_faces"] = list(record.boundary_faces)
    return row


def timing_summary(records: list[PhaseCell], elapsed_seconds: float) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "elapsed_seconds": round(float(elapsed_seconds), 3),
        "reconstructed_volumes": int(len(records)),
    }
    for field in (
        "sample_seconds",
        "skeleton_seconds",
        "graph_seconds",
        "yamada_seconds",
        "surface_seconds",
        "total_seconds",
    ):
        values = np.asarray([getattr(r, field) for r in records], dtype=float)
        summary[field] = {
            "sum": round(float(values.sum()), 3),
            "min": round(float(values.min()), 6),
            "median": round(float(np.median(values)), 6),
            "max": round(float(values.max()), 6),
        }
    return summary


def field_correlation(a: np.ndarray, b: np.ndarray) -> float:
    x = np.asarray(a, dtype=float).ravel()
    y = np.asarray(b, dtype=float).ravel()
    x = x - float(x.mean())
    y = y - float(y.mean())
    denom = float(np.linalg.norm(x) * np.linalg.norm(y))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(x, y) / denom)


def verify_endpoints(
    families: tuple[TPMSFamily, ...],
    *,
    endpoint_threshold: float = 0.0,
) -> list[dict[str, Any]]:
    fields = implicit_fields()
    records: list[dict[str, Any]] = []
    for family in families:
        _, x, y, z = sample_grid(family.span, family.dimension)
        domain_values = family.domain.sample(x, y, z)
        expected_samples = {key: field.sample(x, y, z) for key, field in fields.items()}
        for lam, expected_key in ((0.0, family.start.key), (1.0, family.end.key)):
            raw_values = family.field_at(lam, x, y, z)
            compact_values = family.compact_level_field(
                raw_values,
                domain_values,
                float(endpoint_threshold),
            )
            expected = expected_samples[expected_key]
            region = ImplicitSolidRegion(
                compact_values,
                threshold_c=0.0,
                span=family.span,
            )
            mask = region._interior_mask
            try:
                skeleton = region._skeleton_image
                graph = region.skeleton_graph(smooth_epsilon=0, simplify=True)
                graph_stats = _graph_summary(graph)
                skeleton_voxels = int(np.sum(skeleton))
            except Exception as exc:
                graph_stats = {
                    "nodes": 0,
                    "edges": 0,
                    "components": 0,
                    "cycle_rank": 0,
                    "degree_sequence": tuple(),
                    "total_edge_points": 0,
                }
                skeleton_voxels = 0
                graph_stats["error"] = f"{type(exc).__name__}: {exc}"
            surface_stats, _ = surface_stats_and_sample(
                compact_values,
                0.0,
                spacing=region.spacing * region.axis_scale,
                origin=region.origin,
                max_sample_points=0,
            )
            records.append(
                {
                    "family": family.key,
                    "lambda": lam,
                    "expected_field": expected_key,
                    "self_max_abs_error": float(np.max(np.abs(raw_values - expected))),
                    "correlations": {
                        key: field_correlation(raw_values, sample)
                        for key, sample in expected_samples.items()
                    },
                    "threshold_c": float(endpoint_threshold),
                    "compact_domain": family.domain.key,
                    "volume_fraction": float(mask.mean()),
                    "interior": {
                        **{
                            k: v
                            for k, v in interior_summary(mask).items()
                            if k != "boundary_faces"
                        },
                        "boundary_faces": list(boundary_faces(mask)),
                    },
                    "skeleton_voxels": skeleton_voxels,
                    "graph": {
                        **{
                            k: v
                            for k, v in graph_stats.items()
                            if k != "degree_sequence"
                        },
                        "degree_sequence": list(graph_stats["degree_sequence"]),
                    },
                    "surface": surface_stats,
                }
            )
    return records


def html_escape(text: Any) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def family_hover(
    records: list[PhaseCell],
    lambdas: np.ndarray,
    thresholds: tuple[float, ...],
    stable_grid: np.ndarray,
) -> list[list[str]]:
    lookup = {(round(r.threshold_c, 12), round(r.lam, 12)): r for r in records}
    hover = []
    for row, threshold_c in enumerate(thresholds):
        hover_row = []
        for col, lam in enumerate(lambdas):
            r = lookup[(round(float(threshold_c), 12), round(float(lam), 12))]
            poly = r.polynomial or "(not computed; large core or error)"
            hover_row.append(
                "<br>".join(
                    [
                        f"<b>{html_escape(r.title)}</b>",
                        f"lambda={r.lam:.4f}",
                        f"c={r.threshold_c:.6g}",
                        f"phase id={int(stable_grid[row, col])}",
                        f"source={html_escape(r.source)}",
                        f"Yamada={html_escape(poly)}",
                        f"nodes={r.nodes}, edges={r.edges}, beta={r.cycle_rank}",
                        (
                            "interior components="
                            f"{r.interior_components}, handle rank={r.handle_rank}"
                        ),
                        f"volume fraction={r.interior_fraction:.4f}",
                        (
                            f"surface vertices={r.surface_points}, faces={r.surface_cells}, "
                            f"closed={r.surface_is_closed}"
                        ),
                        f"box boundary touched={r.touches_boundary}",
                        f"geometry={html_escape(Path(r.geometry_path).name)}",
                    ]
                )
            )
        hover.append(hover_row)
    return hover


def discrete_colorscale(n: int) -> list[list[Any]]:
    colors = [PHASE_COLORS[i % len(PHASE_COLORS)] for i in range(max(1, n))]
    if n <= 1:
        return [[0, colors[0]], [1, colors[0]]]
    scale: list[list[Any]] = []
    for idx, color in enumerate(colors):
        lo = idx / n
        hi = (idx + 1) / n
        scale.append([lo, color])
        scale.append([hi, color])
    return scale


def write_plotly_html(
    output_dir: Path,
    families: tuple[TPMSFamily, ...],
    records_by_family: dict[str, list[PhaseCell]],
    lambdas: np.ndarray,
    stable_by_family: dict[str, np.ndarray],
    labels_by_family: dict[str, dict[int, str]],
) -> Path | None:
    try:
        import plotly.graph_objects as go
        import plotly.io as pio
    except Exception:
        return None

    sections = []
    for index, family in enumerate(families, start=1):
        records = records_by_family[family.key]
        stable_grid = stable_by_family[family.key]
        labels = labels_by_family[family.key]
        n = int(stable_grid.max())
        fig = go.Figure(
            data=[
                go.Heatmap(
                    x=lambdas,
                    y=list(family.thresholds),
                    z=stable_grid,
                    zmin=1,
                    zmax=max(1, n),
                    colorscale=discrete_colorscale(n),
                    colorbar={
                        "title": {"text": "signature"},
                        "tickmode": "array",
                        "tickvals": list(range(1, n + 1)),
                    },
                    hoverinfo="text",
                    text=family_hover(records, lambdas, family.thresholds, stable_grid),
                )
            ]
        )
        fig.update_layout(
            title=f"{family.title}: {family.convention}",
            xaxis_title="lambda",
            yaxis_title="level-set threshold c",
            height=470,
            margin={"l": 72, "r": 42, "t": 70, "b": 60},
            font={
                "family": "Times New Roman, Times, serif",
                "size": 16,
                "color": "black",
            },
            paper_bgcolor="white",
            plot_bgcolor="white",
        )
        fig.update_xaxes(
            showline=True,
            linewidth=1.4,
            linecolor="black",
            mirror=True,
            ticks="outside",
        )
        fig.update_yaxes(
            showline=True,
            linewidth=1.4,
            linecolor="black",
            mirror=True,
            ticks="outside",
        )
        fragment = pio.to_html(
            fig,
            include_plotlyjs=(index == 1),
            full_html=False,
            config={"responsive": True},
        )
        legend_items = "".join(
            f"<li><b>{phase_id}</b>: {html_escape(label)}</li>"
            for phase_id, label in sorted(labels.items())
        )
        sections.append(
            f"""
            <section class="panel">
              <div class="panel-header">
                <h2>{html_escape(family.title)}</h2>
                <p><code>F_lambda = (1-lambda) {html_escape(family.start.key)} + lambda {html_escape(family.end.key)}</code></p>
                <p>{html_escape(family.rationale)}</p>
              </div>
              {fragment}
              <details>
                <summary>Yamada-distinguished signature legend</summary>
                <ul>{legend_items}</ul>
              </details>
            </section>
            """
        )

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Compact TPMS Yamada-Distinguished Phase Maps</title>
  <style>
    body {{ margin: 0; font-family: "Times New Roman", Times, serif; background: #ffffff; color: #111827; }}
    main {{ max-width: 1260px; margin: 0 auto; padding: 24px 28px 48px; }}
    h1 {{ font-size: 30px; margin: 0 0 8px; }}
    h2 {{ font-size: 23px; margin: 0; }}
    p {{ margin: 6px 0; font-size: 16px; line-height: 1.35; }}
    code {{ font-family: Menlo, Monaco, Consolas, monospace; font-size: 13px; }}
    .lede {{ margin-bottom: 22px; }}
    .panel {{ border: 1.5px solid #111827; padding: 16px 16px 12px; margin: 22px 0; }}
    .panel-header {{ margin-bottom: 8px; }}
    summary {{ cursor: pointer; font-size: 16px; font-weight: 700; margin-top: 8px; }}
    ul {{ columns: 2; padding-left: 22px; margin-top: 8px; }}
    li {{ break-inside: avoid; margin: 4px 0; font-size: 14px; }}
  </style>
</head>
<body>
  <main>
    <h1>Compact TPMS Yamada-Distinguished Phase Maps</h1>
    <p class="lede">Analytic porous-material TPMS fields are intersected with a finite outer design domain, sampled as compact scalar level-set solids, and passed through the same KnottedGraph volume skeletonization, embedded graph extraction, and Yamada/audit classification route as the Hamiltonian phase maps. Matching Yamada signatures are reported as Yamada-distinguished phases, not as complete topological equivalence claims.</p>
    {"".join(sections)}
  </main>
</body>
</html>
"""
    html_path = output_dir / "tpms_compact_scaffold_phase_maps.html"
    html_path.write_text(html, encoding="utf-8")
    print("wrote:", html_path)
    return html_path


def write_static_overview(
    output_dir: Path,
    families: tuple[TPMSFamily, ...],
    lambdas: np.ndarray,
    stable_by_family: dict[str, np.ndarray],
) -> tuple[Path, Path]:
    fig, axes = plt.subplots(
        1,
        len(families),
        figsize=(4.4 * len(families), 4.25),
        constrained_layout=True,
        facecolor="white",
    )
    if len(families) == 1:
        axes = [axes]
    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    for ax, family, panel in zip(axes, families, panel_labels):
        labels = stable_by_family[family.key]
        n = max(1, int(labels.max()))
        cmap = ListedColormap([PHASE_COLORS[i % len(PHASE_COLORS)] for i in range(n)])
        norm = BoundaryNorm(np.arange(0.5, n + 1.5), cmap.N)
        ax.imshow(
            labels,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            extent=(
                float(lambdas[0]),
                float(lambdas[-1]),
                float(family.thresholds[0]),
                float(family.thresholds[-1]),
            ),
            cmap=cmap,
            norm=norm,
        )
        ax.set_title(family.title, fontsize=14, fontweight="semibold", pad=8)
        ax.set_xlabel(r"$\lambda$", fontsize=16)
        ax.set_ylabel(r"$c$", fontsize=16)
        ax.tick_params(axis="both", labelsize=11, width=1.1, length=4)
        for spine in ax.spines.values():
            spine.set_linewidth(1.3)
        ax.text(
            -0.16,
            1.05,
            panel,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=15,
            fontweight="bold",
            color=TEXT_COLOR,
            clip_on=False,
        )
    fig.suptitle(
        "Compact finite TPMS scaffolds: Yamada-distinguished signatures",
        fontsize=16,
        fontweight="semibold",
    )
    png_path = output_dir / "tpms_compact_scaffold_phase_maps_overview.png"
    pdf_path = output_dir / "tpms_compact_scaffold_phase_maps_overview.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print("wrote:", png_path)
    print("wrote:", pdf_path)
    return png_path, pdf_path


def surface_mesh_for_plot(
    values: np.ndarray,
    level: float,
    *,
    spacing: np.ndarray,
    origin: np.ndarray,
    max_faces: int,
) -> tuple[np.ndarray, np.ndarray]:
    verts, faces, _, _ = marching_cubes(
        values,
        level=float(level),
        spacing=tuple(float(v) for v in spacing),
    )
    verts = verts + origin[None, :]
    if len(faces) > max_faces:
        indices = np.linspace(0, len(faces) - 1, int(max_faces), dtype=int)
        faces = faces[indices]
    return verts, faces


def write_endpoint_geometry_figure(
    output_dir: Path,
    *,
    dimension: int,
    threshold_c: float = 0.0,
    domain_kind: str = "sphere",
    span_half_width: float = 2.25 * math.pi,
    domain_radius_fraction: float = 0.72,
) -> tuple[Path, Path]:
    fields = implicit_fields()
    span = cube_span(float(span_half_width))
    domain = compact_domain(
        kind=domain_kind,
        span=span,
        radius_fraction=float(domain_radius_fraction),
    )
    _, x, y, z = sample_grid(span, dimension)
    domain_values = domain.sample(x, y, z)
    dummy_region = ImplicitSolidRegion(
        np.zeros((dimension, dimension, dimension), dtype=float),
        threshold_c=0.0,
        span=span,
    )
    fig = plt.figure(figsize=(10.2, 3.6), facecolor="white", constrained_layout=True)
    colors = {
        "gyroid": "#0f766e",
        "schwarz_p": "#2563eb",
        "diamond": "#dc2626",
    }
    for index, (key, field) in enumerate(fields.items(), start=1):
        ax = fig.add_subplot(1, 3, index, projection="3d")
        raw_values = field.sample(x, y, z)
        values = np.maximum(raw_values - float(threshold_c), domain_values)
        verts, faces = surface_mesh_for_plot(
            values,
            0.0,
            spacing=dummy_region.spacing,
            origin=dummy_region.origin,
            max_faces=50000,
        )
        mesh = Poly3DCollection(
            verts[faces],
            alpha=0.92,
            linewidths=0.0,
            edgecolors=(1, 1, 1, 0.0),
        )
        mesh.set_facecolor(colors[key])
        ax.add_collection3d(mesh)
        lo = float(span[0][0])
        hi = float(span[0][1])
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_zlim(lo, hi)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=22, azim=36)
        ax.set_title(field.title, fontsize=13, fontweight="semibold", pad=4)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_alpha(0.0)
            axis.line.set_color((0.3, 0.3, 0.3, 0.25))
    fig.suptitle(
        f"Compact finite TPMS scaffold endpoints at c={threshold_c:g}",
        fontsize=15,
        fontweight="semibold",
    )
    png_path = output_dir / "tpms_compact_endpoint_geometries.png"
    pdf_path = output_dir / "tpms_compact_endpoint_geometries.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print("wrote:", png_path)
    print("wrote:", pdf_path)
    return png_path, pdf_path


def write_records(
    output_dir: Path,
    records: list[PhaseCell],
    families: tuple[TPMSFamily, ...],
    lambdas: np.ndarray,
    map_info: dict[str, Any],
    endpoint_verification: list[dict[str, Any]],
    source_data: dict[str, Any],
    scan_parameters: dict[str, Any] | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "tpms_parameter_phase_map_records.json"
    csv_path = output_dir / "tpms_parameter_phase_map_records.csv"
    summary_path = output_dir / "tpms_parameter_phase_map_summary.json"
    source_data_path = output_dir / "tpms_parameter_phase_map_source_data.json"

    rows = [record_dict(r) for r in records]
    json_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "knotted_graph_root": package_location(),
        "scan_parameters": scan_parameters or {},
        "pipeline_reuse": {
            "skeletonization": "knotted_graph.extraction.skeletonize_volume",
            "graph_extraction": "knotted_graph.extraction.skeleton_image_to_graph",
            "graph_postprocessing": [
                "knotted_graph.core.remove_leaf_nodes",
                "knotted_graph.core.simplify_edges",
                "knotted_graph.core.smooth_edges",
            ],
            "yamada": "knotted_graph.applications.phase_maps._compute_yamada",
        },
        "scientific_scope_note": (
            "Cells sharing a Yamada value are described only as "
            "Yamada-distinguished signatures/phases. The run does not claim "
            "complete topological equivalence without stronger validation."
        ),
        "literature_basis": [
            {
                "claim": (
                    "Gyroid, Schwarz primitive, and Schwarz diamond are standard "
                    "TPMS scaffold fields used for porous-material design."
                ),
                "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8156807/",
            },
            {
                "claim": (
                    "Finite TPMS scaffold models should be closed by combining "
                    "the implicit surface with an outer boundary/design domain."
                ),
                "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5943328/",
            },
        ],
        "lambda_samples": len(lambdas),
        "lambdas": [float(v) for v in lambdas],
        "families": [
            {
                "key": f.key,
                "title": f.title,
                "start_field": f.start.key,
                "end_field": f.end.key,
                "start_formula": f.start.formula,
                "end_formula": f.end.formula,
                "thresholds": list(f.thresholds),
                "dimension": f.dimension,
                "span": f.span,
                "compact_domain": {
                    "key": f.domain.key,
                    "title": f.domain.title,
                    "formula": f.domain.formula,
                    "reference_note": f.domain.reference_note,
                },
                "convention": f.convention,
                "rationale": f.rationale,
            }
            for f in families
        ],
        "compactness_audit": {
            "reconstructed_volumes": len(records),
            "box_boundary_touch_count": sum(int(r.touches_boundary) for r in records),
            "closed_surface_count": sum(int(r.surface_is_closed) for r in records),
            "open_surface_count": sum(int(not r.surface_is_closed) for r in records),
            "surface_nonmanifold_edge_max": max(
                (int(r.surface_nonmanifold_edges) for r in records),
                default=0,
            ),
            "surface_open_edge_max": max(
                (int(r.surface_open_edges) for r in records),
                default=0,
            ),
        },
        "endpoint_verification": endpoint_verification,
        "map_info": map_info,
        "source_counts": dict(Counter(r.source for r in records)),
        "errors": [record_dict(r) for r in records if r.error],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    source_data_path.write_text(json.dumps(source_data, indent=2), encoding="utf-8")
    print("wrote:", json_path)
    print("wrote:", csv_path)
    print("wrote:", summary_path)
    print("wrote:", source_data_path)


def build_source_data(
    families: tuple[TPMSFamily, ...],
    lambdas: np.ndarray,
    raw_by_family: dict[str, np.ndarray],
    stable_by_family: dict[str, np.ndarray],
    labels_by_family: dict[str, dict[int, str]],
    signature_ids_by_family: dict[str, dict[str, int]],
    endpoint_verification: list[dict[str, Any]],
    timing: dict[str, Any],
) -> dict[str, Any]:
    source_data: dict[str, Any] = {
        "lambdas": [float(v) for v in lambdas],
        "endpoint_verification": endpoint_verification,
        "timing": timing,
        "phase_maps": {},
    }
    for family in families:
        reverse = {
            str(phase_id): signature
            for signature, phase_id in signature_ids_by_family[family.key].items()
        }
        source_data["phase_maps"][family.key] = {
            "title": family.title,
            "start_field": family.start.key,
            "end_field": family.end.key,
            "thresholds": [float(v) for v in family.thresholds],
            "compact_domain": {
                "key": family.domain.key,
                "title": family.domain.title,
                "formula": family.domain.formula,
            },
            "convention": family.convention,
            "raw_grid": raw_by_family[family.key].astype(int).tolist(),
            "stable_grid": stable_by_family[family.key].astype(int).tolist(),
            "display_reassigned_mask": (
                stable_by_family[family.key] != raw_by_family[family.key]
            ).tolist(),
            "display_filter_role": "visualization only; raw_grid is the classification record",
            "display_filter_rule": (
                "at most 3 synchronous passes; four-connected components; "
                "unique boundary-contact plurality among adjacent large components; "
                "ties unchanged; errors protected"
            ),
            "phase_labels": {
                str(key): value for key, value in labels_by_family[family.key].items()
            },
            "phase_signatures": reverse,
        }
    return source_data


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("_build/new_phase_maps/tpms"),
    )
    parser.add_argument("--dimension", type=int, default=64)
    parser.add_argument("--lambda-count", type=int, default=9)
    parser.add_argument("--threshold-count", type=int, default=9)
    parser.add_argument("--threshold-min", type=float, default=-0.40)
    parser.add_argument("--threshold-max", type=float, default=0.40)
    parser.add_argument(
        "--domain",
        choices=("sphere", "cylinder"),
        default="sphere",
        help="finite outer design domain used to close the scaffold inside the box",
    )
    parser.add_argument(
        "--span-half-width",
        type=float,
        default=2.25 * math.pi,
        help="half-width of the sampled cubic box in TPMS phase coordinates",
    )
    parser.add_argument(
        "--domain-radius-fraction",
        type=float,
        default=0.72,
        help="outer-domain radius as a fraction of the sampled half-width",
    )
    parser.add_argument("--max-exact-yamada-edges", type=int, default=18)
    parser.add_argument("--max-surface-sample-points", type=int, default=256)
    parser.add_argument("--max-edge-points", type=int, default=32)
    parser.add_argument("--min-stable-cells", type=int, default=1)
    parser.add_argument(
        "--only",
        nargs="*",
        choices=[
            family.key
            for family in tpms_families(
                dimension=64,
                thresholds=tuple(np.linspace(-0.40, 0.40, 9)),
            )
        ],
        default=None,
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.dimension < 8:
        raise ValueError("--dimension must be at least 8")
    if args.lambda_count < 2:
        raise ValueError("--lambda-count must be at least 2")
    if args.threshold_count < 2:
        raise ValueError("--threshold-count must be at least 2")
    thresholds = tuple(
        float(v)
        for v in np.linspace(
            float(args.threshold_min),
            float(args.threshold_max),
            int(args.threshold_count),
        )
    )
    lambdas = np.linspace(0.0, 1.0, int(args.lambda_count))
    families = tpms_families(
        dimension=int(args.dimension),
        thresholds=thresholds,
        domain_kind=str(args.domain),
        span_half_width=float(args.span_half_width),
        domain_radius_fraction=float(args.domain_radius_fraction),
    )
    if args.only:
        allowed = set(args.only)
        families = tuple(f for f in families if f.key in allowed)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    all_records: list[PhaseCell] = []
    records_by_family: dict[str, list[PhaseCell]] = {}
    raw_by_family: dict[str, np.ndarray] = {}
    stable_by_family: dict[str, np.ndarray] = {}
    labels_by_family: dict[str, dict[int, str]] = {}
    signature_ids_by_family: dict[str, dict[str, int]] = {}
    map_info: dict[str, Any] = {}

    for family in families:
        print(f"\n[{family.key}] {family.title}", flush=True)
        records = compute_family(
            family,
            lambdas,
            output_dir=args.output_dir,
            max_exact_yamada_edges=int(args.max_exact_yamada_edges),
            max_surface_sample_points=int(args.max_surface_sample_points),
            max_edge_points=int(args.max_edge_points),
        )
        raw_grid, label_lookup, signature_to_id = phase_grid(
            records,
            lambdas,
            family.thresholds,
        )
        stable_grid, changed_cells = stable_labels(
            raw_grid,
            min_cells=int(args.min_stable_cells),
            protected_labels={signature_to_id[r.phase_signature]
                              for r in records if r.source == "error"},
        )
        records_by_family[family.key] = records
        raw_by_family[family.key] = raw_grid
        stable_by_family[family.key] = stable_grid
        labels_by_family[family.key] = label_lookup
        signature_ids_by_family[family.key] = signature_to_id
        all_records.extend(records)
        map_info[family.key] = {
            "raw_phase_count": int(raw_grid.max()),
            "stable_phase_count": int(len(set(int(v) for v in stable_grid.ravel()))),
            "stable_reassigned_cells": int(changed_cells),
            "source_counts": dict(Counter(r.source for r in records)),
            "distinct_signatures": int(len({r.phase_signature for r in records})),
            "reconstructed_volumes": int(len(records)),
        }
        print(
            f"  phases: raw={map_info[family.key]['raw_phase_count']}, "
            f"stable={map_info[family.key]['stable_phase_count']}, "
            f"sources={map_info[family.key]['source_counts']}",
            flush=True,
        )

    png_path, pdf_path = write_static_overview(
        args.output_dir,
        families,
        lambdas,
        stable_by_family,
    )
    endpoint_png, endpoint_pdf = write_endpoint_geometry_figure(
        args.output_dir,
        dimension=int(args.dimension),
        threshold_c=0.0,
        domain_kind=str(args.domain),
        span_half_width=float(args.span_half_width),
        domain_radius_fraction=float(args.domain_radius_fraction),
    )
    html_path = write_plotly_html(
        args.output_dir,
        families,
        records_by_family,
        lambdas,
        stable_by_family,
        labels_by_family,
    )

    elapsed = time.perf_counter() - started
    timing = timing_summary(all_records, elapsed)
    map_info["timing"] = timing
    map_info["static_png"] = str(png_path)
    map_info["static_pdf"] = str(pdf_path)
    map_info["endpoint_png"] = str(endpoint_png)
    map_info["endpoint_pdf"] = str(endpoint_pdf)
    map_info["interactive_html"] = str(html_path) if html_path is not None else None

    endpoint_verification = verify_endpoints(families, endpoint_threshold=0.0)
    source_data = build_source_data(
        families,
        lambdas,
        raw_by_family,
        stable_by_family,
        labels_by_family,
        signature_ids_by_family,
        endpoint_verification,
        timing,
    )
    write_records(
        args.output_dir,
        all_records,
        families,
        lambdas,
        map_info,
        endpoint_verification,
        source_data,
        scan_parameters={
            "domain_kind": args.domain,
            "span_half_width": args.span_half_width,
            "domain_radius_fraction": args.domain_radius_fraction,
        },
    )

    print(
        f"Total reconstructed volumes: {len(all_records)} in {elapsed:.1f}s",
        flush=True,
    )
    print(
        "Yamada note: equal signatures are Yamada-distinguished phase labels, "
        "not complete topological-equivalence proofs.",
        flush=True,
    )


if __name__ == "__main__":
    main()
