#!/usr/bin/env python3
"""Material-Hamiltonian lambda scans for topology/Yamada phase maps.

Each family starts at the parameter set recovered from the old Multiband.ipynb
notebook.  Lambda changes one selected Hamiltonian coefficient, and the second
axis is the material band-gap threshold E.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/codex-matplotlib")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import sympy as sp
from matplotlib.colors import BoundaryNorm, ListedColormap
from skimage.measure import euler_number, label as label_volume


THIS_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(THIS_REPO / "tmp"))

from reproduce_multiband_material_surfaces import (  # noqa: E402
    D6_BASE,
    D6_PANELS,
    H_Co2MnGa_TB6_sympy,
    H_D6_sympy,
    H_Ti3Al_sympy,
    H_YH3_sympy,
    K_SYMBOLS,
)

KG_ROOT = Path(
    "/Users/hakanakgun/Desktop/Projects/ProfLeeProjects/"
    "Knotted_graph_code_paper/KnottedGraph_1earlier"
)
sys.path.insert(0, str(KG_ROOT / "src"))

from knotted_graph.applications.material_surface import MaterialFermiSurface  # noqa: E402
from knotted_graph.applications.phase_maps import (  # noqa: E402
    _compute_yamada,
    _graph_summary,
    _one_vertex_graph,
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
)

DEFAULT_MATERIAL_KEYS = ("tib2_d6_F", "co2mnga_t8")
DEFAULT_MIN_STABLE_CELLS = int(os.environ.get("MATERIAL_PHASE_MIN_ISLAND_CELLS", "8"))
DEFAULT_ISLAND_MERGE_STRATEGY = os.environ.get("MATERIAL_PHASE_ISLAND_MERGE", "below").strip().lower()

PHASE_SIGNATURE_MERGES: dict[str, dict[str, str]] = {
    "tib2_d6_F": {
        "yamada:Mul(Pow(Symbol('Y'), Integer(-7)), Pow(Add(Pow(Symbol('Y'), Integer(2)), Integer(1)), Integer(2)), Add(Pow(Symbol('Y'), Integer(2)), Symbol('Y'), Integer(1)), Add(Pow(Symbol('Y'), Integer(8)), Pow(Symbol('Y'), Integer(7)), Mul(Integer(5), Pow(Symbol('Y'), Integer(6))), Mul(Integer(2), Pow(Symbol('Y'), Integer(5))), Mul(Integer(10), Pow(Symbol('Y'), Integer(4))), Mul(Integer(2), Pow(Symbol('Y'), Integer(3))), Mul(Integer(5), Pow(Symbol('Y'), Integer(2))), Symbol('Y'), Integer(1)))":
            "yamada:Mul(Pow(Symbol('Y'), Integer(-11)), Add(Pow(Symbol('Y'), Integer(2)), Symbol('Y'), Integer(1)), Add(Pow(Symbol('Y'), Integer(20)), Mul(Integer(3), Pow(Symbol('Y'), Integer(19))), Mul(Integer(20), Pow(Symbol('Y'), Integer(18))), Mul(Integer(41), Pow(Symbol('Y'), Integer(17))), Mul(Integer(153), Pow(Symbol('Y'), Integer(16))), Mul(Integer(222), Pow(Symbol('Y'), Integer(15))), Mul(Integer(607), Pow(Symbol('Y'), Integer(14))), Mul(Integer(633), Pow(Symbol('Y'), Integer(13))), Mul(Integer(1367), Pow(Symbol('Y'), Integer(12))), Mul(Integer(1052), Pow(Symbol('Y'), Integer(11))), Mul(Integer(1789), Pow(Symbol('Y'), Integer(10))), Mul(Integer(1052), Pow(Symbol('Y'), Integer(9))), Mul(Integer(1367), Pow(Symbol('Y'), Integer(8))), Mul(Integer(633), Pow(Symbol('Y'), Integer(7))), Mul(Integer(607), Pow(Symbol('Y'), Integer(6))), Mul(Integer(222), Pow(Symbol('Y'), Integer(5))), Mul(Integer(153), Pow(Symbol('Y'), Integer(4))), Mul(Integer(41), Pow(Symbol('Y'), Integer(3))), Mul(Integer(20), Pow(Symbol('Y'), Integer(2))), Mul(Integer(3), Symbol('Y')), Integer(1)))",
        "vertex:nodes=5;components=5;interior_components=5;touches_boundary=1;yamada=Integer(-1)":
            "vertex:nodes=1;components=1;interior_components=1;touches_boundary=1;yamada=Integer(-1)",
    },
    "co2mnga_t8": {
        "core:ic=1;h=40;b=1;nodes=68;edges=107;gcomp=1;beta=40;deg=((3, 67), (13, 1))":
            "core:ic=1;h=40;b=1;nodes=69;edges=108;gcomp=1;beta=40;deg=((3, 68), (12, 1))",
        "core:ic=10;h=13;b=1;nodes=25;edges=36;gcomp=2;beta=13;deg=((0, 1), (3, 24))":
            "core:ic=10;h=13;b=1;nodes=24;edges=36;gcomp=1;beta=13;deg=((3, 24),)",
    },
}


@dataclass(frozen=True)
class MaterialFamily:
    key: str
    title: str
    builder: Callable[[dict[str, float]], sp.Matrix]
    base_params: dict[str, float]
    lambda_param: str
    param_start: float
    param_end: float
    energies: tuple[float, ...]
    span: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    band_pair: tuple[int, int] = (0, 1)
    dimension: int = 80
    rationale: str = ""
    landmark_energies: tuple[float, ...] = ()

    def params_at(self, lam: float) -> dict[str, float]:
        params = dict(self.base_params)
        params[self.lambda_param] = self.param_start + float(lam) * (
            self.param_end - self.param_start
        )
        return params

    def hamiltonian_at(self, lam: float) -> sp.Matrix:
        return self.builder(self.params_at(lam))

    def param_value_at(self, lam: float) -> float:
        return self.param_start + float(lam) * (self.param_end - self.param_start)


@dataclass
class PhaseCell:
    material: str
    title: str
    lam: float
    parameter_name: str
    parameter_value: float
    energy: float
    dimension: int
    span: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    band_pair: tuple[int, int]
    phase_signature: str
    phase_label: str
    source: str
    polynomial: str
    nodes: int
    edges: int
    components: int
    cycle_rank: int
    degree_sequence: tuple[int, ...]
    interior_voxels: int
    interior_components: int
    euler_characteristic: int
    handle_rank: int
    touches_boundary: bool
    boundary_faces: tuple[str, ...]
    exact_yamada_attempted: bool
    error: str | None = None
    classification_computed: bool = True


def cube_span(half_width: float) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    return ((-float(half_width), float(half_width)),) * 3


def energy_axis(start: float, stop: float, step: float = 0.01) -> tuple[float, ...]:
    values = []
    current = float(start)
    while current <= float(stop) + 1e-12:
        values.append(round(current, 12))
        current += float(step)
    if not values or not math.isclose(values[-1], float(stop), rel_tol=0.0, abs_tol=1e-10):
        values.append(round(float(stop), 12))
    return tuple(values)


def material_families(dimension: int) -> tuple[MaterialFamily, ...]:
    d6_base = {**D6_BASE, **D6_PANELS["a"]}
    ti3al_base = {
        "A1": -9.66,
        "A2": 11.37,
        "B1": 36.22,
        "B2": -25.71,
        "M1": 0.12,
        "M2": -0.52,
        "C": 22.34,
    }
    yh3_base = {
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
    co2_base = {
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
    return (
        MaterialFamily(
            key="tib2_d6_F",
            title="TiB2",
            builder=H_D6_sympy,
            base_params=d6_base,
            lambda_param="F",
            param_start=d6_base["F"],
            param_end=D6_PANELS["c"]["F"],
            energies=energy_axis(0.180, 5.20, 0.01),
            span=cube_span(1.5),
            band_pair=(0, 1),
            dimension=dimension,
            rationale=(
                "F is the high-order off-diagonal coupling; the phase map follows the legacy "
                "Multiband.ipynb band pair (0,1), k-window [-1.5,1.5]^3, and the E range "
                "covering the reported D6 transitions near 2.8, 3.9, and 4.2 eV, shown from 0.18 to 5.2 eV."
            ),
            landmark_energies=(0.18, 2.8, 3.65, 3.9, 4.2, 4.25, 4.6, 5.0, 5.2),
        ),
        MaterialFamily(
            key="ti3al_M2",
            title="Ti3Al",
            builder=H_Ti3Al_sympy,
            base_params=ti3al_base,
            lambda_param="M2",
            param_start=ti3al_base["M2"],
            param_end=-0.10,
            energies=energy_axis(0.05, 0.30, 0.01),
            span=((-0.45, 0.45), (-0.45, 0.45), (-0.08, 0.08)),
            dimension=dimension,
            rationale="M2 is a mass offset in the two-band effective model; a focused k-window is needed because the low-energy pocket is tiny inside the old [-pi,pi]^3 box.",
            landmark_energies=(0.25, 0.3),
        ),
        MaterialFamily(
            key="yh3_m1",
            title="YH3",
            builder=H_YH3_sympy,
            base_params=yh3_base,
            lambda_param="m1",
            param_start=yh3_base["m1"],
            param_end=2.75,
            energies=energy_axis(0.002, 0.20, 0.01),
            span=cube_span(0.5),
            dimension=dimension,
            rationale="m1 shifts the first cosine constraint and changes the small-energy nodal network; the E range is extended for diagnostics, but YH3 is omitted from the default paper HTML when not requested.",
            landmark_energies=(0.006, 0.2),
        ),
        MaterialFamily(
            key="co2mnga_t8",
            title="Co2MnGa",
            builder=H_Co2MnGa_TB6_sympy,
            base_params=co2_base,
            lambda_param="t8",
            param_start=co2_base["t8"],
            param_end=-0.70,
            energies=energy_axis(0.125, 1.60, 0.01),
            span=cube_span(2.05 * math.pi),
            band_pair=(0, 1),
            dimension=dimension,
            rationale=(
                "t8 is the d-p hybridization coupling; the phase map follows the legacy "
                "Multiband.ipynb band pair (0,1), full [-2pi,2pi]^3 Brillouin-zone window "
                "with a small padding, and the E range covering the 0.25 and 1.15 eV transitions, "
                "now shown from 0.125 to 1.6 eV."
            ),
            landmark_energies=(0.25, 1.15, 1.4, 1.6),
        ),
    )


def reset_gap_threshold(obj: MaterialFermiSurface, energy: float) -> None:
    obj.gap_tol = float(energy)
    obj.skeleton_graph_cache = None
    obj.skeleton_graph_cache_args = None
    obj._pv_data_args = None
    for name in ("_skeleton_image", "fields_pv", "exceptional_surface_pv"):
        obj.__dict__.pop(name, None)


def interior_summary(mask: np.ndarray) -> dict[str, Any]:
    mask = np.asarray(mask, dtype=bool)
    boundary_faces = []
    for label, touched in (
        ("kx_min", mask[0, :, :].any()),
        ("kx_max", mask[-1, :, :].any()),
        ("ky_min", mask[:, 0, :].any()),
        ("ky_max", mask[:, -1, :].any()),
        ("kz_min", mask[:, :, 0].any()),
        ("kz_max", mask[:, :, -1].any()),
    ):
        if bool(touched):
            boundary_faces.append(label)
    if not mask.any():
        return {
            "interior_voxels": 0,
            "interior_components": 0,
            "euler_characteristic": 0,
            "handle_rank": 0,
            "touches_boundary": False,
            "boundary_faces": tuple(),
        }
    _, component_count = label_volume(mask, connectivity=3, return_num=True)
    euler = int(euler_number(mask, connectivity=3))
    handle_rank = max(0, int(component_count) - euler)
    return {
        "interior_voxels": int(mask.sum()),
        "interior_components": int(component_count),
        "euler_characteristic": euler,
        "handle_rank": int(handle_rank),
        "touches_boundary": bool(boundary_faces),
        "boundary_faces": tuple(boundary_faces),
    }


def graph_signature(summary: dict[str, Any], interior: dict[str, Any]) -> str:
    degree_hist = tuple(sorted(Counter(summary["degree_sequence"]).items()))
    return (
        "core:"
        f"ic={interior['interior_components']};"
        f"h={interior['handle_rank']};"
        f"b={int(interior['touches_boundary'])};"
        f"nodes={summary['nodes']};"
        f"edges={summary['edges']};"
        f"gcomp={summary['components']};"
        f"beta={summary['cycle_rank']};"
        f"deg={degree_hist}"
    )


def genus_zero_vertex_graph(component_count: int) -> nx.MultiGraph:
    graph = nx.MultiGraph()
    for node in range(max(1, int(component_count))):
        graph.add_node(node, pos=(float(node), 0.0, 0.0))
    return graph


def vertex_phase_signature(
    polynomial_expr: sp.Expr,
    summary: dict[str, Any],
    interior: dict[str, Any],
) -> str:
    return (
        "vertex:"
        f"nodes={summary['nodes']};"
        f"components={summary['components']};"
        f"interior_components={interior['interior_components']};"
        f"touches_boundary={int(interior['touches_boundary'])};"
        f"yamada={sp.srepr(polynomial_expr)}"
    )


def boundary_open_signature(summary: dict[str, Any], interior: dict[str, Any]) -> str:
    degree_hist = tuple(sorted(Counter(summary["degree_sequence"]).items()))
    return (
        "boundary-open:"
        f"ic={interior['interior_components']};"
        f"h={interior['handle_rank']};"
        f"faces={','.join(interior['boundary_faces'])};"
        f"nodes={summary['nodes']};"
        f"edges={summary['edges']};"
        f"gcomp={summary['components']};"
        f"beta={summary['cycle_rank']};"
        f"deg={degree_hist}"
    )


def short_signature_label(signature: str, source: str, polynomial: str) -> str:
    if source == "vertex":
        return f"vertex Yamada {polynomial}"
    if source == "yamada":
        return f"Yamada {polynomial}"
    if source == "boundary-open":
        chunks = dict(
            part.split("=", 1)
            for part in signature.removeprefix("boundary-open:").split(";")
            if "=" in part
        )
        faces = chunks.get("faces", "") or "none"
        return (
            f"boundary-open ic={chunks.get('ic','?')}, h={chunks.get('h','?')}, "
            f"faces={faces}, core V={chunks.get('nodes','?')}, E={chunks.get('edges','?')}, "
            f"beta={chunks.get('beta','?')}"
        )
    if source == "large-core":
        chunks = dict(
            part.split("=", 1)
            for part in signature.removeprefix("core:").split(";")
            if "=" in part
        )
        return (
            f"core V={chunks.get('nodes','?')}, E={chunks.get('edges','?')}, "
            f"beta={chunks.get('beta','?')}, h={chunks.get('h','?')}, "
            f"ic={chunks.get('ic','?')}"
        )
    if source == "error":
        return "extraction error"
    return source


def evaluate_cell(
    family: MaterialFamily,
    obj: MaterialFermiSurface,
    lam: float,
    energy: float,
    *,
    max_exact_yamada_edges: int,
) -> PhaseCell:
    reset_gap_threshold(obj, energy)
    mask = np.asarray(obj._interior_mask, dtype=bool)
    interior = interior_summary(mask)

    if interior["handle_rank"] == 0:
        graph = genus_zero_vertex_graph(interior["interior_components"])
        summary = _graph_summary(graph)
        polynomial_expr = _compute_yamada(graph, Y, {"normalize": True})
        polynomial_expr = sp.factor(sp.together(sp.expand(polynomial_expr)))
        polynomial = str(polynomial_expr)
        signature = vertex_phase_signature(polynomial_expr, summary, interior)
        source = "vertex"
        attempted = True
        error = None
    else:
        graph = nx.MultiGraph()
        polynomial = ""
        attempted = False
        error = None
        try:
            graph = obj.skeleton_graph(
                smooth_epsilon=0,
                simplify=True,
                force_small_edge_contraction=True,
                small_edge_limit=math.pi * 0.1,
                previous_n_edgepoint=20,
            )
            summary = _graph_summary(graph)
            if graph.number_of_edges() == 0:
                graph = genus_zero_vertex_graph(interior["interior_components"])
                summary = _graph_summary(graph)
                polynomial_expr = _compute_yamada(graph, Y, {"normalize": True})
                polynomial_expr = sp.factor(sp.together(sp.expand(polynomial_expr)))
                polynomial = str(polynomial_expr)
                signature = vertex_phase_signature(polynomial_expr, summary, interior)
                source = "vertex"
                attempted = True
            elif graph.number_of_edges() <= max_exact_yamada_edges:
                attempted = True
                polynomial_expr = _compute_yamada(graph, Y, {"normalize": True})
                polynomial_expr = sp.factor(sp.together(sp.expand(polynomial_expr)))
                polynomial = str(polynomial_expr)
                signature = "yamada:" + sp.srepr(polynomial_expr)
                source = "yamada"
            else:
                signature = graph_signature(summary, interior)
                source = "large-core"
        except Exception as exc:
            message = str(exc)
            graph = _one_vertex_graph()
            summary = _graph_summary(graph)
            if any(
                text in message
                for text in (
                    "graph has no edges",
                    "skeleton image is empty",
                    "does not contain any True voxels",
                    "Skeletonization produced no points",
                    "collapsed to fewer than two distinct points",
                )
            ):
                if interior["handle_rank"] == 0:
                    graph = genus_zero_vertex_graph(interior["interior_components"])
                    summary = _graph_summary(graph)
                    polynomial_expr = _compute_yamada(graph, Y, {"normalize": True})
                    polynomial_expr = sp.factor(sp.together(sp.expand(polynomial_expr)))
                    polynomial = str(polynomial_expr)
                    signature = vertex_phase_signature(polynomial_expr, summary, interior)
                    source = "vertex"
                    attempted = True
                else:
                    graph = nx.MultiGraph()
                    summary = _graph_summary(graph)
                    polynomial = "boundary-open"
                    signature = boundary_open_signature(summary, interior)
                    source = "boundary-open"
                    attempted = False
                    error = (
                        "Skeleton extraction failed for a boundary-touching nonzero-genus mask. "
                        f"Boundary faces: {', '.join(interior['boundary_faces'])}"
                    )
            else:
                polynomial = ""
                signature = "error:" + type(exc).__name__ + ":" + message[:120]
                source = "error"
                error = type(exc).__name__ + ": " + message[:240]

    phase_label = short_signature_label(signature, source, polynomial)
    if source == "vertex":
        phase_label = f"vertex x{summary['nodes']}; Yamada {polynomial}"

    return PhaseCell(
        material=family.key,
        title=family.title,
        lam=float(lam),
        parameter_name=family.lambda_param,
        parameter_value=float(family.param_value_at(lam)),
        energy=float(energy),
        dimension=family.dimension,
        span=family.span,
        band_pair=family.band_pair,
        phase_signature=signature,
        phase_label=phase_label,
        source=source,
        polynomial=polynomial,
        nodes=int(summary["nodes"]),
        edges=int(summary["edges"]),
        components=int(summary["components"]),
        cycle_rank=int(summary["cycle_rank"]),
        degree_sequence=tuple(int(value) for value in summary["degree_sequence"]),
        exact_yamada_attempted=attempted,
        error=error,
        classification_computed=True,
        **interior,
    )


def adaptive_probe_indices(family: MaterialFamily, adaptive_energy_step: float) -> list[int]:
    energies = [float(value) for value in family.energies]
    if len(energies) <= 2 or adaptive_energy_step <= 0:
        return list(range(len(energies)))
    diffs = [
        abs(energies[index + 1] - energies[index])
        for index in range(len(energies) - 1)
        if not math.isclose(energies[index + 1], energies[index])
    ]
    native_step = min(diffs) if diffs else float(adaptive_energy_step)
    stride = max(1, int(round(float(adaptive_energy_step) / native_step)))
    indices = {0, len(energies) - 1}
    indices.update(range(0, len(energies), stride))
    for landmark in family.landmark_energies:
        if energies[0] - 1e-12 <= float(landmark) <= energies[-1] + 1e-12:
            nearest = min(range(len(energies)), key=lambda idx: abs(energies[idx] - float(landmark)))
            indices.add(nearest)
            if nearest > 0:
                indices.add(nearest - 1)
            if nearest + 1 < len(energies):
                indices.add(nearest + 1)
    return sorted(indices)


def adaptive_fill_record(record: PhaseCell, energy: float) -> PhaseCell:
    return dataclasses.replace(
        record,
        energy=float(energy),
        classification_computed=False,
    )


def evaluate_energy_grid_adaptive(
    family: MaterialFamily,
    obj: MaterialFermiSurface,
    lam: float,
    *,
    max_exact_yamada_edges: int,
    adaptive_energy_step: float,
) -> list[PhaseCell]:
    energies = [float(value) for value in family.energies]
    if adaptive_energy_step <= 0:
        return [
            evaluate_cell(
                family,
                obj,
                float(lam),
                float(energy),
                max_exact_yamada_edges=max_exact_yamada_edges,
            )
            for energy in energies
        ]

    computed: dict[int, PhaseCell] = {}

    def ensure(index: int) -> PhaseCell:
        if index not in computed:
            computed[index] = evaluate_cell(
                family,
                obj,
                float(lam),
                float(energies[index]),
                max_exact_yamada_edges=max_exact_yamada_edges,
            )
        return computed[index]

    probes = adaptive_probe_indices(family, adaptive_energy_step)
    for index in probes:
        ensure(index)

    for left, right in zip(probes, probes[1:]):
        if right - left <= 1:
            continue
        if ensure(left).phase_signature != ensure(right).phase_signature:
            for index in range(left + 1, right):
                ensure(index)

    sorted_computed = sorted(computed)
    records: list[PhaseCell] = []
    for index, energy in enumerate(energies):
        if index in computed:
            records.append(computed[index])
            continue
        left = [probe for probe in sorted_computed if probe < index]
        right = [probe for probe in sorted_computed if probe > index]
        source_index = left[-1] if left else right[0]
        records.append(adaptive_fill_record(ensure(source_index), energy))
    return records


def compute_family(
    family: MaterialFamily,
    lambdas: np.ndarray,
    *,
    chunk_size: int,
    max_exact_yamada_edges: int,
    workers: int = 1,
    adaptive_energy_step: float = 0.05,
) -> list[PhaseCell]:
    if int(workers) > 1:
        return compute_family_parallel(
            family,
            lambdas,
            chunk_size=chunk_size,
            max_exact_yamada_edges=max_exact_yamada_edges,
            workers=int(workers),
            adaptive_energy_step=float(adaptive_energy_step),
        )

    records: list[PhaseCell] = []
    for col, lam in enumerate(lambdas):
        t0 = time.perf_counter()
        obj = MaterialFermiSurface(
            family.hamiltonian_at(float(lam)),
            k_symbols=K_SYMBOLS,
            span=family.span,
            dimension=family.dimension,
            band_pair=family.band_pair,
            gap_tol=float(family.energies[0]),
            check_hermitian=False,
            check_pt_symmetry=False,
            chunk_size=chunk_size,
            force_small_edge_contraction=True,
            small_edge_limit=math.pi * 0.1,
            previous_n_edgepoint=20,
        )
        gap = obj.band_gap
        volume_elapsed = time.perf_counter() - t0
        print(
            f"  {family.key} lambda {col + 1:02d}/{len(lambdas)} "
            f"{family.lambda_param}={family.param_value_at(lam):.6g}: "
            f"gap [{float(np.nanmin(gap)):.5g}, {float(np.nanmax(gap)):.5g}] "
            f"in {volume_elapsed:.2f}s",
            flush=True,
        )
        records.extend(
            evaluate_energy_grid_adaptive(
                family,
                obj,
                float(lam),
                max_exact_yamada_edges=max_exact_yamada_edges,
                adaptive_energy_step=float(adaptive_energy_step),
            )
        )
    return records


def compute_lambda_slice(
    family: MaterialFamily,
    lam: float,
    col: int,
    total: int,
    *,
    chunk_size: int,
    max_exact_yamada_edges: int,
    adaptive_energy_step: float,
) -> tuple[int, list[PhaseCell], dict[str, Any]]:
    t0 = time.perf_counter()
    obj = MaterialFermiSurface(
        family.hamiltonian_at(float(lam)),
        k_symbols=K_SYMBOLS,
        span=family.span,
        dimension=family.dimension,
        band_pair=family.band_pair,
        gap_tol=float(family.energies[0]),
        check_hermitian=False,
        check_pt_symmetry=False,
        chunk_size=chunk_size,
        force_small_edge_contraction=True,
        small_edge_limit=math.pi * 0.1,
        previous_n_edgepoint=20,
    )
    gap = obj.band_gap
    volume_elapsed = time.perf_counter() - t0
    records = evaluate_energy_grid_adaptive(
        family,
        obj,
        float(lam),
        max_exact_yamada_edges=max_exact_yamada_edges,
        adaptive_energy_step=float(adaptive_energy_step),
    )
    elapsed = time.perf_counter() - t0
    return (
        int(col),
        records,
        {
            "col": int(col),
            "total": int(total),
            "lambda": float(lam),
            "parameter_value": float(family.param_value_at(lam)),
            "gap_min": float(np.nanmin(gap)),
            "gap_max": float(np.nanmax(gap)),
            "volume_elapsed": round(volume_elapsed, 3),
            "elapsed": round(elapsed, 3),
            "source_counts": dict(Counter(record.source for record in records)),
            "distinct_signatures": int(len({record.phase_signature for record in records})),
            "computed_cuts": int(sum(record.classification_computed for record in records)),
            "filled_cuts": int(sum(not record.classification_computed for record in records)),
        },
    )


def compute_family_parallel(
    family: MaterialFamily,
    lambdas: np.ndarray,
    *,
    chunk_size: int,
    max_exact_yamada_edges: int,
    workers: int,
    adaptive_energy_step: float,
) -> list[PhaseCell]:
    records_by_col: dict[int, list[PhaseCell]] = {}
    total = len(lambdas)
    max_workers = max(1, int(workers))
    print(f"  using {max_workers} worker processes for lambda columns", flush=True)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                compute_lambda_slice,
                family,
                float(lam),
                col,
                total,
                chunk_size=chunk_size,
                max_exact_yamada_edges=max_exact_yamada_edges,
                adaptive_energy_step=float(adaptive_energy_step),
            )
            for col, lam in enumerate(lambdas)
        ]
        for future in as_completed(futures):
            col, col_records, info = future.result()
            records_by_col[col] = col_records
            print(
                f"  {family.key} lambda {info['col'] + 1:02d}/{info['total']} "
                f"{family.lambda_param}={info['parameter_value']:.6g}: "
                f"gap [{info['gap_min']:.5g}, {info['gap_max']:.5g}] "
                f"volume {info['volume_elapsed']:.2f}s total {info['elapsed']:.2f}s "
                f"computed={info['computed_cuts']} filled={info['filled_cuts']} "
                f"sources={info['source_counts']} phases={info['distinct_signatures']}",
                flush=True,
            )
    records: list[PhaseCell] = []
    for col in range(total):
        records.extend(records_by_col[col])
    return records


def contiguous_energy_runs(energies: tuple[float, ...], native_step: float = 0.01) -> list[tuple[float, ...]]:
    runs: list[list[float]] = []
    for energy in sorted(float(value) for value in energies):
        if not runs or not math.isclose(
            energy - runs[-1][-1], float(native_step), rel_tol=0.0, abs_tol=1e-10
        ):
            runs.append([energy])
        else:
            runs[-1].append(energy)
    return [tuple(run) for run in runs]


def compute_missing_lambda_slice(
    family: MaterialFamily,
    missing_energies: tuple[float, ...],
    lam: float,
    col: int,
    total: int,
    *,
    chunk_size: int,
    max_exact_yamada_edges: int,
    adaptive_energy_step: float,
) -> tuple[int, list[PhaseCell], dict[str, Any]]:
    t0 = time.perf_counter()
    obj = MaterialFermiSurface(
        family.hamiltonian_at(float(lam)),
        k_symbols=K_SYMBOLS,
        span=family.span,
        dimension=family.dimension,
        band_pair=family.band_pair,
        gap_tol=float(min(missing_energies)),
        check_hermitian=False,
        check_pt_symmetry=False,
        chunk_size=chunk_size,
        force_small_edge_contraction=True,
        small_edge_limit=math.pi * 0.1,
        previous_n_edgepoint=20,
    )
    gap = obj.band_gap
    volume_elapsed = time.perf_counter() - t0
    records: list[PhaseCell] = []
    for run in contiguous_energy_runs(missing_energies):
        run_family = dataclasses.replace(
            family,
            energies=run,
            landmark_energies=tuple(
                value for value in family.landmark_energies
                if run[0] - 1e-12 <= float(value) <= run[-1] + 1e-12
            ),
        )
        records.extend(
            evaluate_energy_grid_adaptive(
                run_family,
                obj,
                float(lam),
                max_exact_yamada_edges=max_exact_yamada_edges,
                adaptive_energy_step=float(adaptive_energy_step),
            )
        )
    elapsed = time.perf_counter() - t0
    return (
        int(col),
        records,
        {
            "col": int(col),
            "total": int(total),
            "parameter_value": float(family.param_value_at(lam)),
            "gap_min": float(np.nanmin(gap)),
            "gap_max": float(np.nanmax(gap)),
            "volume_elapsed": round(volume_elapsed, 3),
            "elapsed": round(elapsed, 3),
            "computed_cuts": int(sum(record.classification_computed for record in records)),
            "filled_cuts": int(sum(not record.classification_computed for record in records)),
        },
    )


def compute_missing_family_parallel(
    family: MaterialFamily,
    lambdas: np.ndarray,
    missing_energies: tuple[float, ...],
    *,
    chunk_size: int,
    max_exact_yamada_edges: int,
    workers: int,
    adaptive_energy_step: float,
) -> list[PhaseCell]:
    records_by_col: dict[int, list[PhaseCell]] = {}
    total = len(lambdas)
    max_workers = max(1, int(workers))
    print(
        f"  extending {len(missing_energies)} E rows with {max_workers} worker processes",
        flush=True,
    )
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                compute_missing_lambda_slice,
                family,
                missing_energies,
                float(lam),
                col,
                total,
                chunk_size=chunk_size,
                max_exact_yamada_edges=max_exact_yamada_edges,
                adaptive_energy_step=float(adaptive_energy_step),
            )
            for col, lam in enumerate(lambdas)
        ]
        for future in as_completed(futures):
            col, col_records, info = future.result()
            records_by_col[col] = col_records
            print(
                f"  {family.key} extension lambda {info['col'] + 1:02d}/{info['total']} "
                f"{family.lambda_param}={info['parameter_value']:.6g}: "
                f"gap [{info['gap_min']:.5g}, {info['gap_max']:.5g}] "
                f"volume {info['volume_elapsed']:.2f}s total {info['elapsed']:.2f}s "
                f"computed={info['computed_cuts']} filled={info['filled_cuts']}",
                flush=True,
            )
    records: list[PhaseCell] = []
    for col in range(total):
        records.extend(records_by_col[col])
    return records


def connected_components_for_label(labels: np.ndarray, phase_id: int) -> list[list[tuple[int, int]]]:
    h, w = labels.shape
    seen = np.zeros(labels.shape, dtype=bool)
    components: list[list[tuple[int, int]]] = []
    for i in range(h):
        for j in range(w):
            if seen[i, j] or int(labels[i, j]) != int(phase_id):
                continue
            queue = deque([(i, j)])
            seen[i, j] = True
            comp: list[tuple[int, int]] = []
            while queue:
                ci, cj = queue.popleft()
                comp.append((ci, cj))
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ni, nj = ci + di, cj + dj
                    if 0 <= ni < h and 0 <= nj < w and not seen[ni, nj] and int(labels[ni, nj]) == int(phase_id):
                        seen[ni, nj] = True
                        queue.append((ni, nj))
            components.append(comp)
    return components


def label_components(labels: np.ndarray) -> list[tuple[int, list[tuple[int, int]]]]:
    components: list[tuple[int, list[tuple[int, int]]]] = []
    for phase_id in sorted(set(int(v) for v in labels.ravel())):
        for component in connected_components_for_label(labels, phase_id):
            components.append((phase_id, component))
    return components


def nearest_large_phase(
    component: list[tuple[int, int]],
    phase_id: int,
    large_components: list[tuple[int, list[tuple[int, int]]]],
) -> int | None:
    if not large_components:
        return None
    center = np.asarray(component, dtype=float).mean(axis=0)
    best: tuple[float, int, int] | None = None
    for target_phase, target_component in large_components:
        if int(target_phase) == int(phase_id):
            continue
        target = np.asarray(target_component, dtype=float)
        distance = float(np.min(np.sum((target - center) ** 2, axis=1)))
        candidate = (distance, -len(target_component), int(target_phase))
        if best is None or candidate < best:
            best = candidate
    if best is not None:
        return best[2]

    target_phase, _ = min(
        large_components,
        key=lambda item: (
            float(np.sum((np.asarray(item[1], dtype=float).mean(axis=0) - center) ** 2)),
            -len(item[1]),
            int(item[0]),
        ),
    )
    return int(target_phase)


def closest_lower_large_phase(
    component: list[tuple[int, int]],
    phase_id: int,
    large_components: list[tuple[int, list[tuple[int, int]]]],
) -> int | None:
    if not large_components:
        return None
    component_arr = np.asarray(component, dtype=float)
    bottom_row = int(min(row for row, _ in component))
    best: tuple[float, float, int, int] | None = None
    for target_phase, target_component in large_components:
        if int(target_phase) == int(phase_id):
            continue
        target_arr = np.asarray(target_component, dtype=float)
        below = target_arr[target_arr[:, 0] < bottom_row]
        if below.size == 0:
            continue
        deltas = component_arr[:, None, :] - below[None, :, :]
        distances = np.sum(deltas * deltas, axis=2)
        flat_index = int(np.argmin(distances))
        comp_index, target_index = np.unravel_index(flat_index, distances.shape)
        vertical_gap = float(component_arr[comp_index, 0] - below[target_index, 0])
        distance = float(distances[comp_index, target_index])
        candidate = (distance, vertical_gap, -len(target_component), int(target_phase))
        if best is None or candidate < best:
            best = candidate
    return None if best is None else int(best[3])


def stable_labels(
    labels: np.ndarray,
    *,
    min_cells: int,
    merge_strategy: str = DEFAULT_ISLAND_MERGE_STRATEGY,
) -> tuple[np.ndarray, int]:
    stable = labels.copy()
    changed_cells = 0
    if int(min_cells) <= 1:
        return stable, changed_cells
    strategy = str(merge_strategy or "below").strip().lower()

    for _ in range(max(1, stable.size)):
        components = label_components(stable)
        large_components = [
            (phase_id, component)
            for phase_id, component in components
            if len(component) >= int(min_cells)
        ]
        if not large_components:
            break
        large_phase_ids = {int(phase_id) for phase_id, _ in large_components}
        changed = False
        for phase_id, component in components:
            if len(component) >= int(min_cells):
                continue
            neighbor_counts: Counter[int] = Counter()
            if strategy != "below":
                for i, j in component:
                    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        ni, nj = i + di, j + dj
                        if 0 <= ni < stable.shape[0] and 0 <= nj < stable.shape[1]:
                            neighbor = int(stable[ni, nj])
                            if neighbor != int(phase_id) and neighbor in large_phase_ids:
                                neighbor_counts[neighbor] += 1
            if strategy == "below":
                replacement = closest_lower_large_phase(component, int(phase_id), large_components)
                if replacement is None:
                    replacement = nearest_large_phase(component, int(phase_id), large_components)
            elif neighbor_counts:
                replacement = neighbor_counts.most_common(1)[0][0]
            else:
                replacement = nearest_large_phase(component, int(phase_id), large_components)
            if replacement is None or int(replacement) == int(phase_id):
                continue
            for i, j in component:
                stable[i, j] = int(replacement)
                changed_cells += 1
            changed = True
        if not changed:
            break
    return stable, changed_cells


def apply_signature_phase_merges(
    labels: np.ndarray,
    signature_to_id: dict[str, int],
    family_key: str,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    merged = labels.copy()
    applied: list[dict[str, Any]] = []
    for source_signature, target_signature in PHASE_SIGNATURE_MERGES.get(str(family_key), {}).items():
        source_id = signature_to_id.get(source_signature)
        target_id = signature_to_id.get(target_signature)
        if source_id is None or target_id is None:
            continue
        cell_count = int(np.count_nonzero(merged == int(source_id)))
        if cell_count == 0:
            continue
        merged[merged == int(source_id)] = int(target_id)
        applied.append(
            {
                "source_signature": source_signature,
                "target_signature": target_signature,
                "source_raw_id": int(source_id),
                "target_raw_id": int(target_id),
                "cells": cell_count,
            }
        )
    return merged, applied


def remap_labels_contiguous(
    labels: np.ndarray,
    label_lookup: dict[int, str],
) -> tuple[np.ndarray, dict[int, str]]:
    used = sorted(set(int(value) for value in labels.ravel()))
    old_to_new = {old: new for new, old in enumerate(used, start=1)}
    remapped = np.vectorize(old_to_new.__getitem__)(labels).astype(int)
    remapped_lookup = {
        new: label_lookup.get(old, f"phase {old}")
        for old, new in old_to_new.items()
    }
    return remapped, remapped_lookup


def phase_grid(records: list[PhaseCell], lambdas: np.ndarray, energies: tuple[float, ...]) -> tuple[np.ndarray, dict[int, str], dict[str, int]]:
    signatures: list[str] = []
    seen: set[str] = set()
    for record in records:
        if record.phase_signature not in seen:
            signatures.append(record.phase_signature)
            seen.add(record.phase_signature)
    signature_to_id = {signature: idx + 1 for idx, signature in enumerate(signatures)}
    label_lookup = {idx + 1: next(r.phase_label for r in records if r.phase_signature == signature) for idx, signature in enumerate(signatures)}
    lookup = {(round(r.energy, 12), round(r.lam, 12)): r for r in records}
    grid = np.zeros((len(energies), len(lambdas)), dtype=int)
    for row, energy in enumerate(energies):
        for col, lam in enumerate(lambdas):
            record = lookup[(round(float(energy), 12), round(float(lam), 12))]
            grid[row, col] = signature_to_id[record.phase_signature]
    return grid, label_lookup, signature_to_id


def record_dict(record: PhaseCell) -> dict[str, Any]:
    row = dataclasses.asdict(record)
    row["degree_sequence"] = list(record.degree_sequence)
    row["boundary_faces"] = list(record.boundary_faces)
    return row


def phase_cell_from_dict(row: dict[str, Any]) -> PhaseCell:
    data = dict(row)
    data["span"] = tuple(tuple(float(value) for value in pair) for pair in data["span"])
    data["band_pair"] = tuple(int(value) for value in data["band_pair"])
    data["degree_sequence"] = tuple(int(value) for value in data["degree_sequence"])
    data["boundary_faces"] = tuple(str(value) for value in data["boundary_faces"])
    return PhaseCell(**data)


def write_records(output_dir: Path, records: list[PhaseCell], families: tuple[MaterialFamily, ...], lambdas: np.ndarray, map_info: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "material_parameter_phase_map_records.json"
    csv_path = output_dir / "material_parameter_phase_map_records.csv"
    summary_path = output_dir / "material_parameter_phase_map_summary.json"
    json_path.write_text(json.dumps([record_dict(r) for r in records], indent=2), encoding="utf-8")
    rows = [record_dict(r) for r in records]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "knotted_graph_root": str(KG_ROOT),
        "material_api": "knotted_graph.applications.material_surface.MaterialFermiSurface",
        "lambda_samples": len(lambdas),
        "lambdas": [float(v) for v in lambdas],
        "families": [
            {
                "key": f.key,
                "title": f.title,
                "lambda_param": f.lambda_param,
                "param_start": f.param_start,
                "param_end": f.param_end,
                "dimension": f.dimension,
                "span": f.span,
                "band_pair": f.band_pair,
                "energies": list(f.energies),
                "rationale": f.rationale,
            }
            for f in families
        ],
        "map_info": map_info,
        "source_counts": dict(Counter(r.source for r in records)),
        "errors": [record_dict(r) for r in records if r.error],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("wrote:", json_path)
    print("wrote:", csv_path)
    print("wrote:", summary_path)


def family_hover(records: list[PhaseCell], lambdas: np.ndarray, energies: tuple[float, ...], stable_grid: np.ndarray) -> list[list[str]]:
    lookup = {(round(r.energy, 12), round(r.lam, 12)): r for r in records}
    hover = []
    for row, energy in enumerate(energies):
        hover_row = []
        for col, lam in enumerate(lambdas):
            r = lookup[(round(float(energy), 12), round(float(lam), 12))]
            poly = r.polynomial or "(not computed; large core)"
            hover_row.append(
                "<br>".join(
                    [
                        f"<b>{r.title}</b>",
                        f"lambda={r.lam:.4f}",
                        f"{r.parameter_name}={r.parameter_value:.6g}",
                        f"E={r.energy:.6g} eV",
                        f"phase id={int(stable_grid[row, col])}",
                        f"classification={'computed' if r.classification_computed else 'adaptive fill'}",
                        f"source={r.source}",
                        f"Yamada={poly}",
                        f"nodes={r.nodes}, edges={r.edges}, beta={r.cycle_rank}",
                        f"interior components={r.interior_components}, handle rank={r.handle_rank}",
                        f"touches boundary={r.touches_boundary}",
                    ]
                )
            )
        hover.append(hover_row)
    return hover


def discrete_colorscale(n: int) -> list[list[Any]]:
    colors = [PHASE_COLORS[i % len(PHASE_COLORS)] for i in range(max(1, n))]
    if n <= 1:
        return [[0, colors[0]], [1, colors[0]]]
    scale = []
    for idx, color in enumerate(colors):
        lo = idx / n
        hi = (idx + 1) / n
        scale.append([lo, color])
        scale.append([hi, color])
    return scale


def write_plotly_html(
    output_dir: Path,
    families: tuple[MaterialFamily, ...],
    records_by_family: dict[str, list[PhaseCell]],
    lambdas: np.ndarray,
    stable_by_family: dict[str, np.ndarray],
    labels_by_family: dict[str, dict[int, str]],
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
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
                    y=list(family.energies),
                    z=stable_grid,
                    zmin=1,
                    zmax=max(1, n),
                    colorscale=discrete_colorscale(n),
                    xgap=1,
                    ygap=1,
                    colorbar={"title": {"text": "phase"}, "tickmode": "array", "tickvals": list(range(1, n + 1))},
                    hoverinfo="text",
                    text=family_hover(records, lambdas, family.energies, stable_grid),
                )
            ]
        )
        fig.update_layout(
            title=f"{family.title}: lambda varies {family.lambda_param}",
            xaxis_title="lambda",
            yaxis_title="E (eV)",
            height=480,
            margin={"l": 70, "r": 40, "t": 70, "b": 60},
            font={"family": "Times New Roman, Times, serif", "size": 16, "color": "black"},
            paper_bgcolor="white",
            plot_bgcolor="white",
        )
        fig.update_xaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, ticks="outside")
        fig.update_yaxes(showline=True, linewidth=1.4, linecolor="black", mirror=True, ticks="outside")
        fragment = pio.to_html(fig, include_plotlyjs=(index == 1), full_html=False, config={"responsive": True})
        legend_items = "".join(
            f"<li><b>{phase_id}</b>: {label}</li>" for phase_id, label in sorted(labels.items())
        )
        sections.append(
            f"""
            <section class="panel">
              <div class="panel-header">
                <h2>{family.title}</h2>
                <p><code>{family.lambda_param}(lambda) = {family.param_start:.6g} + lambda * ({family.param_end - family.param_start:.6g})</code></p>
                <p>{family.rationale}</p>
              </div>
              {fragment}
              <details open>
                <summary>Phase legend</summary>
                <ul>{legend_items}</ul>
              </details>
            </section>
            """
        )
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Material Hamiltonian Lambda Phase Maps</title>
  <style>
    body {{ margin: 0; font-family: "Times New Roman", Times, serif; background: #ffffff; color: #111827; }}
    main {{ max-width: 1240px; margin: 0 auto; padding: 24px 28px 48px; }}
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
    <h1>Material Hamiltonian Lambda Phase Maps</h1>
    <p class="lede">Each panel starts from the old Multiband.ipynb material Hamiltonian and varies one coefficient as lambda. The vertical axis is the band-gap threshold E. Hover a cell for the parameter value, graph summary, and exact Yamada polynomial when it was tractable.</p>
    {''.join(sections)}
  </main>
</body>
</html>
"""
    html_path = output_dir / "material_parameter_phase_maps.html"
    html_path.write_text(html, encoding="utf-8")
    print("wrote:", html_path)
    return html_path


def write_static_overview(
    output_dir: Path,
    families: tuple[MaterialFamily, ...],
    lambdas: np.ndarray,
    stable_by_family: dict[str, np.ndarray],
) -> Path:
    family_count = len(families)
    ncols = 3 if family_count <= 3 else 2
    nrows = int(math.ceil(family_count / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.8 * ncols, 4.6 * nrows),
        constrained_layout=True,
        facecolor="white",
        squeeze=False,
    )
    for ax, family, panel in zip(axes.ravel(), families, ["(a)", "(b)", "(c)", "(d)"]):
        labels = stable_by_family[family.key]
        n = int(labels.max())
        cmap = ListedColormap([PHASE_COLORS[i % len(PHASE_COLORS)] for i in range(n)])
        norm = BoundaryNorm(np.arange(0.5, n + 1.5), cmap.N)
        ax.pcolormesh(lambdas, family.energies, labels, shading="nearest", cmap=cmap, norm=norm)
        for phase_id in sorted(set(int(value) for value in labels.ravel())):
            mask = (labels == phase_id).astype(float)
            if mask.min() != mask.max():
                ax.contour(lambdas, family.energies, mask, levels=[0.5], colors="black", linewidths=1.0)
        ax.set_title(
            rf"{family.title}: ${family.lambda_param}$",
            fontsize=17,
            fontweight="semibold",
            pad=10,
        )
        ax.set_xlabel(r"$\lambda$", fontsize=19)
        ax.set_ylabel(r"$E$ (eV)", fontsize=19)
        ax.tick_params(axis="both", labelsize=14, width=1.2, length=5)
        for spine in ax.spines.values():
            spine.set_linewidth(1.4)
        ax.text(
            -0.13,
            1.05,
            panel,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=18,
            fontweight="bold",
            color=TEXT_COLOR,
            clip_on=False,
        )
    for ax in axes.ravel()[family_count:]:
        ax.set_visible(False)
    png_path = output_dir / "material_parameter_phase_maps_overview.png"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print("wrote:", png_path)
    return png_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=THIS_REPO / "tmp" / "material_parameter_phase_maps",
    )
    parser.add_argument("--dimension", type=int, default=140)
    parser.add_argument("--lambda-count", type=int, default=60)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--max-exact-yamada-edges", type=int, default=24)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--adaptive-energy-step",
        type=float,
        default=0.05,
        help=(
            "Exact classify this approximate E stride, then refine every dense E row "
            "inside intervals where adjacent probes disagree. Use <=0 for exact every row."
        ),
    )
    parser.add_argument(
        "--only",
        nargs="*",
        choices=[family.key for family in material_families(80)],
        default=list(DEFAULT_MATERIAL_KEYS),
    )
    parser.add_argument("--min-stable-cells", type=int, default=DEFAULT_MIN_STABLE_CELLS)
    parser.add_argument(
        "--island-merge-strategy",
        choices=("below", "nearest"),
        default=DEFAULT_ISLAND_MERGE_STRATEGY if DEFAULT_ISLAND_MERGE_STRATEGY in {"below", "nearest"} else "below",
        help="How small connected phase islands are absorbed during post-processing.",
    )
    parser.add_argument(
        "--reuse-records",
        action="store_true",
        help="Regenerate HTML/PNG post-processing from existing records without recomputing N-grid classifications.",
    )
    parser.add_argument(
        "--extend-existing-records",
        action="store_true",
        help="Compute only missing E rows for selected families, preserve existing records, and rebuild all outputs.",
    )
    return parser.parse_args()


def extend_existing_dataset(args: argparse.Namespace) -> None:
    started = time.perf_counter()
    records_path = args.output_dir / "material_parameter_phase_map_records.json"
    summary_path = args.output_dir / "material_parameter_phase_map_summary.json"
    records_data = json.loads(records_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    lambdas = np.asarray(summary["lambdas"], dtype=float)

    family_dimensions = {
        str(item["key"]): int(item["dimension"])
        for item in summary["families"]
    }
    family_lookup = {
        family.key: family
        for dimension in sorted(set(family_dimensions.values()))
        for family in material_families(dimension)
    }
    family_keys = [str(item["key"]) for item in summary["families"]]
    families = tuple(family_lookup[key] for key in family_keys)
    selected = set(str(key) for key in (args.only or family_keys))

    records_by_family: dict[str, list[PhaseCell]] = defaultdict(list)
    for row in records_data:
        records_by_family[str(row["material"])].append(phase_cell_from_dict(row))

    for family in families:
        desired_energies = {round(float(value), 12) for value in family.energies}
        records_by_family[family.key] = [
            dataclasses.replace(record, title=family.title)
            for record in records_by_family[family.key]
            if round(float(record.energy), 12) in desired_energies
        ]
        if family.key not in selected:
            continue

        counts_by_energy = Counter(
            round(float(record.energy), 12)
            for record in records_by_family[family.key]
        )
        partial_rows = {
            energy: count for energy, count in counts_by_energy.items()
            if count != len(lambdas)
        }
        if partial_rows:
            raise RuntimeError(
                f"Cannot incrementally extend {family.key}: partial existing E rows {partial_rows}"
            )
        missing_energies = tuple(
            float(value) for value in family.energies
            if counts_by_energy.get(round(float(value), 12), 0) == 0
        )
        if not missing_energies:
            print(f"[{family.key}] no missing E rows", flush=True)
            continue
        new_records = compute_missing_family_parallel(
            family,
            lambdas,
            missing_energies,
            chunk_size=int(args.chunk_size),
            max_exact_yamada_edges=int(args.max_exact_yamada_edges),
            workers=int(args.workers),
            adaptive_energy_step=float(args.adaptive_energy_step),
        )
        records_by_family[family.key].extend(new_records)

    stable_by_family: dict[str, np.ndarray] = {}
    labels_by_family: dict[str, dict[int, str]] = {}
    map_info: dict[str, Any] = {}
    all_records: list[PhaseCell] = []
    for family in families:
        records = records_by_family[family.key]
        raw_grid, label_lookup, signature_to_id = phase_grid(records, lambdas, family.energies)
        stable_grid, changed_cells = stable_labels(
            raw_grid,
            min_cells=int(args.min_stable_cells),
            merge_strategy=str(args.island_merge_strategy),
        )
        stable_grid, manual_merges = apply_signature_phase_merges(
            stable_grid, signature_to_id, family.key
        )
        stable_grid, stable_label_lookup = remap_labels_contiguous(stable_grid, label_lookup)
        stable_by_family[family.key] = stable_grid
        labels_by_family[family.key] = stable_label_lookup
        all_records.extend(records)
        map_info[family.key] = {
            "raw_phase_count": int(len(set(int(value) for value in raw_grid.ravel()))),
            "stable_phase_count": int(len(set(int(value) for value in stable_grid.ravel()))),
            "stable_reassigned_cells": int(changed_cells),
            "stable_min_component_cells": int(args.min_stable_cells),
            "stable_island_merge_strategy": str(args.island_merge_strategy),
            "manual_signature_merges": manual_merges,
            "source_counts": dict(Counter(record.source for record in records)),
            "distinct_signatures": int(len({record.phase_signature for record in records})),
            "classification_computed_cells": int(sum(record.classification_computed for record in records)),
            "classification_adaptive_fill_cells": int(sum(not record.classification_computed for record in records)),
            "adaptive_energy_step": float(args.adaptive_energy_step),
        }
        print(
            f"[{family.key}] E={family.energies[0]:.3g}..{family.energies[-1]:.3g}, "
            f"rows={len(family.energies)}, stable phases={map_info[family.key]['stable_phase_count']}, "
            f"manual merges={sum(item['cells'] for item in manual_merges)} cells",
            flush=True,
        )

    html_path = write_plotly_html(
        args.output_dir, families, records_by_family, lambdas, stable_by_family, labels_by_family
    )
    png_path = write_static_overview(args.output_dir, families, lambdas, stable_by_family)
    map_info["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    map_info["html"] = str(html_path)
    map_info["png"] = str(png_path)
    map_info["incremental_extension"] = True
    write_records(args.output_dir, all_records, families, lambdas, map_info)
    print(f"Total incremental extension elapsed: {time.perf_counter() - started:.1f}s")


def main() -> None:
    args = parse_args()
    if args.reuse_records and args.extend_existing_records:
        raise ValueError("Choose either --reuse-records or --extend-existing-records, not both")
    if args.extend_existing_records:
        extend_existing_dataset(args)
        return
    if args.reuse_records:
        records_path = args.output_dir / "material_parameter_phase_map_records.json"
        summary_path = args.output_dir / "material_parameter_phase_map_summary.json"
        records_data = json.loads(records_path.read_text(encoding="utf-8"))
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        lambdas = np.asarray(summary["lambdas"], dtype=float)
        family_keys = [str(item["key"]) for item in summary["families"]]
        family_dimensions = {str(item["key"]): int(item["dimension"]) for item in summary["families"]}
        family_lookup = {
            family.key: family
            for dimension in sorted(set(family_dimensions.values()))
            for family in material_families(dimension)
        }
        families = tuple(family_lookup[key] for key in family_keys)
        records_by_family: dict[str, list[PhaseCell]] = defaultdict(list)
        for row in records_data:
            records_by_family[str(row["material"])].append(phase_cell_from_dict(row))
        stable_by_family: dict[str, np.ndarray] = {}
        labels_by_family: dict[str, dict[int, str]] = {}
        for family in families:
            raw_grid, label_lookup, signature_to_id = phase_grid(
                records_by_family[family.key], lambdas, family.energies
            )
            stable_grid, _ = stable_labels(
                raw_grid,
                min_cells=int(args.min_stable_cells),
                merge_strategy=str(args.island_merge_strategy),
            )
            stable_grid, _ = apply_signature_phase_merges(
                stable_grid, signature_to_id, family.key
            )
            stable_grid, stable_label_lookup = remap_labels_contiguous(stable_grid, label_lookup)
            stable_by_family[family.key] = stable_grid
            labels_by_family[family.key] = stable_label_lookup
        write_plotly_html(args.output_dir, families, records_by_family, lambdas, stable_by_family, labels_by_family)
        write_static_overview(args.output_dir, families, lambdas, stable_by_family)
        return

    lambdas = np.linspace(0.0, 1.0, int(args.lambda_count))
    families = material_families(int(args.dimension))
    if args.only:
        allowed = set(args.only)
        families = tuple(f for f in families if f.key in allowed)
    started = time.perf_counter()
    all_records: list[PhaseCell] = []
    records_by_family: dict[str, list[PhaseCell]] = {}
    stable_by_family: dict[str, np.ndarray] = {}
    labels_by_family: dict[str, dict[int, str]] = {}
    map_info: dict[str, Any] = {}
    for family in families:
        print(f"\n[{family.key}] {family.title}: varying {family.lambda_param}", flush=True)
        records = compute_family(
            family,
            lambdas,
            chunk_size=int(args.chunk_size),
            max_exact_yamada_edges=int(args.max_exact_yamada_edges),
            workers=int(args.workers),
            adaptive_energy_step=float(args.adaptive_energy_step),
        )
        raw_grid, label_lookup, signature_to_id = phase_grid(records, lambdas, family.energies)
        stable_grid, changed_cells = stable_labels(
            raw_grid,
            min_cells=int(args.min_stable_cells),
            merge_strategy=str(args.island_merge_strategy),
        )
        stable_grid, manual_merges = apply_signature_phase_merges(
            stable_grid, signature_to_id, family.key
        )
        stable_grid, stable_label_lookup = remap_labels_contiguous(stable_grid, label_lookup)
        records_by_family[family.key] = records
        stable_by_family[family.key] = stable_grid
        labels_by_family[family.key] = stable_label_lookup
        all_records.extend(records)
        map_info[family.key] = {
            "raw_phase_count": int(len(set(int(v) for v in raw_grid.ravel()))),
            "stable_phase_count": int(len(set(int(v) for v in stable_grid.ravel()))),
            "stable_reassigned_cells": int(changed_cells),
            "stable_min_component_cells": int(args.min_stable_cells),
            "stable_island_merge_strategy": str(args.island_merge_strategy),
            "manual_signature_merges": manual_merges,
            "source_counts": dict(Counter(r.source for r in records)),
            "distinct_signatures": int(len({r.phase_signature for r in records})),
            "classification_computed_cells": int(sum(r.classification_computed for r in records)),
            "classification_adaptive_fill_cells": int(sum(not r.classification_computed for r in records)),
            "adaptive_energy_step": float(args.adaptive_energy_step),
        }
        print(
            f"  phases: raw={map_info[family.key]['raw_phase_count']}, "
            f"stable={map_info[family.key]['stable_phase_count']}, "
            f"reassigned tiny cells={changed_cells}, sources={map_info[family.key]['source_counts']}",
            flush=True,
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    html_path = write_plotly_html(args.output_dir, families, records_by_family, lambdas, stable_by_family, labels_by_family)
    png_path = write_static_overview(args.output_dir, families, lambdas, stable_by_family)
    map_info["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    map_info["html"] = str(html_path)
    map_info["png"] = str(png_path)
    write_records(args.output_dir, all_records, families, lambdas, map_info)
    print(f"Total elapsed: {time.perf_counter() - started:.1f}s")


if __name__ == "__main__":
    main()
