#!/usr/bin/env python3
"""Execute every original map point without filtering or abstract fallback.

Hamiltonian plan: the notebook's exact np.linspace arrays, 60 lambdas, energy
prefix lengths 16/24/16/24/39, and dimension 120. TPMS plan: original 21x21
parameter arrays at an explicitly recorded voxel resolution. All source masks,
collapse witnesses and graph geometries are archived. Polynomial timeouts are
unavailable values, never phases. Cavity processing and formula code are not
modified. These are finite-domain digital-spine results, not inferred continuum
classifications. A second generic projection checks each subcubic evaluation.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
from functools import lru_cache
import gzip
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import signal
import subprocess
from time import monotonic
import warnings

import networkx as nx
import numpy as np
import sympy as sp

from full_map_kernel import reconstruct
from knotted_graph.applications.nodal.deformation import NodalBlochPath
from knotted_graph.applications.nodal.models import (
    hopf_link_bloch_vector, trefoil_bloch_vector, unknot_bloch_vector,
    solomon_bloch_vector, pq_torus_knot_bloch_vector,
)
from knotted_graph.applications.nodal.skeleton import NodalSkeleton
from knotted_graph.applications.phase_maps import volume_topology, _compute_yamada_audited
from knotted_graph.applications.phase_map_examples._tpms import implicit_fields
from knotted_graph.extraction.cubical_spine import sample_grid_breaks

HAMILTONIAN_FAMILIES = (
    ("hopf_to_trefoil", 16), ("hopf_to_solomon", 24),
    ("unknot_to_trefoil", 16), ("unknot_to_solomon", 24),
    ("trefoil_to_cinquefoil", 39),
)
TPMS_FAMILIES = ("schwarz_p_to_diamond", "gyroid_to_schwarz_p", "gyroid_to_diamond")


def dump(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def plan(mode):
    if mode == "hamiltonian":
        lambdas, levels = np.linspace(0., 1., 60), np.linspace(.30, 5.25, 50)
        families = HAMILTONIAN_FAMILIES
    else:
        lambdas, levels = np.linspace(0., 1., 21), np.linspace(0., .3, 21)
        families = tuple((key, 21) for key in TPMS_FAMILIES)
    return [{"id": f"{mode}_{key}_c{ci:03d}_l{li:03d}", "source": mode,
             "family": key, "lambda_index": li, "level_index": ci,
             "lambda": float(lam), "level": float(level),
             "lambda_hex": float(lam).hex(), "level_hex": float(level).hex()}
            for key, count in families for ci, level in enumerate(levels[:count])
            for li, lam in enumerate(lambdas)]


def cinquefoil(gamma):
    return pq_torus_knot_bloch_vector(2, 5, gamma)


@lru_cache(maxsize=256)
def endpoints(family, level):
    paths = {
        "hopf_to_trefoil": (hopf_link_bloch_vector, trefoil_bloch_vector),
        "hopf_to_solomon": (hopf_link_bloch_vector, solomon_bloch_vector),
        "unknot_to_trefoil": (unknot_bloch_vector, trefoil_bloch_vector),
        "unknot_to_solomon": (unknot_bloch_vector, solomon_bloch_vector),
        "trefoil_to_cinquefoil": (trefoil_bloch_vector, cinquefoil),
    }
    return NodalBlochPath(*paths[family]).endpoints(level)


def sample(case, dimension):
    if case["source"] == "hamiltonian":
        left, right = endpoints(case["family"], case["level"])
        lam = case["lambda"]
        vector = tuple(sp.expand((1 - lam) * a + lam * b) for a, b in zip(left, right))
        model = NodalSkeleton(char=vector, dimension=dimension)
        mask = np.asarray(model._interior_mask, dtype=bool)
        axes = (model.kx_vals, model.ky_vals, model.kz_vals)
        specification = {
            "domain": "closed_finite_box", "periodic_identification": False,
            "bounds": np.asarray(model.span).tolist(), "bloch_vector": [sp.srepr(x) for x in vector],
            "mask_rule": "unchanged NodalSkeleton._interior_mask",
        }
    else:
        fields = implicit_fields()
        names = {
            "schwarz_p_to_diamond": ("schwarz_p", "diamond"),
            "gyroid_to_schwarz_p": ("gyroid", "schwarz_p"),
            "gyroid_to_diamond": ("gyroid", "diamond"),
        }
        first, last = names[case["family"]]
        axis = np.linspace(-2.25 * np.pi, 2.25 * np.pi, dimension)
        axes = (axis, axis, axis)
        x, y, z = np.meshgrid(*axes, indexing="ij")
        radius = .72 * (2.25 * np.pi)
        value = ((1 - case["lambda"]) * fields[first].function(x, y, z)
                 + case["lambda"] * fields[last].function(x, y, z))
        mask = (value <= case["level"]) & (x*x + y*y + z*z <= radius*radius)
        specification = {"domain": "original_spherical_coupon", "radius_hex": float(radius).hex(),
                         "periodic_identification": False, "fields": [first, last]}
    specification["voxel_realization"] = "closed_dual_cells_clipped_at_sample_endpoints"
    return mask, axes, specification


def save_graph(path, graph):
    ids = {n: i for i, n in enumerate(graph)}
    value = {"nodes": [{"id": ids[n], "pos": list(d["pos"])} for n, d in graph.nodes(data=True)],
             "edges": [{"u": ids[u], "v": ids[v], "key": int(k), "pts": d["pts"].tolist()}
                       for u, v, k, d in graph.edges(keys=True, data=True)]}
    raw = gzip.compress(json.dumps(value, separators=(",", ":")).encode(), mtime=0)
    Path(path).write_bytes(raw)
    return hashlib.sha256(raw).hexdigest()


def evaluate(case, dimension, kernel, directory):
    os.setsid()
    directory = Path(directory)
    started = monotonic()
    row = {**case, "dimension": dimension, "status": "started", "yamada": None,
           "analytic_source_correspondence_certified": False, "source_phase_class_certified": False}
    try:
        mask, axes, specification = sample(case, dimension)
        row["source_specification"] = specification
        row["source_mask_sha256"] = hashlib.sha256(mask.astype(np.uint8).tobytes()).hexdigest()
        topology = volume_topology(mask)
        row["volume"] = asdict(topology)
        np.savez_compressed(directory / "source.npz", mask=mask, x=axes[0], y=axes[1], z=axes[2])
        row["status"] = "sampled"
        dump(directory / "record.json", row)
        if topology.enclosed_voids:
            row["status"] = "existing_cavity_route_not_revised"
            return
        if not mask.any():
            row["status"] = "empty_source"
            return
        lines = tuple(sample_grid_breaks(a) for a in axes)
        graph, evidence, _, _ = reconstruct(mask, kernel, gridlines=lines,
                                            archive=directory / "collapse.npz", timeout=120)
        row["reconstruction"] = evidence
        row["status"] = evidence["kind"]
        dump(directory / "record.json", row)
        if graph is None:
            return
        components = nx.number_connected_components(graph)
        row["graph"] = {"vertices": len(graph), "edges": graph.number_of_edges(),
                        "components": components, "cycle_rank": graph.number_of_edges() - len(graph) + components,
                        "max_degree": max(dict(graph.degree()).values(), default=0)}
        if (components, row["graph"]["cycle_rank"]) != (topology.connected_components, topology.handle_rank):
            raise RuntimeError("independent digital homology disagrees with the collapse endpoint")
        row["graph_sha256"] = save_graph(directory / "graph.json.gz", graph)
        row["status"] = "polynomial_pending"
        dump(directory / "record.json", row)
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            polynomial, audit = _compute_yamada_audited(graph, sp.Symbol("A"),
                                                      {"normalize": True, "n_jobs": 1})
            row["projection_checks"] = [{"value": str(polynomial), "audit": audit}]
            if audit["is_subcubic"] and graph.number_of_edges():
                second = None
                errors = []
                for angles in ((17., 31., 53.), (37., 61., 19.), (71., 23., 47.)):
                    try:
                        second, second_audit = _compute_yamada_audited(
                            graph, sp.Symbol("A"), {"normalize": True, "n_jobs": 1, "rotation_angles": angles})
                        break
                    except Exception as exc:
                        errors.append(f"{type(exc).__name__}: {exc}")
                row["secondary_projection_errors"] = errors
                if second is None:
                    raise RuntimeError("no independent generic projection completed")
                row["projection_checks"].append({"value": str(second), "audit": second_audit})
                if sp.expand(polynomial - second) != 0:
                    raise RuntimeError("normalized subcubic polynomial disagrees between generic projections")
            row["warnings"] = [str(w.message) for w in captured]
        row["yamada"] = str(polynomial)
        row["evaluation"] = audit
        row["signature"] = audit["evaluation_kind"] + ":" + sp.srepr(sp.expand(polynomial))
        row["status"] = "evaluated" if audit["is_subcubic"] else "fixed_diagram_evaluated"
    except Exception as exc:
        row["status"] = "unavailable"
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["yamada"] = None
    finally:
        row["seconds"] = monotonic() - started
        dump(directory / "record.json", row)


def run_cases(cases, dimension, kernel, out, seconds):
    if "fork" not in mp.get_all_start_methods():
        raise RuntimeError("this bounded full-map runner currently requires POSIX fork")
    context = mp.get_context("fork")
    rows = []
    for number, case in enumerate(cases):
        directory = out / case["id"]
        directory.mkdir()
        child = context.Process(target=evaluate, args=(case, dimension, kernel, str(directory)))
        child.start()
        child.join(seconds)
        if child.is_alive():
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                child.kill()
            child.join()
            partial = json.loads((directory / "record.json").read_text()) if (directory / "record.json").exists() else case.copy()
            row = {**partial, "status": "time_budget", "yamada": None,
                   "error": "case exceeded its frozen time budget", "case_seconds_limit": seconds}
            dump(directory / "record.json", row)
        elif (directory / "record.json").exists():
            row = json.loads((directory / "record.json").read_text())
            if child.exitcode:
                row.update(status="worker_error", yamada=None, worker_exitcode=child.exitcode)
                dump(directory / "record.json", row)
        else:
            row = {**case, "status": "worker_error", "yamada": None, "worker_exitcode": child.exitcode}
            dump(directory / "record.json", row)
        rows.append(row)
        dump(out / "records.json", rows)
        if number % 10 == 0 or number + 1 == len(cases):
            print("MAP_PROGRESS " + json.dumps({"completed": len(rows), "planned": len(cases),
                  "last": case["id"], "status": row["status"],
                  "counts": dict(Counter(r["status"] for r in rows))}), flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("hamiltonian", "tpms", "preflight"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--dimension", type=int, default=120)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--case-seconds", type=int, default=90)
    args = parser.parse_args()
    if args.out.exists() or args.dimension < 2 or not 0 <= args.shard < args.shards or args.case_seconds < 1:
        parser.error("new output directory and valid dimension/shard/time parameters required")
    args.out.mkdir(parents=True)
    if args.mode == "preflight":
        chosen = [("hopf_to_trefoil", 0, 0), ("hopf_to_trefoil", 30, 0),
                  ("hopf_to_trefoil", 30, 15), ("trefoil_to_cinquefoil", 59, 0)]
        cases = [r for r in plan("hamiltonian")
                 if (r["family"], r["lambda_index"], r["level_index"]) in chosen]
        cases += [r for r in plan("tpms") if r["lambda_index"] == 0 and r["level_index"] == 7][:2]
        all_cases = cases
    else:
        all_cases = plan(args.mode)
        cases = all_cases[args.shard::args.shards]
    dump(args.out / "plan.json", {"mode": args.mode, "all_planned_cells": len(all_cases),
         "shard": args.shard, "shards": args.shards, "cases": cases, "dimension": args.dimension,
         "case_seconds": args.case_seconds, "display_filter": False, "abstract_fallback": False,
         "source_scope": "finite_domain_voxel_complex_and_selected_embedded_graph",
         "all_word_formulas_changed": False, "cavity_treatment_changed": False})
    paths = [Path(__file__), Path(__file__).with_name("full_map_kernel.py"),
             Path(__file__).with_name("full_map_collapse.cpp"), Path("uv.lock")]
    dump(args.out / "provenance.json", {"commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
         "python": platform.python_version(), "platform": platform.platform(),
         "source_hashes": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
         "thread_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}})
    started = monotonic()
    rows = run_cases(cases, args.dimension, str(args.kernel.resolve()), args.out, args.case_seconds)
    summary = {"mode": args.mode, "completed_cells": len(rows), "planned_shard_cells": len(cases),
               "all_planned_cells": len(all_cases), "status_counts": dict(Counter(r["status"] for r in rows)),
               "seconds": monotonic() - started, "complete_shard": len(rows) == len(cases),
               "all_polynomials_completed": all(r["status"] in ("evaluated", "fixed_diagram_evaluated") for r in rows),
               "analytic_source_correspondence_established": False}
    dump(args.out / "summary.json", summary)
    print("MAP_SUMMARY " + json.dumps(summary), flush=True)
    if args.mode == "preflight" and not summary["all_polynomials_completed"]:
        raise SystemExit("preflight did not complete every required evaluation")


if __name__ == "__main__":
    main()
