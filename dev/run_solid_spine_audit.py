#!/usr/bin/env python3
"""Bounded actual-field tests of witnessed voxel-to-spine reconstruction.

These representative checks do not replace either manuscript's full map.
Original files, existing cavity processing and all-word formulas are untouched.
Every executed case stores its exact source mask, grid, collapse witness and
(if reached) the complete graph and audited spatial polynomial. A finite-box
result is never called a periodic Brillouin-zone result.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import gzip
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import platform
import subprocess
from time import monotonic

import numpy as np


def dump(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def sample(case, dimension):
    if case["source"] == "tpms":
        from knotted_graph.applications.phase_map_examples._tpms import implicit_fields
        fields = implicit_fields()
        endpoints = {
            "schwarz_p_to_diamond": ("schwarz_p", "diamond"),
            "gyroid_to_schwarz_p": ("gyroid", "schwarz_p"),
            "gyroid_to_diamond": ("gyroid", "diamond"),
        }
        axes = (np.linspace(-2.25 * np.pi, 2.25 * np.pi, dimension),) * 3
        x, y, z = np.meshgrid(*axes, indexing="ij")
        a, b = endpoints[case["family"]]
        field = (1 - case["lambda"]) * fields[a].function(x, y, z)
        field += case["lambda"] * fields[b].function(x, y, z)
        radius = .72 * (2.25 * np.pi)
        mask = (field <= case["level"]) & (x*x + y*y + z*z <= radius*radius)
        return mask, axes, {"domain": "original_spherical_coupon",
                            "radius_hex": float(radius).hex(), "periodic_identification": False}
    from knotted_graph.applications.nodal.deformation import NodalBlochPath
    from knotted_graph.applications.nodal.models import (
        hopf_link_bloch_vector, trefoil_bloch_vector, pq_torus_knot_bloch_vector,
    )
    from knotted_graph.applications.nodal.skeleton import NodalSkeleton
    if case["family"] == "hopf_to_trefoil":
        path = NodalBlochPath(hopf_link_bloch_vector, trefoil_bloch_vector)
    else:
        path = NodalBlochPath(trefoil_bloch_vector, lambda g: pq_torus_knot_bloch_vector(2, 5, g))
    model = NodalSkeleton(char=path.at(case["level"], case["lambda"]), dimension=dimension)
    return np.asarray(model._interior_mask, dtype=bool), (model.kx_vals, model.ky_vals, model.kz_vals), {
        "domain": "closed_finite_box_with_dual_cells_clipped_at_sample_endpoints",
        "periodic_identification": False, "mask": "unchanged NodalSkeleton._interior_mask",
        "visual_axis_scale_applied_to_topology": False,
    }


def evaluate(case, dimension, directory):
    started = monotonic()
    directory = Path(directory)
    row = {**case, "dimension": dimension, "status": "error", "error": None}
    try:
        from knotted_graph.applications.phase_maps import volume_topology, _compute_yamada_audited
        from knotted_graph.extraction.cubical_spine import (
            cubical_retract, verify_cubical_retract, sample_grid_breaks, certificate_json,
        )
        import sympy as sp
        mask, axes, domain = sample(case, dimension)
        row["domain"] = domain
        topology = volume_topology(mask)
        row["volume"] = asdict(topology)
        np.savez_compressed(directory / "source.npz", mask=mask, x=axes[0], y=axes[1], z=axes[2])
        row["source_npz_sha256"] = hashlib.sha256((directory / "source.npz").read_bytes()).hexdigest()
        if topology.enclosed_voids:
            row["status"] = "not_revised_existing_cavity_path"
            return
        if not topology.interior_voxels:
            row["status"] = "empty"
            return
        gridlines = tuple(sample_grid_breaks(a) for a in axes)
        result = cubical_retract(mask, gridlines=gridlines)
        raw = certificate_json(result.certificate).encode()
        archived = gzip.compress(raw, mtime=0)
        (directory / "collapse.json.gz").write_bytes(archived)
        row["certificate_sha256"] = hashlib.sha256(archived).hexdigest()
        replay = verify_cubical_retract(mask, json.loads(gzip.decompress(archived)), gridlines=gridlines)
        row["replay"] = replay
        if not replay["valid"]:
            raise RuntimeError(f"collapse replay failed: {replay}")
        row["status"] = result.certificate["kind"]
        row["initial_cell_counts"] = result.certificate["initial_cell_counts"]
        row["terminal_cell_counts"] = result.certificate["terminal_cell_counts"]
        row["collapses"] = len(result.certificate["collapses"])
        if result.graph is None:
            return
        graph = result.graph
        summary = result.certificate["graph_summary"]
        row["graph"] = summary
        if (summary["components"] != topology.connected_components
                or summary["cycle_rank"] != topology.handle_rank):
            raise RuntimeError("independent digital Betti check disagrees with certified graph")
        ids = {node: i for i, node in enumerate(graph)}
        dump(directory / "graph.json", {
            "nodes": [{"id": ids[n], "pos": list(d["pos"])} for n, d in graph.nodes(data=True)],
            "edges": [{"u": ids[u], "v": ids[v], "key": k, "pts": d["pts"].tolist()}
                      for u, v, k, d in graph.edges(keys=True, data=True)],
            "scope": "embedded graph retract of specified voxel complex",
        })
        # Expensive or failed spatial calculations are never replaced by an abstract value.
        row["polynomial"] = {"status": "above_explicit_18_edge_limit", "value": None}
        if graph.number_of_edges() <= 18:
            try:
                value, audit = _compute_yamada_audited(graph, sp.Symbol("Y"), {"normalize": True, "n_jobs": 1})
                row["polynomial"] = {"status": "evaluated", "value": str(value), "audit": audit}
            except Exception as exc:
                row["polynomial"] = {"status": "unavailable", "value": None,
                                     "error": f"{type(exc).__name__}: {exc}"}
    except Exception as exc:
        row["status"] = "error"
        row["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        row["seconds"] = monotonic() - started
        row["analytic_source_correspondence_certified"] = False
        dump(directory / "record.json", row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dimension", type=int, default=64)
    parser.add_argument("--seconds", type=int, default=420)
    parser.add_argument("--case-seconds", type=int, default=65)
    args = parser.parse_args()
    if args.out.exists() or args.dimension < 2 or min(args.seconds, args.case_seconds) < 1:
        parser.error("use a new output directory, dimension >= 2 and positive time budgets")
    args.out.mkdir(parents=True)
    cases = [
        {"id": "P_to_D_start", "source": "tpms", "family": "schwarz_p_to_diamond", "lambda": 0., "level": .105},
        {"id": "Hopf_start_low", "source": "hamiltonian", "family": "hopf_to_trefoil", "lambda": 0., "level": .3},
        {"id": "G_to_P_middle", "source": "tpms", "family": "gyroid_to_schwarz_p", "lambda": .5, "level": .105},
        {"id": "Hopf_trefoil_middle_high", "source": "hamiltonian", "family": "hopf_to_trefoil", "lambda": .5, "level": .3 + 99 * 15 / 980},
        {"id": "G_to_D_start", "source": "tpms", "family": "gyroid_to_diamond", "lambda": 0., "level": .105},
        {"id": "Cinquefoil_endpoint", "source": "hamiltonian", "family": "trefoil_to_cinquefoil", "lambda": 1., "level": .3},
    ]
    dump(args.out / "plan.json", {"cases": cases, "dimension": args.dimension,
        "selection": "six explicit source-field regression cases, not a phase-map replacement",
        "full_7140_scan_executed": False, "old_figure_replaced": False,
        "original_hamiltonian_plan": {"lambda_array": [i / 59 for i in range(60)],
            "energy_candidates": [.3 + 99 * i / 980 for i in range(50)],
            "family_energy_prefix_lengths": [16, 24, 16, 24, 39], "dimension": 120}})
    dump(args.out / "provenance.json", {"git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "python": platform.python_version(), "platform": platform.platform(),
        "source_hashes": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in Path("src/knotted_graph").rglob("*.py")},
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    started = monotonic()
    records = []
    context = mp.get_context("spawn")
    for case in cases:
        remaining = args.seconds - (monotonic() - started)
        if remaining <= 0:
            records.append({**case, "status": "not_executed_global_budget"})
            continue
        directory = args.out / case["id"]
        directory.mkdir()
        child = context.Process(target=evaluate, args=(case, args.dimension, str(directory)))
        child.start()
        child.join(min(args.case_seconds, remaining))
        if child.is_alive():
            child.terminate()
            child.join()
            row = {**case, "status": "case_time_budget", "partial_files_are_not_results": True}
            dump(directory / "record.json", row)
        elif (directory / "record.json").exists():
            row = json.loads((directory / "record.json").read_text())
        else:
            row = {**case, "status": "worker_error", "exitcode": child.exitcode}
        records.append(row)
        print("SPINE_CASE " + json.dumps(row, sort_keys=True), flush=True)
        dump(args.out / "records.json", records)
    summary = {"dimension": args.dimension, "planned_cases": len(cases),
        "status_counts": dict(Counter(r["status"] for r in records)),
        "seconds": monotonic() - started, "full_phase_map_rebuilt": False,
        "analytical_solid_to_discrete_solid_correspondence_established": False,
        "cavity_path_changed": False, "all_word_formulas_changed": False}
    dump(args.out / "summary.json", summary)
    print("SPINE_SUMMARY " + json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
