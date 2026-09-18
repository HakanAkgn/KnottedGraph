#!/usr/bin/env python3
"""Certify original TPMS corridors or whole cells; freeze, archive, replay.

Run with the package installed, or PYTHONPATH=src. Always use a new output
folder. Unknown is never converted to inequivalent or assigned a phase label.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import hashlib
import inspect
import json
import math
from pathlib import Path
import platform
from time import perf_counter

import mpmath
import numpy as np
import sympy

from knotted_graph.core.field_isotopy import certify, tpms_problem, verify

BASE = "fab2673f665796c65a32ec42efdc1bf9b6c37678"
FAMILIES = ("schwarz_p_to_diamond", "gyroid_to_schwarz_p", "gyroid_to_diamond")


def write_json(path: Path, data: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def build_plan(mode: str, samples: int) -> list[dict]:
    lambdas, levels = np.linspace(0., 1., samples), np.linspace(0., .3, samples)
    if mode == "historical-nine":
        lambdas, levels = np.linspace(0., 1., 21), np.linspace(0., .3, 21)
        indices = ((0, 10, 19), (1, 12, 19), (1, 11, 19))
        return [{"id": f"{family}_lambda{j:03d}_c007", "family": family,
                 "lambda_index": j, "c_index": 7,
                 "lambda_bounds": [float(lambdas[j]), float(lambdas[j+1])],
                 "c_bounds": [float(levels[7]), float(levels[7])]}
                for family, columns in zip(FAMILIES, indices) for j in columns]
    return [{"id": f"{family}_lambda{j:03d}_c{k:03d}", "family": family,
             "lambda_index": j, "c_index": k,
             "lambda_bounds": [float(lambdas[j]), float(lambdas[j+1])],
             "c_bounds": [float(levels[k]), float(levels[k+1])]}
            for family in FAMILIES for k in range(samples-1) for j in range(samples-1)]


def run_case(payload: tuple) -> dict:
    row, directory, max_boxes, max_depth, radius = payload
    problem = tpms_problem(row["family"], row["lambda_bounds"], row["c_bounds"], radius=radius)
    certificate = certify(problem, max_boxes=max_boxes, max_depth=max_depth)
    raw = json.dumps(certificate, separators=(",", ":"), allow_nan=False).encode()
    archived = gzip.compress(raw, mtime=0)
    filename = "certificates/" + row["id"] + ".json.gz"
    (Path(directory)/filename).write_bytes(archived)
    replay_start = perf_counter()
    replay = (verify(problem, json.loads(gzip.decompress(archived)))
              if certificate["status"] == "certified"
              else {"valid": False, "reason": "unknown_is_not_a_certificate"})
    if certificate["status"] == "certified" and not replay["valid"]:
        raise RuntimeError(f"replay failed: {row['id']}: {replay}")
    return {**row, "status": certificate["status"], "reason": certificate["reason"],
            "checked_boxes": certificate["checked_boxes"], "counts": certificate["counts"],
            "solver_seconds": certificate["elapsed_seconds"],
            "replay_seconds": perf_counter()-replay_start, "replay_valid": replay["valid"],
            "certificate": filename, "certificate_sha256": hashlib.sha256(archived).hexdigest(),
            "problem_sha256": problem.fingerprint, "unresolved": certificate.get("unresolved")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--mode", choices=("historical-nine", "grid"), default="historical-nine")
    parser.add_argument("--samples", type=int, default=21)
    parser.add_argument("--jobs", type=int, choices=range(1, 5), default=4)
    parser.add_argument("--max-boxes", type=int, default=250000)
    parser.add_argument("--max-depth", type=int, default=72)
    args = parser.parse_args()
    if args.samples < 2 or args.max_boxes < 1 or args.max_depth < 1:
        parser.error("samples >= 2 and positive budgets are required")
    if args.out.exists():
        parser.error("output exists; use a new directory to retain earlier attempts")
    args.out.mkdir(parents=True)
    (args.out/"certificates").mkdir()
    cases = build_plan(args.mode, args.samples)
    radius = 0.72*(2.25*math.pi)
    plan = {"base_code_commit": BASE, "mode": args.mode, "samples_per_axis": args.samples,
            "planned_cases": len(cases), "jobs": args.jobs, "max_boxes": args.max_boxes,
            "max_depth": args.max_depth, "radius": radius, "radius_hex": radius.hex(),
            "module_sha256": hashlib.sha256(Path(inspect.getfile(certify)).read_bytes()).hexdigest(),
            "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "environment": {"python": platform.python_version(), "platform": platform.platform(),
                            "numpy": np.__version__, "sympy": sympy.__version__, "mpmath": mpmath.__version__},
            "scope": "continuous_analytic_solids_in_the_original_fixed_ball",
            "voxel_resolution": None, "display_filter": False,
            "selection": ("same nine pairs as audited embedding_contraction_pair_probe.json"
                          if args.mode == "historical-nine" else "every adjacent rectangle in the original parameter grid"),
            "cases": cases}
    write_json(args.out/"plan.json", plan)
    records = []
    started = perf_counter()
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_case, (row, str(args.out), args.max_boxes, args.max_depth, radius))
                   for row in cases]
        for future in as_completed(futures):
            row = future.result()
            records.append(row)
            print(f"{len(records):4d}/{len(cases)} {row['id']} {row['status']} "
                  f"boxes={row['checked_boxes']} time={row['solver_seconds']:.3f}", flush=True)
            write_json(args.out/"records.json", sorted(records, key=lambda r: r["id"]))
            write_json(args.out/"summary.json", {
                "complete": len(records) == len(cases), "completed_cases": len(records),
                "planned_cases": len(cases), "status_counts": dict(Counter(r["status"] for r in records)),
                "by_family": {f: dict(Counter(r["status"] for r in records if r["family"] == f)) for f in FAMILIES},
                "all_successes_replayed": all(r["replay_valid"] for r in records if r["status"] == "certified"),
                "wall_seconds": perf_counter()-started, "scope": plan["scope"],
                "distinct_topology_class_count_established": False, "old_graph_spines_validated": False,
                "plan_sha256": hashlib.sha256((args.out/"plan.json").read_bytes()).hexdigest()})
    print(json.dumps(json.loads((args.out/"summary.json").read_text()), indent=2))


if __name__ == "__main__":
    main()
