"""Execute public reconstruction scans and certificate-only TPMS comparisons.

These checks do not replace the complete historical Hamiltonian map. Every
result records its actual scope and preserves unsuccessful evaluations.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import gzip
from hashlib import sha256
import json
import multiprocessing as mp
from pathlib import Path
import platform
import subprocess
from time import monotonic


def write(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def public_case(case, dimension, directory):
    from knotted_graph.applications.nodal.deformation import NodalBlochPath, NodalPhaseScan
    from knotted_graph.applications.nodal.models import (
        hopf_link_bloch_vector, trefoil_bloch_vector, pq_torus_knot_bloch_vector,
    )
    started = monotonic()
    row = {**case, "dimension": dimension}
    try:
        path = (NodalBlochPath(hopf_link_bloch_vector, trefoil_bloch_vector)
                if case["family"] == "hopf_to_trefoil" else
                NodalBlochPath(trefoil_bloch_vector, lambda g: pq_torus_knot_bloch_vector(2, 5, g)))
        result = NodalPhaseScan(path, lambdas=[case["lambda"]], gammas=[case["level"]],
                               dimension=dimension, reconstruction="cubical").run()
        record = asdict(result.records[0])
        record["yamada"] = None if record["yamada"] is None else str(record["yamada"])
        row["record"] = record
        row["phase_grid"], row["phase_names"] = result.phase_grid()
        row["phase_grid"] = row["phase_grid"].tolist()
        row["transition_intervals"] = result.transition_intervals()
        row["status"] = "evaluated" if result.records[0].available else "unavailable"
    except Exception as exc:
        row.update({"status": "error", "error": f"{type(exc).__name__}: {exc}"})
    row["seconds"] = monotonic() - started
    row["complete_phase_map"] = False
    write(Path(directory) / "record.json", row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tpms", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dimension", type=int, default=120)
    parser.add_argument("--case-seconds", type=int, default=180)
    args = parser.parse_args()
    if args.out.exists() or args.dimension < 2 or args.case_seconds < 1:
        parser.error("use a new output directory and positive valid controls")
    args.out.mkdir(parents=True)
    from knotted_graph.core.continuation_atlas import ContinuationAtlas
    from knotted_graph.core.field_isotopy import tpms_problem
    plan = json.loads((args.tpms / "plan.json").read_text())
    records = json.loads((args.tpms / "records.json").read_text())
    atlases, queries = {}, []
    for row in records:
        if row["status"] != "certified":
            queries.append({"id": row["id"], "status": "unknown"})
            continue
        raw = (args.tpms / row["certificate"]).read_bytes()
        if sha256(raw).hexdigest() != row["certificate_sha256"]:
            raise RuntimeError("archived certificate hash mismatch")
        proof = json.loads(gzip.decompress(raw))
        problem = tpms_problem(row["family"], row["lambda_bounds"], row["c_bounds"], radius=plan["radius"])
        atlas = atlases.setdefault(row["family"], ContinuationAtlas())
        atlas.add(problem, proof)
        first = (row["lambda_bounds"][0], row["c_bounds"][0])
        last = (row["lambda_bounds"][1], row["c_bounds"][1])
        queries.append({"id": row["id"], **atlas.compare(first, last)})
    write(args.out / "tpms_equivalence_paths.json", queries)
    print("ATLAS_SUMMARY " + json.dumps(dict(Counter(r["status"] for r in queries))), flush=True)
    cases = [
        {"id": "hopf_low", "family": "hopf_to_trefoil", "lambda": 0., "level": .3},
        {"id": "hopf_trefoil_low", "family": "hopf_to_trefoil", "lambda": .5, "level": .3},
        {"id": "hopf_trefoil_high", "family": "hopf_to_trefoil", "lambda": .5, "level": .3 + 99 * 15 / 980},
        {"id": "cinquefoil_low", "family": "trefoil_to_cinquefoil", "lambda": 1., "level": .3},
    ]
    write(args.out / "plan.json", {"cases": cases, "dimension": args.dimension,
        "selection": "public-entry-point regression cases; not the full 7140-point map",
        "tpms_input_records_sha256": sha256((args.tpms / "records.json").read_bytes()).hexdigest(),
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "environment": {"python": platform.python_version(), "platform": platform.platform()},
        "source_hashes": {str(path): sha256(path.read_bytes()).hexdigest()
                          for path in Path("src/knotted_graph").rglob("*.py")}})
    context = mp.get_context("spawn")
    outputs = []
    for case in cases:
        directory = args.out / case["id"]
        directory.mkdir()
        worker = context.Process(target=public_case, args=(case, args.dimension, str(directory)))
        worker.start()
        worker.join(args.case_seconds)
        if worker.is_alive():
            worker.terminate()
            worker.join()
            row = {**case, "status": "time_budget", "unfinished_is_not_a_result": True}
        elif (directory / "record.json").is_file():
            row = json.loads((directory / "record.json").read_text())
        else:
            row = {**case, "status": "worker_error", "exitcode": worker.exitcode}
        outputs.append(row)
        print("PUBLIC_SCAN_CASE " + json.dumps(row, sort_keys=True), flush=True)
        write(args.out / "public_scan_records.json", outputs)
    summary = {"tpms_queries": dict(Counter(r["status"] for r in queries)),
               "public_scan_outcomes": dict(Counter(r["status"] for r in outputs)),
               "full_hamiltonian_map_rebuilt": False, "manuscript_rebuilt": False,
               "new_analytic_to_discrete_correspondence_claim": False,
               "all_word_formulas_changed": False, "existing_cavity_route_changed": False}
    write(args.out / "summary.json", summary)
    print("SCOPED_PUBLIC_SUMMARY " + json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
