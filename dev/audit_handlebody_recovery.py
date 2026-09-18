#!/usr/bin/env python3
"""Recheck exact archived embeddings without overwriting benchmark artifacts.

Use ``--case random_4920`` for the recorded mismatch, or ``--case missing``
for every row without a completed Yamada comparison. Run from the intended
checkout with ``uv run --no-sync python /path/to/this/script.py ...``.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time
import warnings

import networkx as nx
import numpy as np
import sympy as sp
from knotted_graph.projection import compute_yamada_polynomial
from knotted_graph.projection.pd_code import PDCode


def payload_hash(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def graph_from_payload(payload):
    graph = nx.MultiGraph()
    for node in payload["nodes"]:
        graph.add_node(node["id"], pos=np.asarray(node["pos"], dtype=float))
    for edge in payload["edges"]:
        graph.add_edge(edge["u"], edge["v"], key=edge["key"], pts=np.asarray(edge["pts"], dtype=float))
    return graph


def read_archive(path):
    path = Path(path)
    if path.is_dir() or path.name == "manifest.json":
        from shard_handlebody_inputs import load_manifest
        return load_manifest(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def evaluate(payload, angles=None):
    start = time.perf_counter()
    try:
        with warnings.catch_warnings(record=True) as caught:
            result = compute_yamada_polynomial(graph_from_payload(payload), sp.Symbol("A"), rotation_angles=angles, num_rotation_samples=6, normalize=True, n_jobs=1, return_result=True)
        return {"polynomial": str(result.polynomial), "crossings": result.projection.num_crossings,
                "angles": result.projection.rotation_angles, "warnings": [str(item.message) for item in caught],
                "seconds": time.perf_counter() - start}
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}", "seconds": time.perf_counter() - start}


def polynomial_equal(first, second):
    if not first.get("polynomial") or not second.get("polynomial"):
        return None
    return sp.expand(sp.sympify(first["polynomial"]) - sp.sympify(second["polynomial"])) == 0


def seed_matches_row(row, seed):
    """Check recorded generation/certification fields before pairing by name.

    These checks can disprove identity; metadata agreement is not an exact
    geometry certificate, so a fresh recovery records the seed payload hash.
    """
    payload = seed["embedding"]
    checks = {"nodes": len(payload["nodes"]) == int(row["seed_nodes"]),
              "edges": len(payload["edges"]) == int(row["seed_edges"]),
              "wl_hash": not row.get("wl_hash") or row["wl_hash"] == seed.get("wl_hash")}
    if row.get("generation_seed"):
        checks["generation_seed"] = int(float(row["generation_seed"])) == seed.get("generation_seed")
    for key in ("radius_world", "spacing_world", "occupied_voxels"):
        value = seed.get("certification", {}).get(key)
        if row.get(key) and value is not None:
            checks[key] = bool(np.isclose(float(row[key]), float(value), rtol=1e-12, atol=1e-12))
    return checks


def audit_case(row, seed, recovered, *, views=False):
    original = seed["embedding"]
    alignment = seed_matches_row(row, seed)
    extracted = recovered["recovered_graph"]
    digest = payload_hash(extracted)
    hashes = [row.get("recovered_graph_sha256"), row.get("recovered_yamada_graph_sha256"), recovered.get("recovered_graph_sha256")]
    expected_hashes = [value for value in hashes if value]
    result = {"case": row["case"], "index": row["index"], "archived_yamada_match": row["yamada_match"],
              "archived_overall_pass": row["overall_pass"], "recovered_payload_sha256": digest,
              "csv_archive_graph_hashes_agree": bool(expected_hashes) and all(value == digest for value in expected_hashes),
              "seed_payload_sha256": payload_hash(original), "seed_csv_alignment_checks": alignment,
              "seed": evaluate(original) if all(alignment.values()) else {"error": "Seed archive does not match CSV generation/certification metadata; comparison withheld."},
              "recovered": evaluate(extracted) if expected_hashes and all(value == digest for value in expected_hashes) else {"error": "Recovered graph archive does not match CSV/archive hashes; comparison withheld."}}
    result["yamada_match"] = polynomial_equal(result["seed"], result["recovered"])
    if views and all(alignment.values()):
        result["explicit_views"] = [
            {"angles": angles, "seed": evaluate(original, angles), "recovered": evaluate(extracted, angles)}
            for angles in ((0., 0., 0.), (17., 31., 47.), (31., 67., 13.), (-32.46117974981078, 23.556464309101237, 0.))
        ]
    return result


def notebook_recovery_functions(notebook_path):
    """Load only the notebook's named pure functions, never its execution cells."""
    import ast
    from scipy.ndimage import generate_binary_structure, label as ndi_label
    from skimage.measure import euler_number
    from knotted_graph.core import contract_short_edges, ensure_embedding, remove_leaf_nodes, simplify_edges, smooth_edges
    from knotted_graph.extraction import skeletonize_volume, topology_aware_skeleton_image_to_graph
    wanted = {"graph_cycle_rank", "graph_stats", "all_edge_points", "voxelize_regular_neighborhood",
              "volume_betti", "expected_handlebody_betti", "graph_voxel_to_world", "recover_graph_from_volume",
              "_json_node_id", "graph_to_serializable_exact"}
    ns = dict(globals(), **locals())
    ns.update(MAX_JUNCTION_DEGREE=3, ADAPTIVE_MAX_HOPS=4, ANOMALY_RATIO=.15,
              CONTRACT_SHORT_EDGE_VOX=1.75, SMOOTH_EPSILON_VOX=1.5)
    notebook = json.loads(Path(notebook_path).read_text())
    definitions = []
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        tree = ast.parse("".join(cell["source"]))
        definitions.extend(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in wanted)
    if {node.name for node in definitions} != wanted:
        raise ValueError("Notebook recovery functions changed; inspect and update the audited whitelist.")
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(notebook_path), "exec"), ns)
    return ns


def rebuild_case(row, seed, notebook_path):
    start = time.perf_counter()
    ns = notebook_recovery_functions(notebook_path)
    original = seed["embedding"]
    graph = graph_from_payload(original)
    grid = int(seed["generation_config"]["GRID_SIZE"])
    volume, origin, spacing, _ = ns["voxelize_regular_neighborhood"](graph, seed["certification"]["radius_world"], grid)
    betti = ns["volume_betti"](volume)
    expected = ns["expected_handlebody_betti"](graph)
    recovered = ns["recover_graph_from_volume"](volume, origin=origin, spacing=spacing)
    payload = ns["graph_to_serializable_exact"](recovered["graph"])
    seed_stats = ns["graph_stats"](graph)
    recovered_stats = ns["graph_stats"](recovered["graph"])
    result = {"mode": "fresh_recovery", "case": seed["case"], "index": row["index"],
              "seed_payload_sha256": payload_hash(original), "recovered_payload_sha256": payload_hash(payload),
              "seed_csv_alignment_checks": seed_matches_row(row, seed), "seed": evaluate(original), "recovered": evaluate(payload),
              "seed_stats": seed_stats, "recovered_stats": recovered_stats,
              "volume_betti": betti, "expected_betti": expected, "occupied_voxels": int(volume.sum()),
              "recovered_graph": payload, "recovery_seconds": recovered["recovery_total_seconds"],
              "skeleton_voxels": int(np.count_nonzero(recovered["skeleton"]))}
    result["yamada_match"] = polynomial_equal(result["seed"], result["recovered"])
    result["yamada_status"] = "unavailable" if result["yamada_match"] is None else "match" if result["yamada_match"] else "mismatch"
    result["input_topology_match"] = betti == expected
    result["overall_pass"] = bool(result["input_topology_match"] and recovered_stats["connected"] and recovered_stats["subcubic"]
                                  and recovered_stats["cycle_rank"] == seed_stats["cycle_rank"] and result["yamada_match"] is True)
    result["seconds"] = time.perf_counter() - start
    return result


def replay_job(job):
    row, seed, recovered, views = job
    return audit_case(row, seed, recovered, views=views)


def ordered_parallel_jobs(jobs, workers):
    from collections import deque
    from concurrent.futures import ProcessPoolExecutor
    iterator = iter(jobs)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = deque()
        for _ in range(workers * 2):
            job = next(iterator, None)
            if job is not None:
                pending.append(pool.submit(replay_job, job))
        while pending:
            yield pending.popleft().result()
            job = next(iterator, None)
            if job is not None:
                pending.append(pool.submit(replay_job, job))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-archive", type=Path, required=True)
    parser.add_argument("--recovered-archive", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--case", default="random_4920", help="Case name, 'missing', or 'all'")
    parser.add_argument("--output", type=Path, required=True, help="New JSONL audit file; existing files are refused")
    parser.add_argument("--workers", type=int, choices=range(1, 5), default=1, help="At most four processes for archived polynomial replay")
    parser.add_argument("--reuse", type=Path, help="Reuse compatible replay JSONL only after input/code/version/hash checks")
    parser.add_argument("--rebuild", action="store_true", help="Fresh volume/extraction from each immutable archived seed, storing each new exact graph")
    parser.add_argument("--notebook", type=Path, help="Notebook04 source for the audited recovery function whitelist (required with --rebuild)")
    parser.add_argument("--limit", type=int, default=None, help="Evaluate at most this many selected cases")
    parser.add_argument("--views", action="store_true", help="Also evaluate four explicit projections")
    args = parser.parse_args()
    if args.rebuild and (args.reuse or args.workers != 1):
        parser.error("Fresh recovery is serial and cannot reuse archived-pair results; omit --reuse and use --workers 1")
    if args.rebuild and args.notebook is None:
        parser.error("--rebuild requires --notebook")
    seeds = read_archive(args.seed_archive)
    recovered = read_archive(args.recovered_archive)
    seed_map = {row["case"]: row for row in seeds["cases"]}
    recovered_map = {row["case"]: row for row in recovered["cases"]}
    with args.csv.open(newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    selected = [row for row in all_rows if args.case == "all" or (args.case == "missing" and not row["yamada_match"]) or row["case"] == args.case]
    selected = selected[:args.limit]
    if not selected:
        parser.error(f"No cases selected by {args.case!r}")
    import inspect
    code_path = Path(inspect.getsourcefile(PDCode))
    versions = {name: importlib.metadata.version(name) for name in ("shapely", "numpy", "sympy", "networkx")}
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=code_path.parent, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    code_hashes = {str(path.relative_to(code_path.parents[2])): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(code_path.parents[1].rglob("*.py"))}
    historical = Counter(row["yamada_match"] or "missing" for row in all_rows)
    provenance = {"type": "provenance", "archive_recovered_metadata": {key: value for key, value in recovered.items() if key != "cases"},
                  "archive_seed_config": seeds.get("config"), "archive_csv_rows": len(all_rows), "archive_comparison_counts": dict(historical),
                  "selected_cases": len(selected), "code_commit": commit, "projection_source_sha256": hashlib.sha256(code_path.read_bytes()).hexdigest(),
                  "versions": versions, "source_hashes": code_hashes, "files": {str(path.resolve()): hashlib.sha256((path / "manifest.json" if path.is_dir() else path).read_bytes()).hexdigest() for path in (args.seed_archive, args.recovered_archive, args.csv)}}
    provenance["mode"] = "fresh_recovery" if args.rebuild else "archived_graph_projection_recheck"
    if args.notebook:
        provenance["notebook_sha256"] = hashlib.sha256(args.notebook.read_bytes()).hexdigest()
    provenance["relevant_source_hashes"] = {key: value for key, value in code_hashes.items() if "/projection/" in key or "/invariants/" in key}
    provenance["workers"] = args.workers
    provenance["runner_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cached = {}
    if args.reuse:
        previous = [json.loads(line) for line in args.reuse.read_text().splitlines()]
        old_provenance = previous[0]
        old_sources = {key: value for key, value in old_provenance["source_hashes"].items() if "/projection/" in key or "/invariants/" in key}
        if old_sources != provenance["relevant_source_hashes"] or old_provenance["versions"] != versions or old_provenance["files"] != provenance["files"]:
            parser.error("Reuse rejected: projection/invariant source, dependencies, or input archive hashes differ")
        for result in previous[1:]:
            name = result.get("case")
            if name in seed_map and name in recovered_map and result.get("seed_payload_sha256") == payload_hash(seed_map[name]["embedding"]) and result.get("recovered_payload_sha256") == payload_hash(recovered_map[name]["recovered_graph"]):
                result["reused_from"] = str(args.reuse.resolve())
                cached[name] = result
    provenance["reused_cases"] = len(cached)
    new_rows = [row for row in selected if row["case"] not in cached]
    if args.rebuild:
        results_iter = (rebuild_case(row, seed_map[row["case"]], args.notebook) for row in new_rows)
    else:
        jobs = ((row, seed_map[row["case"]], recovered_map[row["case"]], args.views) for row in new_rows)
        results_iter = map(replay_job, jobs) if args.workers == 1 else ordered_parallel_jobs(jobs, args.workers)
    counts = Counter()
    with args.output.open("x") as handle:
        handle.write(json.dumps(provenance) + "\n"); handle.flush()
        for i, row in enumerate(selected):
            result = cached[row["case"]] if row["case"] in cached else next(results_iter)
            counts[str(result["yamada_match"])] += 1
            handle.write(json.dumps(result) + "\n"); handle.flush()
            print(json.dumps({"done": i + 1, "case": row["case"], "match": result["yamada_match"], "counts": counts}), flush=True)
        handle.write(json.dumps({"type": "summary", "counts": counts}) + "\n")


if __name__ == "__main__":
    main()
