#!/usr/bin/env python3
"""Refine every unresolved TPMS rectangle, retaining a complete cover ledger.

Each accepted leaf has both full interval replays. Rejected parents are never
silently counted as accepted: their exact dyadic subdivisions are retained.
The remainder encloses all possible discriminant events but is not asserted
to consist entirely of transitions. Disconnected accepted components are not
asserted inequivalent. This operates on the original analytic source solids,
not on newly extracted graphs or their Yamada values.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
from time import monotonic

from certify_tpms_continuation import build_plan, write_json
from knotted_graph.core.field_isotopy import certify, tpms_problem, verify
from knotted_graph.core.field_isotopy_rational import verify_rational


def subdivide(row):
    lb, cb = row["lambda_bounds"], row["c_bounds"]
    lm, cm = lb[0] / 2 + lb[1] / 2, cb[0] / 2 + cb[1] / 2
    if not lb[0] < lm < lb[1] or not cb[0] < cm < cb[1]:
        raise ValueError("parameter subdivision exhausted floating precision")
    return [{**row, "id": row["id"] + f"_{j}{k}", "depth": row["depth"] + 1,
             "lambda_bounds": [lb[0], lm] if j == 0 else [lm, lb[1]],
             "c_bounds": [cb[0], cm] if k == 0 else [cm, cb[1]],
             "tile_x": 2 * row["tile_x"] + j, "tile_y": 2 * row["tile_y"] + k}
            for k in range(2) for j in range(2)]


def evaluate(row, out, max_boxes):
    problem = tpms_problem(row["family"], row["lambda_bounds"], row["c_bounds"])
    certificate = certify(problem, max_boxes=max_boxes, max_depth=72)
    raw = gzip.compress(json.dumps(certificate, separators=(",", ":"), allow_nan=False).encode(), mtime=0)
    name = "certificates/" + row["id"] + ".json.gz"
    (out / name).write_bytes(raw)
    primary = rational = {"valid": False, "reason": "unresolved"}
    if certificate["status"] == "certified":
        loaded = json.loads(gzip.decompress(raw))
        primary = verify(problem, loaded)
        rational = verify_rational(problem, loaded)
        if not primary["valid"] or not rational["valid"]:
            raise RuntimeError("primary/rational full cover replay disagreement")
    return {**row, "status": certificate["status"], "reason": certificate["reason"],
            "checked_boxes": certificate["checked_boxes"], "certificate": name,
            "certificate_sha256": hashlib.sha256(raw).hexdigest(), "problem_sha256": problem.fingerprint,
            "primary_replay": primary, "rational_replay": rational,
            "unresolved_spatial_box": certificate.get("unresolved")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=20)
    parser.add_argument("--depth", type=int, default=2)
    args = parser.parse_args()
    if args.out.exists() or not 0 <= args.shard < args.shards or not 0 <= args.depth <= 4:
        parser.error("new output directory and valid shard/depth required")
    (args.out / "certificates").mkdir(parents=True)
    roots = [{**r, "depth": 0, "tile_x": r["lambda_index"], "tile_y": r["c_index"]}
             for r in build_plan("grid", 21)][args.shard::args.shards]
    write_json(args.out / "plan.json", {"shard": args.shard, "shards": args.shards,
               "max_parameter_depth": args.depth, "root_count_total": 1200, "roots": roots,
               "base_max_boxes": 250000, "refined_max_boxes": 100000,
               "source": "original analytic TPMS fields in fixed ball", "display_filter": False,
               "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()})
    started = monotonic()
    records = []
    completed_roots = 0
    for root in roots:
        stack = [root]
        while stack:
            row = stack.pop()
            result = evaluate(row, args.out, 250000 if row["depth"] == 0 else 100000)
            result["children"] = []
            if result["status"] != "certified" and row["depth"] < args.depth:
                children = subdivide(row)
                result["children"] = [r["id"] for r in children]
                stack.extend(reversed(children))
            records.append(result)
        completed_roots += 1
        write_json(args.out / "records.json", records)
        leaves = [r for r in records if not r["children"]]
        accepted_area = sum(4 ** (-r["depth"]) for r in leaves if r["status"] == "certified")
        unknown_area = sum(4 ** (-r["depth"]) for r in leaves if r["status"] != "certified")
        if accepted_area + unknown_area != completed_roots:
            raise RuntimeError("adaptive leaf partition lost or double-counted parameter area")
        summary = {"completed_roots": completed_roots, "planned_roots": len(roots),
                   "complete_shard": completed_roots == len(roots), "leaf_counts": dict(Counter(r["status"] for r in leaves)),
                   "certified_root_cell_area": accepted_area, "unresolved_root_cell_area": unknown_area,
                   "seconds": monotonic() - started, "every_success_doubly_replayed": True,
                   "unresolved_means_inequivalent": False, "distinct_isotopy_classes_proved": False}
        write_json(args.out / "summary.json", summary)
        if completed_roots % 5 == 0 or completed_roots == len(roots):
            print("ADAPTIVE_TPMS " + json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
