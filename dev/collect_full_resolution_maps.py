#!/usr/bin/env python3
"""Collect exact full-grid runs and reject missing, duplicated or mixed records."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from run_full_resolution_maps import dump, plan


def collect(root: Path, out: Path, *, shards: int = 20):
    out.mkdir(parents=True, exist_ok=False)
    manifest = {}
    for mode in ("hamiltonian", "tpms"):
        expected = {r["id"]: r for r in plan(mode)}
        rows = {}
        commits = set()
        plans = []
        for shard in range(shards):
            directory = root / f"full-map-shard-{shard}" / mode
            spec = json.loads((directory / "plan.json").read_text())
            provenance = json.loads((directory / "provenance.json").read_text())
            summary = json.loads((directory / "summary.json").read_text())
            records = json.loads((directory / "records.json").read_text())
            if (spec["mode"] != mode or spec["shard"] != shard or spec["shards"] != shards
                    or spec["dimension"] != 120 or spec["display_filter"] or spec["abstract_fallback"]):
                raise ValueError("mixed map plan")
            intended = list(expected)[shard::shards]
            if ([r["id"] for r in spec["cases"]] != intended or not summary["complete_shard"]
                    or len(records) != len(intended)):
                raise ValueError("missing shard records")
            commits.add(provenance["commit"])
            plans.append({"shard": shard, "plan_sha256": hashlib.sha256((directory / "plan.json").read_bytes()).hexdigest(),
                          "records_sha256": hashlib.sha256((directory / "records.json").read_bytes()).hexdigest()})
            for record in records:
                key = record["id"]
                if key in rows or key not in intended:
                    raise ValueError("duplicate or off-grid record")
                if any(record.get(k) != v for k, v in expected[key].items()):
                    raise ValueError("parameter value or identity changed")
                case_dir = directory / key
                if record["status"] in ("evaluated", "fixed_diagram_evaluated"):
                    evidence = record["reconstruction"]
                    if not evidence["replay_valid"] or evidence["analytic_source_correspondence_certified"]:
                        raise ValueError("unsupported reconstruction claim")
                    archive = case_dir / "collapse.npz"
                    if hashlib.sha256(archive.read_bytes()).hexdigest() != evidence["archive_sha256"]:
                        raise ValueError("collapse archive hash mismatch")
                    if hashlib.sha256((case_dir / "graph.json.gz").read_bytes()).hexdigest() != record["graph_sha256"]:
                        raise ValueError("graph archive hash mismatch")
                    if record["yamada"] is None:
                        raise ValueError("successful record lacks a polynomial")
                elif record.get("yamada") is not None:
                    raise ValueError("unavailable calculation has a polynomial label")
                record["artifact_relative_path"] = str(case_dir.relative_to(root))
                rows[key] = record
        if set(rows) != set(expected) or len(commits) != 1:
            raise ValueError("incomplete or mixed-commit full map")
        ordered = [rows[k] for k in expected]
        dump(out / f"{mode}_records.json", ordered)
        status = Counter(r["status"] for r in ordered)
        manifest[mode] = {"requested": len(expected), "recorded": len(rows), "dimension": 120,
                          "status_counts": dict(status), "commit": next(iter(commits)),
                          "all_requested_points_executed": True,
                          "all_spatial_evaluations_completed": all(r["status"] == "evaluated" for r in ordered),
                          "by_family": {f: dict(Counter(r["status"] for r in ordered if r["family"] == f))
                                        for f in sorted({r["family"] for r in ordered})},
                          "source_phase_equivalence_classification_complete": False,
                          "analytic_to_voxel_correspondence_established": False, "shards": plans}
        print("COMPLETE_MAP " + json.dumps({mode: manifest[mode]}, sort_keys=True), flush=True)
    dump(out / "full_map_manifest.json", manifest)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=20)
    args = parser.parse_args()
    collect(args.root, args.out, shards=args.shards)
