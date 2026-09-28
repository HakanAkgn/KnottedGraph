"""Check and summarize committed paper records; never run a scientific calculation.

Usage from a source checkout: uv run --no-project python scripts/inspect_paper_data.py
Only the Python standard library is required. Exit 1 means a missing/changed
file, malformed record, or mismatch between the saved case identifiers.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path("User_guide/benchmarks/results/manifest.json")


def inspect(root: Path) -> dict:
    manifest = json.loads((root / MANIFEST).read_text(encoding="utf-8"))
    tables = {}
    failures = []
    for entry in manifest["files"]:
        path = root / entry["path"]
        if not path.is_file():
            failures.append(f"Missing file: {entry['path']}")
            continue
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            failures.append(f"SHA-256 differs: {entry['path']}")
        with path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            rows = list(reader)
            if reader.fieldnames != entry["columns"]:
                failures.append(f"CSV columns differ: {entry['path']}")
        if len(rows) != entry["rows"]:
            failures.append(f"CSV row count differs: {entry['path']}")
        tables[entry["id"]] = rows
        key = entry.get("case_key")
        if key and (any(not row.get(key) for row in rows)
                    or len({row.get(key) for row in rows}) != len(rows)):
            failures.append(f"Missing/duplicate {key}: {entry['path']}")

    aligned = ["preservation", "timings", "timing_plot", "input_manifest",
               "recovered_manifest", "checkpoint"]
    if all(name in tables for name in aligned):
        cases = {row["case"] for row in tables["preservation"]}
        for name in aligned[1:]:
            if {row["case"] for row in tables[name]} != cases:
                failures.append(f"Handlebody case identifiers differ: {name}")

    if failures:
        return {"ok": False, "failures": failures}

    scaling = tables["scaling"]
    preservation = tables["preservation"]
    plot = tables["timing_plot"]

    def counts(rows, field):
        return dict(sorted(Counter(row[field] for row in rows).items()))

    def measurements(crossings):
        return [{key: row[key] for key in (
            "case_key", "crossings", "knottedgraph_s", "topoly_s",
            "knottedgraph_result", "topoly_status", "topoly_result",
            "topoly_over_knottedgraph",
        )} for row in scaling if int(row["crossings"]) == crossings]

    return {
        "ok": True,
        "record_snapshot_commit": manifest["record_snapshot_commit"],
        "verified_files": len(tables),
        "scaling": {
            "rows": len(scaling),
            "max_crossings": max(int(row["crossings"]) for row in scaling),
            "knottedgraph_result": counts(scaling, "knottedgraph_result"),
            "topoly_status": counts(scaling, "topoly_status"),
            "topoly_result": counts(scaling, "topoly_result"),
            "seven_crossings": measurements(7),
            "five_hundred_crossings": measurements(500),
        },
        "handlebody": {
            "rows": len(preservation),
            "yamada_match": counts(preservation, "yamada_match"),
            "primary_topology_pass": counts(preservation, "primary_topology_pass"),
            "overall_pass": counts(preservation, "overall_pass"),
            "abstract_isomorphic": counts(preservation, "abstract_isomorphic"),
            "topoly_status": counts(plot, "topoly_status"),
        },
        "scope": "Checks saved records and their integrity, not independent polynomial recomputation.",
        "availability_notes": manifest["availability_notes"],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="Repository checkout to inspect")
    parser.add_argument("--json", action="store_true", help="Print a machine-readable report")
    args = parser.parse_args(argv)
    try:
        report = inspect(args.root.resolve())
    except (OSError, ValueError, KeyError, csv.Error) as error:
        report = {"ok": False, "failures": [str(error)]}
    if args.json:
        print(json.dumps(report, indent=2))
    elif not report["ok"]:
        print("Saved-record check failed:", file=sys.stderr)
        print("\n".join(report["failures"]), file=sys.stderr)
    else:
        print(f"PASS: {report['verified_files']} CSV files match the recorded hashes, columns and row counts.")
        s, h = report["scaling"], report["handlebody"]
        print(f"Scaling: {s['rows']} cases, maximum {s['max_crossings']} crossings; KnottedGraph {s['knottedgraph_result']}.")
        print(f"Scaling Topoly execution status: {s['topoly_status']}")
        print(f"Handlebody: {h['rows']} aligned cases; Yamada match {h['yamada_match']}; overall pass {h['overall_pass']}.")
        print(f"Handlebody abstract isomorphism (separate diagnostic): {h['abstract_isomorphic']}")
        print(f"Handlebody Topoly status: {h['topoly_status']}")
        print(report["scope"])
        for note in report["availability_notes"]:
            print(f"Availability: {note}")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
