"""Check file integrity and summarize committed paper records.

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
FIGURE_DATA = Path("User_guide/applications/results/figure4_suppfig10")


def inspect_figure_data(root: Path) -> dict:
    """Check supplied file integrity and agreement between saved CSV/JSON records."""
    directory = root / FIGURE_DATA
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    failures = []
    for entry in manifest["files"]:
        path = root / entry["path"]
        if not path.is_file():
            failures.append(f"Missing file: {entry['path']}")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            failures.append(f"SHA-256 differs: {entry['path']}")
    if failures:
        return {"ok": False, "failures": failures}

    def read_json(relative):
        return json.loads((directory / relative).read_text(encoding="utf-8"))

    mixed = read_json("audit/mixed_family/records.json")
    pure = read_json("audit/discovery/short_exact_records.json")["records"]
    long_words = read_json("audit/discovery_long_words/records.json")
    mixed_certificate = read_json("audit/mixed_family/certificate.json")
    pure_certificate = read_json("audit/discovery/certificate.json")
    long_certificate = read_json("audit/discovery_long_words/certificate.json")

    def indexed(rows, name):
        result = {row["word"]: row for row in rows}
        if len(result) != len(rows):
            failures.append(f"Duplicate word in {name}")
        return result

    mixed_by_word = indexed(mixed, "mixed JSON")
    pure_by_word = indexed(pure, "pure-braid JSON")

    def check_csv(name, reference, expected_words, extra_fields=()):
        with (directory / "data" / name).open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        by_word = indexed(rows, name)
        if set(by_word) != set(expected_words):
            failures.append(f"CSV word identifiers differ: {name}")
        for word, row in by_word.items():
            original = reference.get(word)
            if original is None:
                continue
            for key, value in row.items():
                if key in extra_fields:
                    continue
                expected = original.get(key)
                expected = "" if expected is None else str(expected)
                if value != expected:
                    failures.append(f"CSV/JSON differs: {name}, word={word!r}, {key}")
        return rows

    check_csv("figure4_mixed_family_459.csv", mixed_by_word, mixed_by_word)
    check_csv("figure4_pure_braid_324.csv", pure_by_word, pure_by_word)
    homogeneous = check_csv(
        "figure4_homogeneous_m1_m2.csv", mixed_by_word,
        ("L", "LL", "D", "DD", "B", "BB"), ("family", "figure4_panels"),
    )
    panels = {"L": "b", "LL": "c", "D": "e", "DD": "f", "B": "h", "BB": "i"}
    for row in homogeneous:
        if row["figure4_panels"] != panels.get(row["word"]):
            failures.append(f"Homogeneous panel mapping differs: {row['word']}")
    check_csv(
        "figure4_order_sensitive_AAB_ABA_BAA.csv", pure_by_word,
        ("AAB", "ABA", "BAA"),
    )
    aaaba = read_json("data/suppfig10_AAABA.json")
    if aaaba["direct_record"] != pure_by_word.get("AAABA"):
        failures.append("Supplementary Figure 10 differs from its short-word record")
    if aaaba["is_independent_holdout"] is not False:
        failures.append("AAABA metadata must identify its role as a worked example")
    if len(mixed) != mixed_certificate["completed_cases"]:
        failures.append("Mixed-family count differs from its certificate")
    if len(pure) != pure_certificate["short_exact_words"]:
        failures.append("Pure-braid count differs from its certificate")
    if len(long_words) != long_certificate["completed_cases"]:
        failures.append("Completed long-word count differs from its certificate")
    for row in mixed:
        if (json.loads(row["actual_coefficients_json"])
                != json.loads(row["predicted_coefficients_json"])
                or row["master_formula_pass"] is not True):
            failures.append(f"Stored mixed-family coefficients disagree: {row['word']}")
    for index, row in enumerate(long_words, 1):
        prediction = read_json(f"audit/discovery_long_words/prediction_{index:02d}.json")
        if (prediction["word"] != row["word"]
                or prediction["raw_coefficients"]
                != json.loads(row["actual"]["raw_coefficients_json"])
                or row["raw_coefficient_identity"] is not True):
            failures.append(f"Stored long-word coefficients disagree: case {index}")

    return {
        "ok": not failures,
        "failures": failures,
        "verified_files": len(manifest["files"]),
        "mixed_family_cases": len(mixed),
        "short_pure_braid_cases": len(pure),
        "homogeneous_examples": len(homogeneous),
        "supplementary_figure_10_word": aaaba["word"],
        "completed_long_words": len(long_words),
        "planned_long_words": long_certificate["planned_cases"],
        "record_kind": manifest["record_kind"],
    }


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

    figure_data = inspect_figure_data(root)
    if not figure_data["ok"]:
        return figure_data

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
        "formula_discovery": figure_data,
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
        "scope": "Checks saved-record integrity, case alignment, and CSV/JSON agreement.",
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
        f = report["formula_discovery"]
        print(f"Figure 4 / S10: {f['verified_files']} files verified; "
              f"{f['mixed_family_cases']} mixed-family and "
              f"{f['short_pure_braid_cases']} short pure-braid records; "
              f"{f['completed_long_words']} completed long-word comparisons.")
        print("Figure 4 CSV exports and S10 worked example match their source JSON records.")
        print(report["scope"])
        for note in report["availability_notes"]:
            print(f"Availability: {note}")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
