#!/usr/bin/env python3
"""Build compact browsing tables from the copied retrospective audit JSON.

Run from any directory with: python3 REGENERATE_EXTRACTS.py
No Yamada calculations are performed here. The original JSON files remain the
authoritative records for these retrospective audits.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
AUDIT = ROOT / "audit"


def read_json(path: Path):
    return json.loads(path.read_text())


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    DATA.mkdir(exist_ok=True)
    mixed = read_json(AUDIT / "mixed_family" / "records.json")
    pure = read_json(AUDIT / "discovery" / "short_exact_records.json")["records"]

    assert len(mixed) == 459 and all(r.get("master_formula_pass") for r in mixed)
    assert len(pure) == 324 and all(r.get("status") == "success" for r in pure)

    mixed_fields = [
        "word", "split", "word_length", "n_L", "n_D", "n_B",
        "actual_coefficients_json", "predicted_coefficients_json",
        "actual_polynomial_sha256", "predicted_polynomial_sha256",
        "master_formula_pass", "status",
    ]
    write_csv(DATA / "figure4_mixed_family_459.csv", mixed, mixed_fields)

    homogeneous = [next(r for r in mixed if r["word"] == word)
                   for word in ("L", "LL", "D", "DD", "B", "BB")]
    homogeneous_rows = []
    for r in homogeneous:
        row = {k: r.get(k) for k in mixed_fields}
        row["figure4_panels"] = {
            "L": "b", "LL": "c", "D": "e", "DD": "f", "B": "h", "BB": "i"
        }[r["word"]]
        row["family"] = {"L": "cross-linked", "D": "lower-paired", "B": "two-sided"}[r["word"][0]]
        homogeneous_rows.append(row)
    write_csv(DATA / "figure4_homogeneous_m1_m2.csv", homogeneous_rows,
              ["figure4_panels", "family", *mixed_fields])

    pure_fields = [
        "word", "word_length", "n_A", "n_B", "vertices", "edges",
        "selected_crossings", "raw_coefficients_json",
        "normalized_coefficients_json", "raw_polynomial_sha256",
        "normalized_polynomial_sha256", "status",
    ]
    write_csv(DATA / "figure4_pure_braid_324.csv", pure, pure_fields)
    examples = [next(r for r in pure if r["word"] == word)
                for word in ("AAB", "ABA", "BAA")]
    write_csv(DATA / "figure4_order_sensitive_AAB_ABA_BAA.csv", examples, pure_fields)

    aaaba = next(r for r in pure if r["word"] == "AAABA")
    plans = read_json(AUDIT / "discovery" / "word_plans.json")
    assert "AAABA" in plans["required_short_words"]
    assert "AAABA" not in plans["basis"]
    record = {
        "paper_figure": "Supplementary Figure 10",
        "word": "AAABA",
        "ordered_runs": "A^3 B A",
        "hankel_input_entry": "H[AA,ABA]",
        "is_coordinate_basis_word": False,
        "is_independent_holdout": False,
        "source": "audit/discovery/short_exact_records.json",
        "audit_status": "retrospective, not an original 05c ledger",
        "direct_record": aaaba,
    }
    (DATA / "suppfig10_AAABA.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
