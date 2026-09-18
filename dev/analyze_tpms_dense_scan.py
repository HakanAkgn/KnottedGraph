"""Read-only audit of corrected TPMS records and the existing display filter.

No figures or manuscript files are written. The output is operational evidence,
not a topology classification or an embedding-equivalence certificate.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict, deque
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_actual_filter(source):
    tree = ast.parse(source.read_text())
    wanted = {"stable_labels", "connected_components_for_label"}
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
    if {n.name for n in nodes} != wanted:
        raise ValueError("Cannot locate the two original filter functions")
    namespace = {"np": np, "Counter": Counter, "deque": deque}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"), namespace)
    return namespace["stable_labels"]


def cell_key(r):
    return {"family": r["family"], "lambda": r["lam"], "c": r["threshold_c"]}


def scope(r):
    if r["source"] == "error":
        return "unresolved-error"
    if r["source"] == "large-core":
        return "abstract-graph-summary"
    if r["source"] == "diagram-yamada":
        return "fixed-diagram-polynomial"
    if r["source"] in ("yamada", "vertex"):
        return "subcubic-normalized-yamada"
    return "unrecognized-" + r["source"]


def summarize(records):
    groups = defaultdict(list)
    for r in records:
        groups[r["source"]].append(r)
    classification_ok = [r for r in records if r["source"] != "error"]
    fully_ok = [r for r in classification_ok if not r.get("error")]
    exact = [r for r in fully_ok if r["source"] in ("yamada", "vertex")]
    diagram = [r for r in fully_ok if r["source"] == "diagram-yamada"]
    voids = [r for r in records if r.get("void_components", 0)]
    failed = [r for r in records if r["source"] == "error"]
    return {
        "cells": len(records),
        "source_counts": dict(sorted(Counter(r["source"] for r in records).items())),
        "scope_counts": dict(sorted(Counter(scope(r) for r in records).items())),
        "source_unique_signature_counts": {k: len({r["phase_signature"] for r in v}) for k, v in sorted(groups.items())},
        "raw_signature_count_including_errors": len({r["phase_signature"] for r in records}),
        "raw_signature_count_excluding_source_errors": len({r["phase_signature"] for r in classification_ok}),
        "raw_signature_count_excluding_any_error": len({r["phase_signature"] for r in fully_ok}),
        "subcubic_normalized_polynomial_cells": len(exact),
        "subcubic_normalized_distinct_polynomials": len({r["polynomial"] for r in exact}),
        "fixed_diagram_polynomial_cells": len(diagram),
        "fixed_diagram_distinct_polynomials": len({r["polynomial"] for r in diagram}),
        "classification_error_cells": len(failed),
        "any_error_cells": sum(bool(r.get("error")) for r in records),
        "error_messages": dict(sorted(Counter(r["error"] for r in records if r.get("error")).items())),
        "cavity_positive_cells": len(voids),
        "cavity_positive_total_voids": sum(r["void_components"] for r in voids),
        "cavity_cells": [dict(cell_key(r), void_components=r["void_components"], error=r.get("error")) for r in voids],
        "remaining_nonvoid_classification_errors": [dict(cell_key(r), error=r.get("error")) for r in failed if not r.get("void_components", 0)],
        "b0_mismatches_in_returned_graphs": [dict(cell_key(r), source=r["source"], mask=r["interior_components"], graph=r["components"]) for r in records if r["interior_components"] != r["components"]],
        "b1_mismatches_in_returned_graphs": [dict(cell_key(r), source=r["source"], mask=r["handle_rank"], graph=r["cycle_rank"]) for r in records if r["handle_rank"] != r["cycle_rank"]],
        "successful_cells_b0_mismatch_count": sum(r["interior_components"] != r["components"] for r in classification_ok),
        "successful_cells_b1_mismatch_count": sum(r["handle_rank"] != r["cycle_rank"] for r in classification_ok),
        "component_validation_unreached_cells": sum(r.get("component_count_matches") is None for r in records),
        "classification_computed_flag_false_cells": sum(r.get("classification_computed") is False for r in records),
        "closed_surface_cells": sum(bool(r["surface_is_closed"]) for r in records),
        "open_surface_cells": sum(not r["surface_is_closed"] for r in records),
        "box_boundary_touch_cells": sum(bool(r["touches_boundary"]) for r in records),
        "evaluation_kind_counts": dict(sorted(Counter(r.get("evaluation_kind", "missing") for r in records).items())),
        "exact_yamada_attempted_cells": sum(r["exact_yamada_attempted"] for r in records),
        "zero_polynomial_cells": [dict(cell_key(r), source=r["source"]) for r in records if r.get("polynomial") == "0"],
    }


def filter_sensitivity(records, stable_labels):
    signatures = list(dict.fromkeys(r["phase_signature"] for r in records))
    ids = {s: i + 1 for i, s in enumerate(signatures)}
    by_id = {ids[r["phase_signature"]]: r for r in records}
    lambdas = sorted({r["lam"] for r in records})
    thresholds = sorted({r["threshold_c"] for r in records})
    lookup = {(round(r["threshold_c"], 12), round(r["lam"], 12)): r for r in records}
    grid_records = [[lookup[round(c, 12), round(lam, 12)] for lam in lambdas] for c in thresholds]
    raw = np.array([[ids[r["phase_signature"]] for r in row] for row in grid_records], dtype=int)
    output = {
        "axis_order": "rows=threshold_c, columns=lambda",
        "lambdas": lambdas,
        "thresholds": thresholds,
        "raw_label_grid": raw.tolist(),
        "raw_source_grid": [[r["source"] for r in row] for row in grid_records],
        "raw_void_count_grid": [[r["void_components"] for r in row] for row in grid_records],
        "label_lookup": {str(i): {"signature": r["phase_signature"], "source": r["source"], "scope": scope(r)} for i, r in by_id.items()},
        "thresholds_min_cells": {},
    }
    for threshold in (1, 2, 4, 8):
        stable, edit_count = stable_labels(raw, min_cells=threshold)
        changed = stable != raw
        error_mask = np.array([[r["source"] == "error" for r in row] for row in grid_records])
        new_error_mask = np.vectorize(lambda i: by_id[int(i)]["source"] == "error")(stable)
        reversed_stable, _ = stable_labels(len(ids) + 1 - raw, min_cells=threshold)
        reversed_restored = len(ids) + 1 - reversed_stable
        changed_cells = []
        for i, j in zip(*np.nonzero(changed)):
            r = grid_records[i][j]
            replacement = by_id[int(stable[i, j])]
            changed_cells.append(dict(cell_key(r), row=int(i), col=int(j), original_id=int(raw[i,j]), final_id=int(stable[i,j]), original_source=r["source"], final_source=replacement["source"], original_void_components=r["void_components"]))
        final_ids = {int(v) for v in stable.ravel()}
        output["thresholds_min_cells"][str(threshold)] = {
            "final_label_grid": stable.tolist(),
            "final_relabel_mask": changed.tolist(),
            "changed_cells": changed_cells,
            "unique_final_changed_cells": int(changed.sum()),
            "algorithm_edit_count_including_repeated_edits": edit_count,
            "final_label_count_including_errors": len(final_ids),
            "final_label_count_excluding_source_errors": sum(by_id[i]["source"] != "error" for i in final_ids),
            "final_label_counts_by_source": dict(sorted(Counter(by_id[i]["source"] for i in final_ids).items())),
            "final_cell_counts_by_source_of_assigned_label": dict(sorted(Counter(by_id[int(i)]["source"] for i in stable.ravel()).items())),
            "changed_cells_crossing_signature_scope": sum(scope(grid_records[i][j]) != scope(by_id[int(stable[i,j])]) for i,j in zip(*np.nonzero(changed))),
            "original_error_cells_given_nonerror_label": int(np.logical_and(error_mask, ~new_error_mask).sum()),
            "original_nonerror_cells_given_error_label": int(np.logical_and(~error_mask, new_error_mask).sum()),
            "same_filter_after_reversing_label_ids_different_final_cells": int(np.sum(reversed_restored != stable)),
        }
    return output


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--records", type=Path, required=True)
    p.add_argument("--filter-source", type=Path, required=True)
    p.add_argument("--output-prefix", type=Path, required=True)
    args = p.parse_args()
    records = json.loads(args.records.read_text())
    stable = load_actual_filter(args.filter_source)
    family_records = defaultdict(list)
    for r in records:
        family_records[r["family"]].append(r)
    provenance = {
        "records_path": str(args.records.resolve()),
        "records_sha256": sha256(args.records),
        "filter_source_path": str(args.filter_source.resolve()),
        "filter_source_sha256": sha256(args.filter_source),
        "scope": "Operational raw signatures; equality is not embedding equivalence. Beta matching is necessary, not a spine proof. Filter is a display operation and cannot resolve failures.",
    }
    summary = dict(provenance, aggregate=summarize(records), families={k: summarize(v) for k,v in family_records.items()})
    sensitivity = dict(provenance, filter_rule="The source's original in-place 4-neighbor filter, up to 3 passes, with sequential sorted label processing. Threshold 1 is unfiltered. Label IDs follow first occurrence in the records.", families={k: filter_sensitivity(v, stable) for k,v in family_records.items()})
    summary_path = Path(str(args.output_prefix) + "_summary.json")
    sensitivity_path = Path(str(args.output_prefix) + "_filter_sensitivity.json")
    summary_path.write_text(json.dumps(summary, indent=2))
    sensitivity_path.write_text(json.dumps(sensitivity, indent=2))
    compact = {
        "summary": str(summary_path), "sensitivity": str(sensitivity_path),
        "families": {k: {"sources": v["source_counts"], "valid_raw_signatures": v["raw_signature_count_excluding_source_errors"], "subcubic_polynomials": v["subcubic_normalized_distinct_polynomials"], "diagram_polynomials": v["fixed_diagram_distinct_polynomials"], "void_cells": v["cavity_positive_cells"], "other_errors": len(v["remaining_nonvoid_classification_errors"])} for k,v in summary["families"].items()},
        "filter": {k: {t: {key: row[key] for key in ("unique_final_changed_cells", "algorithm_edit_count_including_repeated_edits", "final_label_count_excluding_source_errors", "original_error_cells_given_nonerror_label", "same_filter_after_reversing_label_ids_different_final_cells")} for t,row in v["thresholds_min_cells"].items()} for k,v in sensitivity["families"].items()},
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
