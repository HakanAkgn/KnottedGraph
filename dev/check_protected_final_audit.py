"""Protect pinned publication implementation and notebook source content.

Git object IDs identify the exact accepted implementation trees and files.
Benchmark notebooks use a canonical source fingerprint that ignores
only explicitly allowlisted bootstrap cells, cell IDs, and transient execution
state.  No baseline commit needs to remain in the current history, so this check
continues to work after a rebase or squash merge.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

# Content-addressed publication implementation objects. Directory
# values are Git tree IDs; file values are Git blob IDs. They do not require the
# commit from which they were recorded to remain in repository history.
EDITORIAL_REFERENCE_REVISION = "9a0f379c38050e80e57967841904917da86c7528"
# The 2026-10-05 prose refresh uses the published main above as its comparison
# snapshot. Computation, notebook code/metadata, and saved data were checked
# against that snapshot before these object IDs and fingerprints were refreshed.
EXPECTED_GIT_OBJECTS = {
    'CMakeLists.txt': 'dc23f7e8ef4456eb92229d9ee93f0c9fb83d546c',
    'dev/final_performance_audit.py': 'd1d6d2568aa77c0b16ad9afc4872efe2daec6ada',
    'dev/final_performance_audit_medium.py': '227e88fbe16c05b7e7c09d1e09cb6fe9323c36c6',
    'src/knotted_graph/extraction': '74661fcc38d0f9bb1e94a16f8543a91697638866',
    'src/knotted_graph/invariants/yamada': '81da4a1d0bdae0388848f8e04a8c6a2687561c90',
    'src/knotted_graph/projection': '3bb3df7404d62b0648ffb4b435a43b408dab6371',
}

# These cells contain environment discovery/reporting rather than the accepted
# mathematical constructions, benchmark cases, or result interpretation.
ALLOWED_BOOTSTRAP_CELLS = {
    "User_guide/benchmarks/01_yamada_sanity_checks.ipynb": {"published"},
    # These notebooks need no source-cell exception. They are listed so that
    # cell IDs and transient execution state can still be normalized without
    # weakening the comparison of their scientific source and metadata.
    "User_guide/benchmarks/02_application_regression_checks.ipynb": set(),
    "User_guide/benchmarks/03_knottedgraph_vs_topoly_scaling.ipynb": set(),
    "User_guide/benchmarks/04_thick_handlebody_validation.ipynb": {"af32d418"},
}

EXPECTED_BENCHMARK_FINGERPRINTS = {
    'User_guide/benchmarks/01_yamada_sanity_checks.ipynb': 'fa593f21b9f3495f43ed49982543573a7f0fe336b9af597c102d5bf4d84db92e',
    'User_guide/benchmarks/02_application_regression_checks.ipynb': 'e3b2c6e00b7bf70dee067c3576556da207d99189c7671df58439c06aa1313aee',
    'User_guide/benchmarks/03_knottedgraph_vs_topoly_scaling.ipynb': '3a02b617e2d52ae77244926091d77693111b31126c7c3d3de56adfa344149354',
    'User_guide/benchmarks/04_thick_handlebody_validation.ipynb': '23740782d8a0d687a110a21b82eafee0a1002f5888cba0e36a42bafa3f60cb98',
}


def _git(*args: str) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode:
        raise RuntimeError(
            f"git {' '.join(args)} failed\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return proc.stdout


def _notebook_at(revision: str, path: str) -> dict:
    return json.loads(_git("show", f"{revision}:{path}"))


def _without_bootstrap_cells(notebook: dict, allowed_ids: set[str]) -> dict:
    normalized = json.loads(json.dumps(notebook))
    normalized["cells"] = [
        cell
        for cell in normalized.get("cells", [])
        if cell.get("id") not in allowed_ids
    ]
    # Cell IDs and transient Jupyter execution state are editor/runtime
    # bookkeeping, not accepted scientific source. Accepted benchmark results
    # live in separately tracked CSV/cache/ground-truth artifacts.
    for cell in normalized["cells"]:
        cell.pop("id", None)
        if cell.get("cell_type") == "code":
            cell["execution_count"] = None
            cell["outputs"] = []
    return normalized


def _benchmark_fingerprint(notebook: dict, allowed_ids: set[str]) -> str:
    canonical = _without_bootstrap_cells(notebook, allowed_ids)
    payload = json.dumps(
        canonical,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def check(revision: str = "HEAD") -> list[str]:
    failures: list[str] = []

    for path, expected in EXPECTED_GIT_OBJECTS.items():
        try:
            actual = _git("rev-parse", f"{revision}:{path}").strip()
        except Exception as exc:
            failures.append(f"could not read protected object {path}: {exc}")
            continue
        if actual != expected:
            failures.append(
                f"protected implementation changed: {path} "
                f"(expected {expected}, found {actual})"
            )

    for path, expected in EXPECTED_BENCHMARK_FINGERPRINTS.items():
        try:
            notebook = _notebook_at(revision, path)
            actual = _benchmark_fingerprint(
                notebook,
                ALLOWED_BOOTSTRAP_CELLS[path],
            )
        except Exception as exc:
            failures.append(f"could not fingerprint protected notebook {path}: {exc}")
            continue
        if actual != expected:
            failures.append(
                f"protected benchmark scientific content changed: {path} "
                f"(expected {expected}, found {actual})"
            )

    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--revision",
        default="HEAD",
        help="Committed tree to audit (default: HEAD).",
    )
    args = parser.parse_args()

    failures = check(args.revision)
    if failures:
        print("PROTECTED FINAL-AUDIT FAILURES:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(
        "PASS: protected implementation objects and benchmark scientific-source "
        "fingerprints match the accepted content-addressed baseline."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
