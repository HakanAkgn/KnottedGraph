"""Read saved phase records without re-running extraction or evaluating formulas."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any


def _boolean(value: Any, *, field: str, row: int) -> bool:
    if value is True or value == "True" or value == "true" or value == "1":
        return True
    if value is False or value == "False" or value == "false" or value == "0":
        return False
    raise ValueError(f"Row {row}: {field} must be a boolean, got {value!r}.")


def read_phase_map_records(path: str | Path) -> tuple[dict[str, Any], ...]:
    """Load a material/TPMS CSV or JSON record list and validate coordinates.

    Polynomial and signature fields remain text: no ``eval`` or ``sympify`` is
    used on files. CSV is preferable for browsing the large reference scans.
    Duplicate coordinates, non-finite coordinates and incomplete required
    fields are errors. Missing grid cells are allowed and plotted as missing.
    """
    path = Path(path)
    if path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
    elif path.suffix.lower() == ".json":
        rows = json.loads(path.read_text(encoding="utf-8"))
    else:
        raise ValueError(
            "Use a *_records.csv or *_records.json file, not an HTML or summary file."
        )
    if not isinstance(rows, list) or not rows:
        raise ValueError("Expected a non-empty list of phase-map records.")
    records: list[dict[str, Any]] = []
    seen: set[tuple[str, str, float, float]] = set()
    for index, source in enumerate(rows, 1):
        if not isinstance(source, dict):
            raise ValueError(f"Row {index}: expected a record object.")
        record = dict(source)
        family_field, level_field = (
            ("material", "energy")
            if "material" in record
            else ("family", "threshold_c")
        )
        for name in (family_field, "lam", level_field, "phase_signature", "source"):
            if name not in record or record[name] in (None, ""):
                raise ValueError(f"Row {index}: missing required field {name!r}.")
        for name in ("lam", level_field):
            try:
                record[name] = float(record[name])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Row {index}: {name} must be a finite number."
                ) from exc
            if not math.isfinite(record[name]):
                raise ValueError(f"Row {index}: {name} must be a finite number.")
        for name in (family_field, "phase_signature", "source"):
            if not isinstance(record[name], str):
                raise ValueError(f"Row {index}: {name} must be text.")
        if "classification_computed" in record:
            record["classification_computed"] = _boolean(
                record["classification_computed"],
                field="classification_computed",
                row=index,
            )
        key = (family_field, record[family_field], record["lam"], record[level_field])
        if key in seen:
            raise ValueError(f"Row {index}: duplicate phase-map coordinate {key!r}.")
        seen.add(key)
        records.append(record)
    return tuple(records)


@dataclass(frozen=True)
class PhaseMapData:
    """One family of sampled records; labels are categories, not invariant values."""

    family: str
    level_field: str
    records: tuple[dict[str, Any], ...]

    @classmethod
    def from_records(
        cls, records: tuple[dict[str, Any], ...], family: str | None = None
    ) -> PhaseMapData:
        """Select a family from records validated by :func:`read_phase_map_records`."""
        keys = sorted({str(r.get("material", r.get("family"))) for r in records})
        if family is None:
            if len(keys) != 1:
                raise ValueError(
                    "Select one family with --family (or family=): " + ", ".join(keys)
                )
            family = keys[0]
        if family not in keys:
            raise ValueError(f"Unknown family {family!r}. Available: {', '.join(keys)}")
        subset = tuple(
            r for r in records if r.get("material", r.get("family")) == family
        )
        fields = {"energy" if "material" in r else "threshold_c" for r in subset}
        if len(fields) != 1:
            raise ValueError(f"Family {family!r} mixes material and TPMS schemas.")
        return cls(family, fields.pop(), subset)

    def summary(self) -> dict[str, Any]:
        """Return JSON-serializable counts, preserving error and adaptive-fill status."""
        nx = len({r["lam"] for r in self.records})
        ny = len({r[self.level_field] for r in self.records})
        return {
            "family": self.family,
            "level_field": self.level_field,
            "records": len(self.records),
            "lambda_samples": nx,
            "level_samples": ny,
            "missing_grid_cells": nx * ny - len(self.records),
            "distinct_signatures": len({r["phase_signature"] for r in self.records}),
            "source_counts": dict(
                sorted(Counter(r["source"] for r in self.records).items())
            ),
            "error_records": sum(
                bool(r.get("error")) or r["source"] == "error" for r in self.records
            ),
            "adaptive_fill_records": sum(
                r.get("classification_computed") is False for r in self.records
            ),
            "classification_status_unrecorded": sum(
                "classification_computed" not in r for r in self.records
            ),
            "display_processing": "raw saved signatures; no smoothing or signature merging",
        }


def load_phase_map(path: str | Path, *, family: str | None = None) -> PhaseMapData:
    """Read and select a family; use ``family=`` for files containing several scans."""
    return PhaseMapData.from_records(read_phase_map_records(path), family)
