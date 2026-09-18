#!/usr/bin/env python3
"""Losslessly package the immutable historical handlebody inputs below GitHub's file limit.

No embeddings are generated or reassociated. The archive content digest hashes
canonical metadata JSON followed by a newline, then each canonical case JSON
and a newline in original order. Reassembly preserves the JSON structure, not
the original gzip header or whitespace. Original compressed-file SHA256 is
retained separately in the manifest.
"""

from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil


def canonical_json(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def content_sha256(archive):
    metadata = {key: value for key, value in archive.items() if key != "cases"}
    digest = hashlib.sha256(canonical_json(metadata) + b"\n")
    for case in archive["cases"]:
        digest.update(canonical_json(case))
        digest.update(b"\n")
    return digest.hexdigest()


def write_gzip_json(path, value):
    # Empty filename and zero modification time eliminate variable gzip headers.
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(
            fileobj=raw, filename="", mode="wb", mtime=0, compresslevel=9
        ) as compressed:
            compressed.write(canonical_json(value))


def read_gzip_json(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def checked_child(root, name):
    candidate = (root / name).resolve()
    if not candidate.is_relative_to(root.resolve()):
        raise ValueError("Manifest path escapes its input directory")
    return candidate


def load_manifest(path):
    path = Path(path)
    if path.is_dir():
        path = path / "manifest.json"
    manifest = json.loads(path.read_text())
    if manifest.get("format") != "knotted-graph-handlebody-input-shards-v1":
        raise ValueError("Unsupported shard manifest")
    cases = []
    for item in manifest["seed_archive"]["shards"]:
        source = checked_child(path.parent, item["file"])
        if sha256_file(source) != item["sha256"]:
            raise ValueError(f"Shard SHA256 mismatch: {item['file']}")
        chunk = read_gzip_json(source)["cases"]
        if len(chunk) != item["count"] or item["start_index"] != len(cases):
            raise ValueError("Shard count/order mismatch")
        cases.extend(chunk)
    archive = dict(manifest["seed_archive"]["metadata"], cases=cases)
    if len(cases) != manifest["seed_archive"]["count"]:
        raise ValueError("Seed archive count mismatch")
    if content_sha256(archive) != manifest["seed_archive"]["content_sha256"]:
        raise ValueError("Seed archive content digest mismatch")
    for item in manifest["files"].values():
        if sha256_file(checked_child(path.parent, item["file"])) != item["sha256"]:
            raise ValueError(f"Companion file SHA256 mismatch: {item['file']}")
    return archive


def pack(
    seed_path,
    recovered_path,
    csv_path,
    output_dir,
    *,
    cases_per_shard=500,
    association_report=None,
):
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError("Refusing a nonempty output directory")
    if cases_per_shard < 1:
        raise ValueError("cases_per_shard must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)
    source = read_gzip_json(seed_path)
    metadata = {key: value for key, value in source.items() if key != "cases"}
    shards = []
    for start in range(0, len(source["cases"]), cases_per_shard):
        chunk = source["cases"][start : start + cases_per_shard]
        name = f"seed_cases_{start:05d}_{start + len(chunk) - 1:05d}.json.gz"
        target = output_dir / name
        write_gzip_json(target, {"cases": chunk})
        if target.stat().st_size >= 50_000_000:
            raise ValueError(
                "Shard exceeds 50 MB; rerun in a new directory with fewer cases per shard"
            )
        shards.append(
            {
                "file": name,
                "start_index": start,
                "count": len(chunk),
                "sha256": sha256_file(target),
                "bytes": target.stat().st_size,
            }
        )
    files = {}
    for kind, path in [("recovered_archive", recovered_path), ("csv", csv_path)]:
        path = Path(path)
        target = output_dir / path.name
        shutil.copyfile(path, target)
        files[kind] = {
            "file": target.name,
            "sha256": sha256_file(target),
            "bytes": target.stat().st_size,
        }
    conflicts = []
    if association_report:
        report = json.loads(Path(association_report).read_text())
        conflicts = [row for row in report if not row["consistent"]]
        target = output_dir / "seed_association_conflicts.json"
        target.write_bytes(canonical_json(conflicts) + b"\n")
        files["association_conflicts"] = {
            "file": target.name,
            "sha256": sha256_file(target),
            "bytes": target.stat().st_size,
        }
    manifest = {
        "format": "knotted-graph-handlebody-input-shards-v1",
        "seed_archive": {
            "original_filename": Path(seed_path).name,
            "original_file_sha256": sha256_file(seed_path),
            "metadata": metadata,
            "count": len(source["cases"]),
            "content_sha256": content_sha256(source),
            "content_digest_scheme": "sha256(canonical_metadata_json + LF + ordered_canonical_case_json_and_LF)",
            "shards": shards,
        },
        "files": files,
        "known_seed_association_conflicts": len(conflicts),
        "association_note": "Historical case names are not unique input identities. Conflicting seed/CSV associations are preserved and must be withheld, not guessed or repaired. Historical rows lack seed payload hashes; metadata compatibility does not independently certify original seed coordinates.",
    }
    (output_dir / "manifest.json").write_bytes(canonical_json(manifest) + b"\n")
    # Read every generated shard and companion file back, checking semantic identity.
    rebuilt = load_manifest(output_dir)
    if rebuilt != source:
        raise AssertionError(
            "Reassembled JSON structure differs from the original archive"
        )
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    packing = sub.add_parser("pack")
    for name in ("seed-archive", "recovered-archive", "csv", "output-dir"):
        packing.add_argument("--" + name, type=Path, required=True)
    packing.add_argument("--cases-per-shard", type=int, default=500)
    packing.add_argument("--association-report", type=Path)
    verify = sub.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    rebuild = sub.add_parser("reassemble")
    rebuild.add_argument("--manifest", type=Path, required=True)
    rebuild.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.operation == "pack":
        manifest = pack(
            args.seed_archive,
            args.recovered_archive,
            args.csv,
            args.output_dir,
            cases_per_shard=args.cases_per_shard,
            association_report=args.association_report,
        )
        print(
            json.dumps(
                {
                    "cases": manifest["seed_archive"]["count"],
                    "shards": len(manifest["seed_archive"]["shards"]),
                    "maximum_shard_bytes": max(
                        item["bytes"] for item in manifest["seed_archive"]["shards"]
                    ),
                    "conflicts": manifest["known_seed_association_conflicts"],
                }
            )
        )
    else:
        archive = load_manifest(args.manifest)
        if args.operation == "reassemble":
            if args.output.exists():
                parser.error("Refusing to overwrite an existing reassembled archive")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            write_gzip_json(args.output, archive)
        print(
            json.dumps(
                {
                    "verified_cases": len(archive["cases"]),
                    "content_sha256": content_sha256(archive),
                }
            )
        )


if __name__ == "__main__":
    main()
