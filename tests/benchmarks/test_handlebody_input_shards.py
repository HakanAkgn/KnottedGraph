import importlib.util
from pathlib import Path

import pytest


def utility():
    root = Path(__file__).resolve().parents[2]
    path = root / "dev/shard_handlebody_inputs.py"
    if not path.exists():
        path = Path(__file__).with_name("shard_handlebody_inputs.py")
    spec = importlib.util.spec_from_file_location("shard_handlebody_inputs", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def inputs(tmp_path):
    module = utility()
    original = {
        "schema_version": 2,
        "config": {"RANDOM_SEED": 20260823, "GRID_SIZE": 200},
        "generation_state": {"seed": 42, "rng_state": [3, [1, 2, 3], None]},
        "cases": [
            {
                "case": f"case_{i}",
                "embedding": {
                    "nodes": [{"id": 0, "pos": [float(i), 0.12345678901234567, -0.25]}],
                    "edges": [],
                },
                "certification": {"radius": 0.001234567890123456},
            }
            for i in range(3)
        ],
    }
    seed = tmp_path / "seed.json.gz"
    recovered = tmp_path / "recovered.json.gz"
    csv = tmp_path / "rows.csv"
    module.write_gzip_json(seed, original)
    module.write_gzip_json(recovered, {"cases": []})
    csv.write_text("index,case\n0,case_0\n")
    return module, original, seed, recovered, csv


def test_sharding_is_deterministic_and_lossless(tmp_path):
    module, original, seed, recovered, csv = inputs(tmp_path)
    one, two = tmp_path / "one", tmp_path / "two"
    first = module.pack(seed, recovered, csv, one, cases_per_shard=2)
    second = module.pack(seed, recovered, csv, two, cases_per_shard=2)
    assert first == second
    assert module.load_manifest(one) == original
    assert module.load_manifest(two / "manifest.json") == original
    for name in [item.name for item in one.iterdir()]:
        assert (one / name).read_bytes() == (two / name).read_bytes()
    assert first["seed_archive"]["metadata"] == {
        key: value for key, value in original.items() if key != "cases"
    }


def test_tampered_shard_is_rejected(tmp_path):
    module, original, seed, recovered, csv = inputs(tmp_path)
    dest = tmp_path / "out"
    manifest = module.pack(seed, recovered, csv, dest, cases_per_shard=2)
    shard = dest / manifest["seed_archive"]["shards"][0]["file"]
    shard.write_bytes(shard.read_bytes() + b"corrupted")
    with pytest.raises(ValueError, match="Shard SHA256 mismatch"):
        module.load_manifest(dest)


def test_companion_hashes_and_existing_outputs_are_protected(tmp_path):
    module, original, seed, recovered, csv = inputs(tmp_path)
    dest = tmp_path / "out"
    module.pack(seed, recovered, csv, dest)
    with pytest.raises(FileExistsError):
        module.pack(seed, recovered, csv, dest)
    (dest / csv.name).write_text("changed\n")
    with pytest.raises(ValueError, match="Companion file SHA256 mismatch"):
        module.load_manifest(dest)
