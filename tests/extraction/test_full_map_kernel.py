"""Differential checks before allowing a fast kernel to produce map data."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from knotted_graph.extraction.cubical_spine import cubical_retract, verify_cubical_retract

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("full_map_kernel", ROOT / "dev/full_map_kernel.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.fixture(scope="module")
def kernel(tmp_path_factory):
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("optional full-map runner requires a C++17 compiler")
    binary = tmp_path_factory.mktemp("full-map-build") / "collapse"
    subprocess.run([compiler, "-O3", "-std=c++17", str(ROOT / "dev/full_map_collapse.cpp"),
                    "-o", str(binary)], check=True)
    return binary


@pytest.mark.parametrize("seed", range(30))
def test_identical_reference_collapse_sequence(kernel, seed):
    rng = np.random.default_rng(seed)
    mask = rng.random((3 + seed % 3, 4, 3)) < .2 + .02 * seed
    expected = cubical_retract(mask)
    graph, metadata, moves, cells = MODULE.reconstruct(mask, kernel)
    certificate = MODULE.expand_reference_certificate(mask, None, metadata, moves, cells)
    assert certificate["collapses"] == expected.certificate["collapses"]
    assert certificate["terminal_cells"] == expected.certificate["terminal_cells"]
    assert certificate["kind"] == expected.certificate["kind"]
    assert verify_cubical_retract(mask, certificate)["valid"]
    assert (graph is None) == (expected.graph is None)


@pytest.mark.parametrize("kind", ["block", "loop", "shell", "isolates", "boundary", "empty"])
def test_structured_inputs_and_archive(kernel, tmp_path, kind):
    mask = np.zeros((7, 7, 7), dtype=bool)
    if kind in ("block", "loop", "shell"):
        mask[1:6, 1:6, 1:6] = True
    if kind == "loop":
        mask[2:5, 2:5, :] = False
    if kind == "shell":
        mask[2:5, 2:5, 2:5] = False
    if kind == "isolates":
        mask[1, 1, 1] = mask[5, 5, 5] = True
    if kind == "boundary":
        mask[:3, :3, :3] = True
    archive = tmp_path / "witness.npz"
    graph, metadata, moves, cells = MODULE.reconstruct(mask, kernel, archive=archive)
    certificate = MODULE.expand_reference_certificate(mask, None, metadata, moves, cells)
    assert verify_cubical_retract(mask, certificate)["valid"]
    with np.load(archive, allow_pickle=False) as saved:
        recovered = np.cumsum(saved["move_delta"].astype(np.int64), axis=0).astype(np.uint32)
        assert np.array_equal(recovered, moves)
        assert np.array_equal(saved["mask"], mask)
    if kind == "shell":
        assert graph is None
    assert metadata["analytic_source_correspondence_certified"] is False
