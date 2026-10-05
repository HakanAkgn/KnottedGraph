"""Verify the copyable 3D example and its invariant across viewing directions."""

import re
import runpy
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from knotted_graph.core.embedding import validate_embedding
from knotted_graph.projection import compute_yamada_polynomial


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_PATH = ROOT / "examples" / "quickstart.py"


def test_quickstart_uses_a_3d_closed_curve_and_preserves_its_invariant():
    namespace = runpy.run_path(str(EXAMPLE_PATH))
    graph, result = namespace["compute_quickstart"]()
    points = graph.edges["loop_anchor", "loop_anchor", "curve"]["pts"]

    assert validate_embedding(graph) == []
    assert np.linalg.matrix_rank(points - points.mean(axis=0)) == 3
    np.testing.assert_array_equal(points[0], points[-1])
    assert result.projection.num_crossings == 3

    Y = sp.Symbol("Y")
    expected = -Y**11 + Y**9 + Y**8 + Y**7 - Y**4 - Y**3 - Y**2 - Y - 1
    assert sp.expand(result.polynomial - expected) == 0

    other_view = compute_yamada_polynomial(
        graph, Y, rotation_angles=(25, 10, 15), n_jobs=1, return_result=True,
    )
    assert sp.expand(other_view.polynomial - result.polynomial) == 0


def test_volume_quickstart_extracts_a_branched_graph_and_preserves_its_invariant():
    pytest.importorskip("skimage")
    from skimage.measure import euler_number

    namespace = runpy.run_path(str(ROOT / "examples" / "volume_to_yamada.py"))
    volume, graph, result = namespace["compute_volume_quickstart"]()

    assert volume.dtype == np.bool_
    assert euler_number(volume, connectivity=3) == -1
    assert graph.number_of_nodes() == 2
    assert graph.number_of_edges() == 3
    assert sorted(degree for _, degree in graph.degree()) == [3, 3]
    assert validate_embedding(graph) == []
    assert result.projection.num_crossings == 4
    assert result.projection.pd_code.count("V[") == 2
    assert result.projection.pd_code.count("X[") == 4

    Y = sp.Symbol("Y")
    expected = Y**12 - Y**8 - Y**6 - Y**4 - Y**3 - Y**2 - Y - 1
    assert sp.expand(result.polynomial - expected) == 0
    other_view = compute_yamada_polynomial(
        graph, Y, rotation_angles=(25, 10, 15), n_jobs=1, return_result=True,
    )
    assert sp.expand(other_view.polynomial - result.polynomial) == 0


def test_readme_volume_quickstart_runs_as_copied(capsys):
    pytest.importorskip("skimage")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    section = readme.split("## Five-minute Quick Start", 1)[1].split("\n## ", 1)[0]
    code = re.search(r"```python\n(.*?)```", section, re.DOTALL).group(1)
    expected_output = re.search(r"```text\n(.*?)```", section, re.DOTALL).group(1)
    exec(compile(code, "README quickstart", "exec"), {})
    assert capsys.readouterr().out == expected_output

    runpy.run_path(str(ROOT / "examples" / "volume_to_yamada.py"), run_name="__main__")
    assert capsys.readouterr().out == expected_output


def test_website_volume_quickstart_runs_as_copied(capsys):
    pytest.importorskip("skimage")
    page = (ROOT / "doc" / "quickstart.md").read_text(encoding="utf-8")
    section = page.split("## Explore the graph in 3D", 1)[0]
    code = "\n".join(re.findall(r"```python\n(.*?)```", section, re.DOTALL))
    expected_output = re.search(r"```text\n(.*?)```", section, re.DOTALL).group(1)
    exec(compile(code, "website quickstart", "exec"), {})
    assert capsys.readouterr().out == expected_output

    curve_section = page.split("## Coordinate curves: a trefoil", 1)[1]
    curve_code = "\n".join(re.findall(r"```python\n(.*?)```", curve_section, re.DOTALL))
    curve_output = re.search(r"```text\n(.*?)```", curve_section, re.DOTALL).group(1)
    exec(compile(curve_code, "website curve example", "exec"), {})
    assert capsys.readouterr().out == curve_output
