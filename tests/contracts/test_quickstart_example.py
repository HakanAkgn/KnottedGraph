"""Verify the copyable 3D example and its invariant across viewing directions."""

import re
import runpy
from pathlib import Path

import numpy as np
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


def test_readme_quickstart_runs_as_copied(capsys):
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    section = readme.split("## Five-minute Quick Start", 1)[1].split("\n## ", 1)[0]
    code = re.search(r"```python\n(.*?)```", section, re.DOTALL).group(1)
    expected_output = re.search(r"```text\n(.*?)```", section, re.DOTALL).group(1)
    exec(compile(code, "README quickstart", "exec"), {})
    assert capsys.readouterr().out == expected_output

    runpy.run_path(str(EXAMPLE_PATH), run_name="__main__")
    assert capsys.readouterr().out == expected_output
