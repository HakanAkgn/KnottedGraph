"""User-facing contracts for reading and running the new phase-map examples."""

import csv
import importlib.util
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest

from knotted_graph.applications.phase_map_examples import (
    load_phase_map,
    plot_phase_map,
    read_phase_map_records,
)
from knotted_graph.applications.phase_map_examples.cli import main, parser, scan_plan


def row(**updates):
    result = {
        "family": "gyroid_to_diamond",
        "lam": 0,
        "threshold_c": 0,
        "source": "yamada",
        "phase_signature": "yamada:-Y-1",
        "polynomial": "-Y-1",
        "error": "",
    }
    result.update(updates)
    return result


def write_records(tmp_path, records, suffix=".json"):
    path = tmp_path / ("records" + suffix)
    if suffix == ".csv":
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    else:
        path.write_text(json.dumps(records))
    return path


@pytest.mark.parametrize("suffix", [".json", ".csv"])
def test_read_and_summarize_preserves_classification(tmp_path, suffix):
    rows = [
        row(classification_computed=True),
        row(
            lam=1,
            classification_computed=False,
            source="large-core",
            phase_signature="core:n=30",
        ),
        row(
            threshold_c=1,
            classification_computed=True,
            source="error",
            error="extraction failed",
            phase_signature="error:failed",
            polynomial="",
        ),
    ]
    data = load_phase_map(write_records(tmp_path, rows, suffix))
    summary = data.summary()
    assert summary["records"] == 3
    assert summary["missing_grid_cells"] == 1
    assert summary["adaptive_fill_records"] == 1
    assert summary["error_records"] == 1
    assert summary["distinct_signatures"] == 3
    assert summary["source_counts"] == {"error": 1, "large-core": 1, "yamada": 1}


@pytest.mark.parametrize(
    "updates, message",
    [
        ({"lam": "nan"}, "finite"),
        ({"threshold_c": "inf"}, "finite"),
        ({"lam": "not a number"}, "finite"),
        ({"phase_signature": ""}, "required"),
        ({"source": 2}, "text"),
        ({"classification_computed": "perhaps"}, "boolean"),
        ({"resolution_calibration_energy": "nan"}, "finite"),
        ({"resolution_calibration_energy": "anchor"}, "finite"),
    ],
)
def test_invalid_records_have_actionable_errors(tmp_path, updates, message):
    with pytest.raises(ValueError, match=message):
        read_phase_map_records(write_records(tmp_path, [row(**updates)]))


def test_duplicate_empty_wrong_schema_and_extension(tmp_path):
    for records, message in [
        ([row(), row()], "duplicate"),
        ([], "non-empty"),
        ({"families": []}, "non-empty"),
        ([1], "record object"),
    ]:
        with pytest.raises(ValueError, match=message):
            read_phase_map_records(write_records(tmp_path, records))
    with pytest.raises(ValueError, match="records.csv"):
        read_phase_map_records(tmp_path / "result.html")


def test_family_selection_and_material_axis(tmp_path):
    material = row(material="tib2_d6_F", energy=0.35)
    del material["family"]
    del material["threshold_c"]
    path = write_records(tmp_path, [row(), material])
    with pytest.raises(ValueError, match="Select one family"):
        load_phase_map(path)
    with pytest.raises(ValueError, match="Unknown family"):
        load_phase_map(path, family="missing")
    assert load_phase_map(path, family="tib2_d6_F").level_field == "energy"


@pytest.mark.parametrize("suffix", [".json", ".csv"])
def test_calibration_is_distinct_from_adaptive_fill_and_zero_anchor_is_valid(
    tmp_path, suffix
):
    rows = [
        row(classification_computed=True, resolution_calibration_energy=None),
        row(lam=1, classification_computed=False, resolution_calibration_energy=None),
        row(
            threshold_c=1,
            classification_computed=False,
            resolution_calibration_energy=0.0,
            source="yamada-set",
            phase_signature="yamada-set:outer[Y]inner[-1]",
            polynomial="{Y,-1}",
        ),
    ]
    data = load_phase_map(write_records(tmp_path, rows, suffix))
    summary = data.summary()
    assert summary["directly_classified_records"] == 1
    assert summary["adaptive_fill_records"] == 1
    assert summary["resolution_calibration_records"] == 1
    assert summary["classification_status_unrecorded"] == 0
    assert data.records[2]["resolution_calibration_energy"] == 0.0
    fig = plot_phase_map(data, output=tmp_path / "provenance")
    plt.close(fig)
    exported = json.loads((tmp_path / "provenance.json").read_text())
    assert exported["classification_status"] == [
        ["computed", "adaptive_fill"],
        ["resolution_calibration", "missing"],
    ]
    assert exported["resolution_calibration_energy"] == [[None, None], [0.0, None]]
    assert data.records[2]["polynomial"] == "{Y,-1}"


def test_plot_is_raw_and_math_style_scoped(tmp_path):
    expression = "__import__('os').system('this must never execute')"
    path = write_records(
        tmp_path, [row(polynomial=expression), row(lam=1, threshold_c=1)]
    )
    before = mpl.rcParams["mathtext.fontset"]
    fig = plot_phase_map(load_phase_map(path), output=tmp_path / "figure.png")
    assert fig.axes[0].xaxis.label.get_math_fontfamily() == "cm"
    assert mpl.rcParams["mathtext.fontset"] == before
    plt.close(fig)
    assert (tmp_path / "figure.png").stat().st_size > 1000
    assert (tmp_path / "figure.pdf").read_bytes().startswith(b"%PDF")
    metadata = json.loads((tmp_path / "figure.json").read_text())
    assert metadata["phase_key"]["1"]["polynomials"] == ["-Y-1", expression]
    assert metadata["summary"]["missing_grid_cells"] == 2
    assert metadata["signature_ids"] == [[1, None], [None, 1]]
    assert metadata["phase_key"]["1"]["color"].startswith("#")


def test_inspect_and_plot_cli(tmp_path, capsys):
    path = write_records(tmp_path, [row(), row(lam=1)])
    main(["inspect", str(path)])
    assert json.loads(capsys.readouterr().out)[0]["records"] == 2
    main(["plot", str(path), "--output", str(tmp_path / "cli")])
    assert (tmp_path / "cli.pdf").exists()


@pytest.mark.parametrize("kind", ["materials", "tpms"])
def test_dry_run_is_bounded_and_has_no_side_effects(tmp_path, capsys, kind):
    output = tmp_path / "unused"
    main(["scan", kind, "--output-dir", str(output), "--dry-run"])
    plan = json.loads(capsys.readouterr().out)
    assert plan["dimension"] == 24
    assert plan["lambda_count"] == plan["level_count"] == 3
    assert plan["workers"] == 1
    assert not output.exists()


@pytest.mark.parametrize(
    "arguments, message",
    [
        (["--family", "wrong"], "Unknown"),
        (["--dimension", "7"], "dimension"),
        (["--lambda-count", "1"], "samples"),
        (["--level-count", "1"], "samples"),
        (["--workers", "2"], "serial"),
    ],
)
def test_scan_plan_rejects_invalid_configuration(arguments, message):
    args = parser().parse_args(["scan", "tpms", "--output-dir", "unused", *arguments])
    with pytest.raises(ValueError, match=message):
        scan_plan(args)


def test_scan_refuses_existing_results_before_optional_imports(tmp_path):
    path = write_records(tmp_path, [row()])
    original = path.read_bytes()
    with pytest.raises(SystemExit) as exc:
        main(["scan", "tpms", "--output-dir", str(tmp_path)])
    assert exc.value.code == 2
    assert path.read_bytes() == original


def test_missing_scan_dependencies_explain_base_route(monkeypatch):
    from knotted_graph.applications.phase_map_examples._runtime import (
        require_scan_dependencies,
    )

    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
    with pytest.raises(ImportError, match="inspect and plot"):
        require_scan_dependencies()


@pytest.mark.parametrize("archived", [True, False])
def test_geometry_paths_follow_selected_dataset(tmp_path, archived):
    from knotted_graph.applications.phase_map_examples._runtime import (
        resolve_geometry_path,
    )

    scan = tmp_path / "data" if archived else tmp_path
    scan.mkdir(exist_ok=True)
    geometry = tmp_path / "geometry"
    geometry.mkdir()
    expected = geometry / "cell.json.gz"
    expected.write_bytes(b"example")
    assert (
        resolve_geometry_path("/old/machine/tmp/geometry/cell.json.gz", scan)
        == expected
    )
    assert resolve_geometry_path(r"C:\old\geometry\cell.json.gz", scan) == expected
    with pytest.raises(FileNotFoundError, match="Missing geometry"):
        resolve_geometry_path("/old/machine/tmp/geometry/missing.json.gz", scan)


def test_missing_new_scan_geometry_does_not_borrow_another_scan(tmp_path):
    from knotted_graph.applications.phase_map_examples._runtime import (
        resolve_geometry_path,
    )

    scan = tmp_path / "new_scan"
    scan.mkdir()
    (tmp_path / "geometry").mkdir()
    (tmp_path / "geometry/cell.json.gz").write_bytes(b"a different dataset")
    with pytest.raises(FileNotFoundError, match="Missing geometry"):
        resolve_geometry_path("geometry/cell.json.gz", scan)


def test_cli_families_match_engines():
    pytest.importorskip("pyvista")
    pytest.importorskip("poly2graph")
    pytest.importorskip("plotly")
    from knotted_graph.applications.phase_map_examples import _materials, _tpms
    from knotted_graph.applications.phase_map_examples.cli import FAMILIES

    assert set(FAMILIES["materials"]) == {
        f.key for f in _materials.material_families(8)
    }
    assert set(FAMILIES["tpms"]) == {
        f.key for f in _tpms.tpms_families(dimension=8, thresholds=(0.0, 0.3))
    }
    defaults = _materials.parse_args([])
    assert defaults.adaptive_energy_step == 0
    assert defaults.min_stable_cells == 1
    assert defaults.apply_signature_merges is False
    assert defaults.apply_c6_review is False
    assert defaults.apply_resolution_calibration is False


@pytest.mark.parametrize(
    "kind, stem",
    [
        ("materials", "material_parameter_phase_map"),
        ("tpms", "tpms_parameter_phase_map"),
    ],
)
def test_quick_scan_round_trip_on_optional_stack(tmp_path, kind, stem):
    """Exercise real extraction/serialization, not a mocked successful CLI call."""
    for dependency in ("pyvista", "poly2graph", "plotly", "skimage"):
        pytest.importorskip(dependency)
    output = tmp_path / kind
    main(["scan", kind, "--output-dir", str(output)])
    data = load_phase_map(output / (stem + "_records.csv"))
    assert data.summary()["records"] == 9
    # Coarse scans can be unresolved. Failed geometry/projection must remain
    # explicit instead of being replaced by a ball or an abstract polynomial.
    for record in data.records:
        if record.get("error"):
            assert record["source"] == "error"
            assert not record.get("polynomial")
            assert "error" in record["phase_signature"]
    assert data.summary()["missing_grid_cells"] == 0
    assert data.summary()["adaptive_fill_records"] == 0
    assert json.loads((output / "run_plan.json").read_text())["workers"] == 1
    summary = json.loads((output / (stem + "_summary.json")).read_text())
    assert summary["families"][0]["dimension"] == 24
    assert str(output) in str(summary)
    if kind == "tpms":
        assert len(list((output / "geometry").glob("*.json.gz"))) == 9
        assert all(r["surface_is_closed"] == "True" for r in data.records)
        assert summary["scan_parameters"]["domain_kind"] == "sphere"
    else:
        assert all(r["classification_computed"] is True for r in data.records)
        assert data.summary()["resolution_calibration_records"] == 0
        from knotted_graph.applications.phase_map_examples import _materials

        records_path = output / (stem + "_records.json")
        before = json.loads(records_path.read_text())
        _materials.main(["--output-dir", str(output), "--reuse-records"])
        assert json.loads(records_path.read_text()) == before
        reused = json.loads((output / (stem + "_summary.json")).read_text())
        assert (
            reused["map_info"][data.family][
                "classification_resolution_calibration_cells"
            ]
            == 0
        )
        assert reused["map_info"][data.family]["c6_symmetry_merges"] == []
        assert reused["map_info"]["reprocessed_from_existing_records"] is True


def test_c6_display_grouping_requires_explicit_research_option():
    for dependency in ("pyvista", "poly2graph", "plotly", "skimage"):
        pytest.importorskip(dependency)
    import numpy as np
    from knotted_graph.applications.phase_map_examples import _materials

    invalid = sorted(_materials.C6_INVALID_PHASE_SIGNATURES["tib2_d6_F"])[0]
    signatures = {"yamada:Integer(-1)": 1, invalid: 2}
    grid = np.array([[1, 1], [2, 2]])
    raw, manual, c6 = _materials._reviewed_display_merges(
        grid.copy(), signatures, "tib2_d6_F", _materials.parse_args([])
    )
    assert np.array_equal(raw, grid)
    assert manual == c6 == []
    displayed, _, c6 = _materials._reviewed_display_merges(
        grid.copy(),
        signatures,
        "tib2_d6_F",
        _materials.parse_args(["--apply-c6-review"]),
    )
    assert np.all(displayed == 1)
    assert sum(item["cells"] for item in c6) == 2
    assert np.array_equal(grid, [[1, 1], [2, 2]])


@pytest.mark.parametrize("adaptive_step", [0.0, 0.4])
def test_energy_sampling_preserves_results_until_explicit_calibration(
    monkeypatch, adaptive_step
):
    for dependency in ("pyvista", "poly2graph", "plotly", "skimage"):
        pytest.importorskip(dependency)
    from dataclasses import dataclass, replace
    from knotted_graph.applications.phase_map_examples import _materials

    @dataclass(frozen=True)
    class Reading:
        energy: float
        phase_signature: str
        lam: float = 0.0
        classification_computed: bool = True
        resolution_calibration_energy: float | None = None
        boundary_components: int = 1
        error: str | None = None
        cycle_rank: int = 1
        removed_component_voxels: int = 0
        filled_void_voxels: int = 0

    family = replace(
        _materials.material_families(8)[0],
        energies=(0.2, 0.5, 0.75, 1.0),
        landmark_energies=(),
    )
    assert family.key == "tib2_d6_F"

    def classify(family, obj, lam, energy, **kwargs):
        return Reading(energy, "thin" if energy < 0.5 else "resolved")

    monkeypatch.setattr(_materials, "evaluate_cell", classify)
    records = _materials.evaluate_energy_grid_adaptive(
        family,
        None,
        0.0,
        max_exact_yamada_edges=3,
        adaptive_energy_step=adaptive_step,
    )
    assert records[0].phase_signature == "thin"
    assert all(r.resolution_calibration_energy is None for r in records)
    assert all(r.classification_computed for r in records)
    calibrated = _materials.stabilize_tib2_records(records)
    assert calibrated[0].phase_signature == "resolved"
    assert calibrated[0].classification_computed is False
    assert calibrated[0].resolution_calibration_energy == 0.75
    assert records[0].phase_signature == "thin"


def test_tpms_viewer_reconstructs_exact_recorded_domain():
    for dependency in ("pyvista", "poly2graph", "skimage"):
        pytest.importorskip(dependency)
    from knotted_graph.applications.phase_map_examples import _tpms

    script = (
        Path(__file__).resolve().parents[2]
        / "User_guide/applications/NewPhaseMapPlots"
        / "TPMS/scripts/tpms_plotly_region_geometry.py"
    )
    spec = importlib.util.spec_from_file_location("tpms_viewer", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    options = {
        "domain_kind": "cylinder",
        "span_half_width": 5.0,
        "domain_radius_fraction": 0.6,
    }
    family = _tpms.tpms_families(dimension=12, thresholds=(0.0, 0.3), **options)[0]
    saved = {
        "dimension": 12,
        "thresholds": [0.0, 0.3],
        "key": family.key,
        "span": family.span,
        "compact_domain": {"key": family.domain.key, "formula": family.domain.formula},
    }
    restored = module.family_from_summary(_tpms, saved, options)
    assert restored.domain.formula == family.domain.formula
    assert restored.span == family.span
    with pytest.raises(ValueError, match="differs"):
        module.family_from_summary(_tpms, saved, {})


def test_demo_separation_is_lossless_and_hash_guarded(tmp_path):
    import hashlib

    script = Path(__file__).resolve().parents[2] / "dev/build_phase_map_demos.py"
    spec = importlib.util.spec_from_file_location("build_phase_map_demos", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    payload = {
        "transitions": [{"key": "test"}],
        "regions": {
            "classic:test:1": {"polynomial": "Y+1", "surface": {"x": [1, 2, 3]}},
            "classic:test:2": {"polynomial": "Y+1", "surface": {"x": [1, 2, 3]}},
        },
    }
    text = (
        "<html><body><script>const payload = " + json.dumps(payload) + ";\n"
        "const regions = payload.regions;\n"
        "function selectRegion(regionKey) {\n  const region = regions[regionKey];\n}\n"
        'Plotly.react("phaseMap", [trace, boundaryTrace], layout, { responsive: true, displaylogo: false });\n'
        'Plotly.react("geometry", geometryTraces(region), layout, { responsive: true, displaylogo: false });\n'
        "buildButtons();\nplotPhaseMap(activeTransition, activeMode);</script></body></html>"
    )
    source = tmp_path / "original.html"
    source.write_text(text)
    before = source.read_bytes()
    with pytest.raises(ValueError, match="Source changed"):
        module.split_demo(source, tmp_path / "wrong", expected_hash="wrong")
    dest = tmp_path / "demo"
    manifest = module.split_demo(
        source, dest, expected_hash=hashlib.sha256(before).hexdigest()
    )
    rebuilt, _, _ = module.extract_payload((dest / "index.html").read_text())
    for key, entry in rebuilt["regions"].items():
        rebuilt["regions"][key] = json.loads((dest / entry["asset"]).read_text())
    assert rebuilt == payload
    assert len(list((dest / "regions").glob("*.json"))) == 1
    assert len(manifest["regions"]) == 2
    assert source.read_bytes() == before
    assert "request !== regionSelection" in (dest / "index.html").read_text()
