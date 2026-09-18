"""Runtime helpers without process-wide environment or import-path changes."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def package_location() -> str:
    """Record the installed package location, not a developer's checkout path."""
    return str(Path(__file__).resolve().parents[2])


def geometry_directory(scan_dir: Path) -> Path:
    """Support both fresh scan/geometry and archived data/../geometry layouts."""
    direct = scan_dir / "geometry"
    if direct.is_dir() or scan_dir.name != "data":
        return direct
    return scan_dir.parent / "geometry"


def resolve_geometry_path(path_text: str, scan_dir: Path) -> Path:
    """Relocate a legacy record inside the selected dataset; reject missing files."""
    name = Path(path_text.replace("\\", "/")).name
    candidate = geometry_directory(scan_dir) / name
    if not candidate.is_file():
        raise FileNotFoundError(f"Missing geometry {name!r} beside dataset {scan_dir}")
    return candidate.resolve()


def require_scan_dependencies() -> None:
    """Fail before a scan starts with one actionable optional-install message."""
    missing = [
        name
        for name in ("pyvista", "skimage", "plotly", "poly2graph")
        if importlib.util.find_spec(name) is None
    ]
    if missing:
        raise ImportError(
            "Phase-map scans need optional dependencies: "
            + ", ".join(missing)
            + ". From this development checkout run `uv sync --extra nodal --extra viz`, "
            'or `pip install ".[nodal,viz]"`. The legacy PyPI release is not sufficient. '
            "The inspect and plot commands only need the base package."
        )
