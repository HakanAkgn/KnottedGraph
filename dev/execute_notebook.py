from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("notebook")
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Save the executed notebook, HTML and execution status under this directory.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    notebook = (root / args.notebook).resolve()
    if root not in notebook.parents:
        raise ValueError(f"Notebook must be inside repository: {notebook}")
    if not notebook.exists():
        raise FileNotFoundError(notebook)
    output = None
    if args.output_dir is not None:
        from nbconvert import HTMLExporter

        output = args.output_dir.resolve() / notebook.relative_to(root)
        if output == notebook:
            raise ValueError("Execution output must not overwrite the source notebook.")
        for target in (
            output,
            output.with_suffix(".html"),
            output.with_suffix(".status.json"),
        ):
            if target.exists():
                raise FileExistsError(f"Choose a fresh output directory: {target}")
        output.parent.mkdir(parents=True, exist_ok=True)

    # Execute with the active environment exactly as a user would.  The CI and
    # PBS gates install this checkout before calling the helper; injecting the
    # raw ``src`` directory here would bypass editable-install hooks and could
    # hide packaging or native-extension failures.
    os.environ.pop("PYTHONPATH", None)
    os.environ["PYTHONNOUSERSITE"] = "1"
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
    os.environ.setdefault("PLOTLY_RENDERER", "json")

    with notebook.open("r", encoding="utf-8") as handle:
        nb = nbformat.read(handle, as_version=4)
    # A failed rerun must not retain old output from cells it never reached.
    for cell in nb.cells:
        if cell.cell_type == "code":
            cell["execution_count"] = None
            cell["outputs"] = []

    client = NotebookClient(
        nb,
        timeout=args.timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(root)}},
        allow_errors=False,
    )
    started = time.monotonic()
    execution_error = None
    try:
        client.execute(cwd=str(root))
    except Exception as exc:
        execution_error = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        if output is not None:
            nbformat.write(nb, output)
            status = {
                "source": str(notebook.relative_to(root)),
                "execution": "passed" if execution_error is None else "failed",
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "error": execution_error,
                "executed_code_cells": sum(
                    cell.cell_type == "code" and cell.get("execution_count") is not None
                    for cell in nb.cells
                ),
                "total_code_cells": sum(cell.cell_type == "code" for cell in nb.cells),
                "visual_review": "pending",
                "executed_notebook": str(output),
            }
            try:
                html, _ = HTMLExporter(template_name="lab").from_notebook_node(nb)
                output.with_suffix(".html").write_text(html, encoding="utf-8")
                status["html"] = str(output.with_suffix(".html"))
                status["html_export"] = "passed"
            except Exception as exc:
                status["html_export"] = "failed"
                status["html_error"] = f"{type(exc).__name__}: {exc}"
                if execution_error is None:
                    raise
            finally:
                output.with_suffix(".status.json").write_text(
                    json.dumps(status, indent=2) + "\n", encoding="utf-8"
                )
    print(f"PASS {notebook.relative_to(root)}")


if __name__ == "__main__":
    main()
