"""Separate reading/plotting saved results from explicitly requested computation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .records import PhaseMapData, load_phase_map, read_phase_map_records

FAMILIES = {
    "materials": ("tib2_d6_F", "co2mnga_t8", "ti3al_M2", "yh3_m1"),
    "tpms": ("gyroid_to_diamond", "gyroid_to_schwarz_p", "schwarz_p_to_diamond"),
}


def _positive(value: str) -> int:
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser(
        "inspect", help="summarize saved CSV/JSON; no extraction"
    )
    inspect.add_argument("records", type=Path)
    plot = commands.add_parser("plot", help="plot raw saved signatures; no extraction")
    plot.add_argument("records", type=Path)
    plot.add_argument("--family", help="family key printed by inspect")
    plot.add_argument(
        "--output", type=Path, required=True, help="filename prefix for PNG/PDF/JSON"
    )
    scan = commands.add_parser(
        "scan", help="compute a new optional-dependency scan (use a compute node)"
    )
    scan.add_argument("kind", choices=FAMILIES)
    scan.add_argument(
        "--profile",
        choices=("quick", "paper"),
        default="quick",
        help="quick is a coarse smoke example, not scientifically converged",
    )
    scan.add_argument(
        "--family",
        action="append",
        help="repeat to select families; default: first family only",
    )
    scan.add_argument("--dimension", type=_positive)
    scan.add_argument("--lambda-count", type=_positive)
    scan.add_argument(
        "--level-count",
        type=_positive,
        help="energy (materials) or threshold (TPMS) sample count",
    )
    scan.add_argument(
        "--workers",
        type=_positive,
        default=1,
        help="material processes; default 1; TPMS is serial",
    )
    scan.add_argument("--max-exact-yamada-edges", type=_positive)
    scan.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="new/empty output directory; never the checked-in reference data",
    )
    scan.add_argument(
        "--dry-run",
        action="store_true",
        help="show plan without imports, output files or computation",
    )
    return result


def scan_plan(args: argparse.Namespace) -> dict:
    """Resolve bounded defaults without importing an optional scientific engine."""
    quick = args.profile == "quick"
    families = args.family or [FAMILIES[args.kind][0]]
    unknown = set(families) - set(FAMILIES[args.kind])
    if unknown:
        raise ValueError(
            f"Unknown {args.kind} families {sorted(unknown)}; choose from {FAMILIES[args.kind]}."
        )
    dimension = args.dimension or (
        24 if quick else (140 if args.kind == "materials" else 64)
    )
    lambdas = args.lambda_count or (
        3 if quick else (60 if args.kind == "materials" else 21)
    )
    levels = args.level_count or (
        3 if quick else (0 if args.kind == "materials" else 21)
    )
    if dimension < 8 or lambdas < 2 or levels == 1:
        raise ValueError("Use dimension >= 8 and at least 2 samples on each scan axis.")
    if args.kind == "tpms" and args.workers != 1:
        raise ValueError("TPMS currently runs serially; use --workers 1.")
    return {
        "kind": args.kind,
        "profile": args.profile,
        "families": list(dict.fromkeys(families)),
        "dimension": dimension,
        "lambda_count": lambdas,
        "level_count": levels,
        "workers": args.workers,
        "max_exact_yamada_edges": args.max_exact_yamada_edges or (8 if quick else 18),
        "output_dir": str(args.output_dir),
        "processing": "all requested cells classified; no resolution calibration, island smoothing, manual signature merges or C6 display grouping",
        "scope": "finite-grid application example; not a convergence result or exact phase-boundary proof",
    }


def main(argv: list[str] | None = None) -> None:
    """Run the selected application command, reporting input errors without a traceback."""
    argparser = parser()
    args = argparser.parse_args(argv)
    try:
        if args.command == "inspect":
            records = read_phase_map_records(args.records)
            families = sorted({r.get("material", r.get("family")) for r in records})
            print(
                json.dumps(
                    [PhaseMapData.from_records(records, f).summary() for f in families],
                    indent=2,
                )
            )
        elif args.command == "plot":
            from .plotting import plot_phase_map
            import matplotlib.pyplot as plt

            data = load_phase_map(args.records, family=args.family)
            fig = plot_phase_map(data, output=args.output)
            plt.close(fig)
            print(json.dumps(data.summary(), indent=2))
        else:
            plan = scan_plan(args)
            print(json.dumps(plan, indent=2), flush=True)
            if args.dry_run:
                return
            if args.output_dir.exists() and (
                not args.output_dir.is_dir() or any(args.output_dir.iterdir())
            ):
                raise ValueError(
                    "Choose a new or empty --output-dir; existing results will not be overwritten."
                )
            from ._runtime import require_scan_dependencies

            require_scan_dependencies()
            import matplotlib.pyplot as plt

            common = [
                "--output-dir",
                str(args.output_dir),
                "--dimension",
                str(plan["dimension"]),
                "--lambda-count",
                str(plan["lambda_count"]),
                "--max-exact-yamada-edges",
                str(plan["max_exact_yamada_edges"]),
                "--min-stable-cells",
                "1",
                "--only",
                *plan["families"],
            ]
            with plt.rc_context(
                {"font.family": "serif", "mathtext.fontset": "cm", "pdf.fonttype": 42}
            ):
                if args.kind == "materials":
                    from . import _materials

                    _materials.main(
                        common
                        + [
                            "--energy-count",
                            str(plan["level_count"]),
                            "--workers",
                            str(plan["workers"]),
                            "--adaptive-energy-step",
                            "0",
                        ]
                    )
                else:
                    from . import _tpms

                    _tpms.main(
                        common
                        + [
                            "--threshold-count",
                            str(plan["level_count"]),
                            "--threshold-min",
                            "0",
                            "--threshold-max",
                            "0.3",
                        ]
                    )
            (args.output_dir / "run_plan.json").write_text(
                json.dumps(plan, indent=2) + "\n", encoding="utf-8"
            )
    except (OSError, ValueError, ImportError) as exc:
        argparser.error(str(exc))
