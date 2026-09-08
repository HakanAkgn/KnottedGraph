"""Raw categorical phase maps with explicit provenance and Computer Modern math."""

from __future__ import annotations

import json
import colorsys
from pathlib import Path
from typing import Any

from .records import PhaseMapData


def plot_phase_map(data: PhaseMapData, *, output: str | Path | None = None) -> Any:
    """Plot one raw saved grid and optionally save PNG, PDF and a phase-key JSON.

    ``output`` is a filename prefix (not a directory). Missing cells stay white,
    adaptive fills are shaded, and errors are marked with crosses. Category
    numbers only identify saved signatures; they are not numerical Yamada
    values and do not establish equivalence under contraction.

    Returns a Matplotlib Figure; the caller owns/should close it. Matplotlib
    configuration is scoped to this call, not applied globally on import.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import BoundaryNorm, ListedColormap, to_hex

    if not data.records:
        raise ValueError("Cannot plot an empty phase map.")
    xs = sorted({r["lam"] for r in data.records})
    ys = sorted({r[data.level_field] for r in data.records})
    xi, yi = {v: i for i, v in enumerate(xs)}, {v: i for i, v in enumerate(ys)}
    signatures = sorted({r["phase_signature"] for r in data.records})
    ids = {signature: index for index, signature in enumerate(signatures, 1)}
    grid = np.full((len(ys), len(xs)), np.nan)
    inferred = np.full_like(grid, np.nan)
    errors = []
    for record in data.records:
        x, y = record["lam"], record[data.level_field]
        grid[yi[y], xi[x]] = ids[record["phase_signature"]]
        if record.get("classification_computed") is False:
            inferred[yi[y], xi[x]] = 1.0
        if record.get("error") or record["source"] == "error":
            errors.append((x, y))

    style = {
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "text.usetex": False,
        "font.size": 11,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 1.1,
    }
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(7.4, 5.8))
        colors = (
            [plt.get_cmap("tab20")(i) for i in range(len(signatures))]
            if len(signatures) <= 20
            else [
                colorsys.hsv_to_rgb((i * 0.61803398875) % 1, 0.62, 0.82)
                for i in range(len(signatures))
            ]
        )
        cmap = ListedColormap(colors).with_extremes(bad="white")
        norm = BoundaryNorm(np.arange(0.5, len(signatures) + 1.5), cmap.N)
        mesh = ax.pcolormesh(
            xs, ys, grid, cmap=cmap, norm=norm, shading="nearest", rasterized=True
        )
        ax.pcolormesh(
            xs,
            ys,
            inferred,
            cmap=ListedColormap(["#111111"]),
            alpha=0.28,
            shading="nearest",
            rasterized=True,
        )
        if errors:
            ex, ey = zip(*errors)
            ax.scatter(ex, ey, marker="x", c="black", s=12, linewidths=0.8)
        ax.set_xlabel(r"$\lambda$")
        ax.set_ylabel(r"$E$" if data.level_field == "energy" else r"$c$")
        title = data.records[0].get("title") or data.family.replace("_", " ")
        ax.set_title(f"{title}: raw phase signatures", pad=12)
        if len(signatures) <= 20:
            bar = fig.colorbar(mesh, ax=ax, ticks=list(ids.values()), pad=0.03)
            bar.set_label("Saved signature ID (not a polynomial value)")
        summary = data.summary()
        fig.text(
            0.5,
            0.04,
            "Raw signatures • white: missing • shading: adaptive fill • ×: error\n"
            f"{len(signatures)} categories; phase-key JSON lists full signatures and sources.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.09, 1, 1))
        if output is not None:
            prefix = Path(output)
            if prefix.suffix.lower() in (".png", ".pdf", ".json"):
                prefix = prefix.with_suffix("")
            prefix.parent.mkdir(parents=True, exist_ok=True)
            for suffix in (".png", ".pdf"):
                fig.savefig(str(prefix) + suffix, dpi=200, facecolor="white")
            key = {
                str(ids[s]): {
                    "signature": s,
                    "color": to_hex(colors[ids[s] - 1]),
                    "sources": sorted(
                        {r["source"] for r in data.records if r["phase_signature"] == s}
                    ),
                    "polynomials": sorted(
                        {
                            str(r.get("polynomial", ""))
                            for r in data.records
                            if r["phase_signature"] == s
                        }
                    ),
                }
                for s in signatures
            }
            Path(str(prefix) + ".json").write_text(
                json.dumps(
                    {
                        "summary": summary,
                        "phase_key": key,
                        "axes": {"lambda": xs, data.level_field: ys},
                        "grid_order": f"rows={data.level_field}, columns=lambda; null=missing",
                        "signature_ids": [
                            [None if np.isnan(v) else int(v) for v in row]
                            for row in grid
                        ],
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
    return fig
