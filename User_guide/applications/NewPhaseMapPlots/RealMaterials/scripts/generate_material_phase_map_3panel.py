import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap


FIGURE_DIR = Path(
    "/Users/hakanakgun/Desktop/Projects/ProfLeeProjects/Knotted_graph_code_paper/"
    "FigureGeneration/figures/07_hamiltonian_yamada_phase_maps"
)
HTML_PATH = FIGURE_DIR / "07_hamiltonian_yamada_plotly_region_geometry.html"
OUTPUT_DIR = Path(__file__).resolve().parent / "material_parameter_phase_maps"
TMP_PDF = OUTPUT_DIR / "07_hamiltonian_yamada_material_phase_maps_3panel_like_porous.pdf"
TMP_PNG = OUTPUT_DIR / "07_hamiltonian_yamada_material_phase_maps_3panel_like_porous.png"
FINAL_PDF = FIGURE_DIR / TMP_PDF.name
FINAL_PNG = FIGURE_DIR / TMP_PNG.name

MATERIAL_KEYS = ("tib2_d6_F", "co2mnga_t8")
PANEL_LABELS = ("(a)", "(b)")
TEXT_COLOR = "#111827"


def configure_paper_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "text.usetex": False,
            "figure.dpi": 220,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.linewidth": 1.5,
            "axes.edgecolor": "black",
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def payload_from_html(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    start = text.index("const payload = ") + len("const payload = ")
    end = text.index("\nconst transitions = payload.transitions;", start)
    raw = text[start:end].strip()
    if raw.endswith(";"):
        raw = raw[:-1]
    return json.loads(raw)


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.22,
        1.055,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=22,
        fontweight="bold",
        family="serif",
        color=TEXT_COLOR,
        clip_on=False,
        zorder=20,
    )


def draw_phase_boundaries(ax, lambdas: np.ndarray, energies: np.ndarray, labels: np.ndarray) -> None:
    for phase_id in sorted(set(int(value) for value in labels.ravel())):
        mask = (labels == phase_id).astype(float)
        if mask.min() == mask.max():
            continue
        ax.contour(
            lambdas,
            energies,
            mask,
            levels=[0.5],
            colors="black",
            linewidths=1.05,
            alpha=1.0,
        )


def panel_title(item: dict) -> str:
    key = item["key"]
    if key == "tib2_d6_F":
        return r"TiB$_2$: $F=(1-\lambda)1.83+\lambda\,7.60$"
    if key == "ti3al_M2":
        return r"Ti$_3$Al: $M_2=(1-\lambda)(-0.52)+\lambda(-0.10)$"
    if key == "co2mnga_t8":
        return r"Co$_2$MnGa: $t_8=(1-\lambda)(-0.34)+\lambda(-0.70)$"
    return str(item.get("title", key))


def style_axis(ax, item: dict, *, show_xlabel: bool = True) -> None:
    energies = np.asarray(item.get("parameters") or item.get("gammas"), dtype=float)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(float(energies[0]), float(energies[-1]))
    ax.set_title(panel_title(item), fontsize=18, fontweight="semibold", pad=11)
    if show_xlabel:
        ax.set_xlabel(r"$\lambda$", fontsize=24, labelpad=7)
    ax.set_ylabel(r"$E$", fontsize=24, labelpad=8)
    ax.tick_params(axis="both", which="major", labelsize=17, width=1.5, length=6, direction="out")
    ax.minorticks_on()
    ax.tick_params(axis="both", which="minor", width=1.0, length=3, direction="out")
    ax.grid(which="major", color="white", linewidth=0.7, alpha=0.22)
    ax.grid(which="minor", color="white", linewidth=0.45, alpha=0.16)
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)
        spine.set_color("black")


def plot_material_panel(ax, item: dict, label: str):
    mode_data = item["modes"]["classic"]
    labels = np.asarray(mode_data["z"], dtype=int)
    lambdas = np.asarray(item["lambdas"], dtype=float)
    energies = np.asarray(item.get("parameters") or item.get("gammas"), dtype=float)
    colors = mode_data["colors"]
    cmap = ListedColormap(colors, name=f"{item['key']}_classic_yamada")
    norm = BoundaryNorm(np.arange(0.5, len(colors) + 1.5, 1.0), cmap.N)
    ax.pcolormesh(
        lambdas,
        energies,
        labels,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        rasterized=False,
    )
    draw_phase_boundaries(ax, lambdas, energies, labels)
    style_axis(ax, item)
    add_panel_label(ax, label)
    return cmap, norm, len(colors)


def add_upsilon_colorbar(fig, ax, cmap, norm, phase_count: int) -> None:
    cax = ax.inset_axes([1.045, 0.0, 0.04, 1.0], transform=ax.transAxes)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(
        sm,
        cax=cax,
        boundaries=np.arange(0.5, phase_count + 1.5, 1.0),
        ticks=[],
        spacing="proportional",
    )
    cbar.ax.set_title(r"$\Upsilon$", fontsize=18, fontweight="semibold", pad=7)
    cbar.ax.tick_params(width=0, length=0)
    cbar.outline.set_linewidth(1.2)


def make_figure(payload: dict) -> None:
    transitions = {item["key"]: item for item in payload["transitions"]}
    materials = [transitions[key] for key in MATERIAL_KEYS]

    fig = plt.figure(figsize=(13.17, 7.60), constrained_layout=False, facecolor="white")
    left = fig.add_axes([0.075, 0.17, 0.34, 0.63])
    right = fig.add_axes([0.57, 0.17, 0.34, 0.63])

    for ax, item, label in zip((left, right), materials, PANEL_LABELS):
        colorbar_data = plot_material_panel(ax, item, label)
        add_upsilon_colorbar(fig, ax, *colorbar_data)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(TMP_PNG, dpi=450, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    fig.savefig(TMP_PDF, dpi=700, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    plt.close(fig)
    print("saved:", TMP_PNG)
    print("saved:", TMP_PDF)
    print("copy targets:")
    print(FINAL_PNG)
    print(FINAL_PDF)


def main() -> None:
    configure_paper_style()
    make_figure(payload_from_html(HTML_PATH))


if __name__ == "__main__":
    main()
