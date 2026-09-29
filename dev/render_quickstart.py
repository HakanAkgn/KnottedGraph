"""Render the README/website figure from the maintained 3D Quick Start."""

from pathlib import Path
import runpy

import numpy as np
from shapely.ops import substring


ROOT = Path(__file__).resolve().parents[1]


def main():
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.patches import Circle

    example = runpy.run_path(str(ROOT / "examples" / "quickstart.py"))
    graph, result = example["compute_quickstart"]()
    points = graph.edges["loop_anchor", "loop_anchor", "curve"]["pts"]
    projection = result.projection
    color = "#236f91"

    with plt.rc_context({"font.size": 12, "font.family": "DejaVu Sans"}):
        figure = plt.figure(figsize=(11, 4.6), facecolor="white")
        left = figure.add_axes((0.035, 0.09, 0.43, 0.79), projection="3d")
        left.plot(*points.T, color=color, linewidth=3.0)
        left.set_box_aspect((1, 1, 0.65))
        left.view_init(elev=29, azim=-58)
        left.set_xlabel("x", labelpad=3)
        left.set_ylabel("y", labelpad=3)
        left.set_zlabel("z", labelpad=3)
        left.set_xticks([-2, 0, 2])
        left.set_yticks([-2, 0, 2])
        left.set_zticks([-1, 0, 1])
        left.tick_params(labelsize=9, pad=0)
        left.set_title("1. A curve embedded in 3D", pad=16, fontweight="bold")

        right = figure.add_axes((0.58, 0.13, 0.38, 0.72))
        right.set_aspect("equal")
        for arc in projection.arcs:
            xyz = np.asarray(arc.line.coords)
            right.plot(xyz[:, 0], xyz[:, 1], color=color, linewidth=3)

        # Use the computed crossing records and their endpoint depths to draw
        # gaps under the upper strand, preserving the actual over/under order.
        arc_by_id = {arc.id: arc for arc in projection.arcs}
        for crossing in projection.crossings:
            endpoints = []
            for arc_id, _ in crossing.incident_arcs:
                arc = arc_by_id[arc_id]
                at_start = arc.start_type == "x" and arc.start_id == crossing.id
                z = arc.line.coords[0 if at_start else -1][2]
                endpoints.append((z, at_start, arc.line))
            top_z = max(z for z, _, _ in endpoints)
            right.add_patch(Circle(
                (crossing.point.x, crossing.point.y), 0.10,
                facecolor="white", edgecolor="none", zorder=3,
            ))
            for z, at_start, line in endpoints:
                if np.isclose(z, top_z):
                    start, end = ((0, 0.22) if at_start
                                  else (line.length - 0.22, line.length))
                    stub = np.asarray(substring(line, start, end).coords)
                    right.plot(stub[:, 0], stub[:, 1], color=color,
                               linewidth=3, zorder=4)

        right.margins(0.12)
        right.set_axis_off()
        right.set_title("2. Project and resolve crossings", pad=18,
                        fontweight="bold")
        figure.text(0.50, 0.50, "→", fontsize=28, ha="center", color="#55616c")
        figure.text(0.25, 0.025, "120 coordinate samples", ha="center",
                    color="#55616c")
        figure.text(0.77, 0.025,
                    f"{projection.num_crossings} crossings · ready for Yamada evaluation",
                    ha="center", color="#55616c")
        destination = ROOT / "doc" / "assets" / "site_figures" / "quickstart-trefoil.png"
        figure.savefig(destination, dpi=180, facecolor="white")
        plt.close(figure)
        print(destination)


if __name__ == "__main__":
    main()
