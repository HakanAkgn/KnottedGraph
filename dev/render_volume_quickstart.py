"""Render the input volume, extracted spatial graph and computed projection."""

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
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage.measure import marching_cubes

    example = runpy.run_path(str(ROOT / "examples" / "volume_to_yamada.py"))
    volume, graph, result = example["compute_volume_quickstart"]()
    projection = result.projection
    vertices, faces, _, _ = marching_cubes(volume.astype(np.float32), level=0.5)
    color, node_color = "#236f91", "#b83e52"
    low, high = vertices.min(axis=0), vertices.max(axis=0)
    center = (low + high) / 2
    radius = (high - low).max() / 2 + 3

    with plt.rc_context({"font.size": 12, "font.family": "DejaVu Sans"}):
        figure = plt.figure(figsize=(13.5, 4.6), facecolor="white")
        left = figure.add_axes((0.005, 0.10, 0.30, 0.78), projection="3d")
        middle = figure.add_axes((0.35, 0.10, 0.30, 0.78), projection="3d")
        right = figure.add_axes((0.71, 0.14, 0.27, 0.69))
        surface = Poly3DCollection(vertices[faces], facecolors="#9cc7cf",
                                   edgecolors="#9cc7cf", linewidths=0, alpha=0.85, shade=True,
                                   lightsource=matplotlib.colors.LightSource(azdeg=315, altdeg=35))
        left.add_collection3d(surface)
        for _, _, data in graph.edges(data=True):
            points = np.asarray(data["pts"])
            middle.plot(*points.T, color=color, linewidth=2.7)
        for _, data in graph.nodes(data=True):
            middle.scatter(*data["pos"], color=node_color, s=48, depthshade=False)
        for axis in (left, middle):
            axis.set(xlim=(center[0] - radius, center[0] + radius),
                     ylim=(center[1] - radius, center[1] + radius),
                     zlim=(center[2] - radius, center[2] + radius))
            axis.set_box_aspect((1, 1, 1), zoom=1.4)
            axis.set_proj_type("ortho")
            axis.view_init(elev=25, azim=-62)
            axis.set_axis_off()
        left.set_title("1. Surface of the input volume", fontweight="bold", pad=12)
        middle.set_title("2. Extracted spatial graph", fontweight="bold", pad=12)

        right.set_aspect("equal")
        for arc in projection.arcs:
            points = np.asarray(arc.line.coords)
            right.plot(points[:, 0], points[:, 1], color=color, linewidth=2.7)
        bounds = np.concatenate([np.asarray(a.line.coords)[:, :2] for a in projection.arcs])
        gap = np.ptp(bounds, axis=0).max() * 0.016
        arc_by_id = {arc.id: arc for arc in projection.arcs}
        for crossing in projection.crossings:
            endpoints = []
            for arc_id, _ in crossing.incident_arcs:
                arc = arc_by_id[arc_id]
                at_start = arc.start_type == "x" and arc.start_id == crossing.id
                z = arc.line.coords[0 if at_start else -1][2]
                endpoints.append((z, at_start, arc.line))
            top_z = max(z for z, _, _ in endpoints)
            right.add_patch(Circle((crossing.point.x, crossing.point.y), gap,
                                   facecolor="white", edgecolor="none", zorder=3))
            for z, at_start, line in endpoints:
                if np.isclose(z, top_z):
                    start, end = ((0, 2.3 * gap) if at_start
                                  else (max(0, line.length - 2.3 * gap), line.length))
                    stub = np.asarray(substring(line, start, end).coords)
                    right.plot(stub[:, 0], stub[:, 1], color=color, linewidth=2.7, zorder=4)
        for vertex in projection.vertices:
            right.scatter(vertex.point.x, vertex.point.y, color=node_color, s=48, zorder=5)
        for arc in projection.arcs:
            midpoint = arc.line.interpolate(0.5, normalized=True)
            offset = {3: (11, -10), 4: (5, -13), 5: (6, 7), 7: (-13, 5)}.get(arc.id, (6, 6))
            right.annotate(str(arc.id), (midpoint.x, midpoint.y), xytext=offset,
                           textcoords="offset points", fontsize=9, color="#33434d")
        right.margins(0.12)
        right.set_axis_off()
        right.set_title("3. Projection and PD encoding", fontweight="bold", pad=18)
        for x in (0.325, 0.68):
            figure.text(x, 0.51, "→", fontsize=28, ha="center", color="#55616c")
        figure.text(0.16, 0.045, "96 × 96 × 96 voxels", ha="center", color="#55616c")
        figure.text(0.50, 0.045, "2 vertices · 3 edges", ha="center", color="#55616c")
        figure.text(0.85, 0.045, "4 crossings · arc labels match the PD code", ha="center",
                    color="#55616c", fontsize=10)
        destination = ROOT / "doc" / "assets" / "site_figures" / "quickstart-volume.png"
        figure.savefig(destination, dpi=180, facecolor="white")
        plt.close(figure)
        print(destination)


if __name__ == "__main__":
    main()
