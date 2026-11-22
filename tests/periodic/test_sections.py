import numpy as np
import pytest

from knotted_graph.periodic import periodic_plane_contours


def test_plane_winding_along_oblique_lattice_direction():
    t = (np.arange(56) + .3) / 56
    values = np.sin(2 * np.pi * (t[:, None] + 2 * t[None, :]) + .19)
    for field in (values, np.roll(values, (17, 9), axis=(0, 1))):
        contours = periodic_plane_contours(field, 0)
        assert len(contours) == 2
        for orbit in contours:
            assert orbit.winding in ((2, -1), (-2, 1))
            assert orbit.area is None
            assert np.allclose(orbit.points[-1] - orbit.points[0], orbit.winding)


def test_plane_closed_area_across_axis_cut():
    t = (np.arange(96) + .34) / 96
    edge_distance = np.minimum(t, 1 - t)
    values = edge_distance[:, None]**2 + edge_distance[None, :]**2
    radius = .17
    expected = np.pi * radius**2
    first, = periodic_plane_contours(values, radius**2)
    moved, = periodic_plane_contours(np.roll(values, (-13, 22), axis=(0, 1)), radius**2)
    assert first.winding == moved.winding == (0, 0)
    assert first.area == pytest.approx(moved.area, abs=1e-13)
    assert first.area == pytest.approx(expected, abs=1.4e-4)


def test_no_plane_contours_outside_energy_range():
    assert not periodic_plane_contours(np.zeros((3, 4)), 1)
