from __future__ import annotations

import numpy as np
import pytest
from shapely import LineString, Point

from knotted_graph.projection.geom import _projected_endpoint_angle


@pytest.mark.parametrize("scale", [1e-12, 1.0, 1e12])
def test_projected_endpoint_angle_skips_machine_negligible_start_step(scale):
    line = LineString(
        [
            (0.0, 0.0, 2.0 * scale),
            (1e-17 * scale, 0.0, 2.0 * scale),
            (0.0, -1.0 * scale, 2.0 * scale),
        ]
    )
    angle = _projected_endpoint_angle(
        line,
        Point(0.0, 0.0, 2.0 * scale),
        start=True,
    )
    assert angle == pytest.approx(-np.pi / 2.0)


@pytest.mark.parametrize("scale", [1e-12, 1.0, 1e12])
def test_projected_endpoint_angle_skips_machine_negligible_end_step(scale):
    line = LineString(
        [
            (0.0, -1.0 * scale, 2.0 * scale),
            (1e-17 * scale, 0.0, 2.0 * scale),
            (0.0, 0.0, 2.0 * scale),
        ]
    )
    angle = _projected_endpoint_angle(
        line,
        Point(0.0, 0.0, 2.0 * scale),
        start=False,
    )
    assert angle == pytest.approx(-np.pi / 2.0)


def test_projected_endpoint_angle_keeps_resolvable_small_step():
    line = LineString(
        [
            (0.0, 0.0, 2.0),
            (1e-12, 0.0, 2.0),
            (0.0, -1.0, 2.0),
        ]
    )
    angle = _projected_endpoint_angle(
        line,
        Point(0.0, 0.0, 2.0),
        start=True,
    )
    assert angle == pytest.approx(0.0)


def test_projected_endpoint_angle_rejects_zero_projected_extent():
    line = LineString(
        [
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 2.0),
        ]
    )
    with pytest.raises(ValueError, match="no resolvable projected endpoint direction"):
        _projected_endpoint_angle(
            line,
            Point(0.0, 0.0, 0.0),
            start=True,
        )
