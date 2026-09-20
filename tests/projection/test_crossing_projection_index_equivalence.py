from __future__ import annotations

import numpy as np
import pytest
from shapely import LineString, MultiLineString, Point
from shapely.strtree import STRtree

from knotted_graph.projection.pd_code import PDCode


def test_indexed_crossing_projection_matches_legacy_scan_randomized():
    rng = np.random.default_rng(20260817)

    for _ in range(100):
        count = int(rng.integers(2, 18))
        xyz = np.cumsum(rng.normal(size=(count, 3)), axis=0)
        edge = LineString(xyz.tolist())

        crossings = []
        # Generate points exactly on randomly selected segments, then append
        # unrelated points. This exercises ordering, duplicate suppression, and
        # the spatial-index candidate filter against the legacy O(S*C) oracle.
        for _ in range(int(rng.integers(0, 20))):
            segment_index = int(rng.integers(0, count - 1))
            t = float(rng.uniform(0.05, 0.95))
            a = xyz[segment_index]
            b = xyz[segment_index + 1]
            point = a + t * (b - a)
            crossings.append(Point(point))

        for _ in range(int(rng.integers(0, 20))):
            crossings.append(Point(rng.normal(size=3) * 20.0))

        expected = PDCode._project_crossings_on_edge(edge, crossings, tolerance=1e-8)
        if crossings:
            tree = STRtree(crossings)
            actual = PDCode._project_crossings_on_edge_indexed(
                edge,
                crossings,
                tree,
                tolerance=1e-8,
            )
        else:
            actual = []

        assert actual == expected


def test_indexed_crossing_projection_preserves_tolerance_boundary_behavior():
    edge = LineString([(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)])
    tolerance = 1e-8
    crossings = [
        Point(2.0, 0.0, 0.0),
        Point(4.0, 0.5 * tolerance, 0.0),
        Point(6.0, 0.999 * tolerance, 0.0),
        Point(8.0, 1.001 * tolerance, 0.0),
    ]
    tree = STRtree(crossings)

    expected = PDCode._project_crossings_on_edge(edge, crossings, tolerance=tolerance)
    actual = PDCode._project_crossings_on_edge_indexed(
        edge,
        crossings,
        tree,
        tolerance=tolerance,
    )

    assert actual == expected



def test_fused_crossing_incidence_matches_independent_legacy_pipeline_randomized():
    rng = np.random.default_rng(20260920)

    for _ in range(100):
        lines = []
        for edge_index in range(int(rng.integers(2, 7))):
            count = int(rng.integers(2, 15))
            xyz = np.cumsum(rng.normal(size=(count, 3)), axis=0)
            xyz[:, 2] += 0.17 * edge_index
            lines.append(LineString(xyz))

        geometry = MultiLineString(lines)
        crossing_points, fused = PDCode._find_crossings_with_incidences(
            geometry,
            tolerance=1e-8,
        )
        legacy_points = PDCode._find_all_crossings(
            geometry,
            tolerance=1e-8,
        )

        assert len(crossing_points) == len(legacy_points)
        for actual, expected in zip(crossing_points, legacy_points, strict=True):
            assert actual.distance(expected) < 1e-9

        if legacy_points:
            tree = STRtree(legacy_points)
            legacy = [
                PDCode._project_crossings_on_edge_indexed(
                    edge,
                    legacy_points,
                    tree,
                    tolerance=1e-8,
                )
                for edge in geometry.geoms
            ]
        else:
            legacy = [[] for _ in geometry.geoms]

        assert len(fused) == len(legacy)
        for actual, expected in zip(fused, legacy, strict=True):
            assert len(actual) == len(expected)
            for (actual_distance, actual_id), (expected_distance, expected_id) in zip(
                actual,
                expected,
                strict=True,
            ):
                assert actual_id == expected_id
                assert actual_distance == pytest.approx(expected_distance, abs=1e-9)
