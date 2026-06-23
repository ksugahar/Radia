"""Curved-boundary mesh geometry gates.

These tests quantify the linear chord error that high-order curved exports are
supposed to remove.  They are analytical and do not require a Cubit session.
"""

import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.mesh_geometry import (
    circle_polyline_area_metrics,
    circular_arc_chord_area_deficit,
    regular_polygon_circle_area,
    regular_polygon_circle_area_metrics,
    required_regular_polygon_sides_for_area_error,
)


def test_arc_chord_deficit_matches_sector_minus_triangle():
    radius = 2.0
    angle = math.pi / 3.0
    got = circular_arc_chord_area_deficit(radius, angle)
    expected = 0.5 * radius * radius * (angle - math.sin(angle))
    assert math.isclose(got, expected, rel_tol=1e-15)

    metrics = circle_polyline_area_metrics(radius, [angle] * 6)
    assert math.isclose(metrics["exact_area"], math.pi * radius * radius, rel_tol=1e-15)
    assert math.isclose(metrics["faceted_area"], 6.0 * 0.5 * radius * radius * math.sin(angle),
                        rel_tol=1e-15)
    assert math.isclose(metrics["area_deficit"], 6.0 * got, rel_tol=1e-15)


def test_regular_polygon_area_known_cases():
    assert math.isclose(regular_polygon_circle_area(1.0, 4), 2.0, rel_tol=1e-15)

    metrics = regular_polygon_circle_area_metrics(1.0, 6)
    assert metrics["n_sides"] == 6
    assert math.isclose(metrics["faceted_area"], 3.0 * math.sqrt(3.0) / 2.0, rel_tol=1e-15)
    assert metrics["relative_deficit"] > 0.0


def test_fine_polygon_matches_asymptotic_area_deficit():
    n = 96
    metrics = regular_polygon_circle_area_metrics(1.0, n)
    asymptotic = 2.0 * math.pi * math.pi / (3.0 * n * n)
    assert math.isclose(metrics["relative_deficit"], asymptotic, rel_tol=2e-3)


def test_required_side_count_is_minimal():
    target = 1.0e-3
    n = required_regular_polygon_sides_for_area_error(target)
    assert regular_polygon_circle_area_metrics(1.0, n)["relative_deficit"] <= target
    assert regular_polygon_circle_area_metrics(1.0, n - 1)["relative_deficit"] > target


def test_invalid_inputs():
    with pytest.raises(ValueError):
        circular_arc_chord_area_deficit(0.0, 0.1)
    with pytest.raises(ValueError):
        circular_arc_chord_area_deficit(1.0, 0.0)
    with pytest.raises(ValueError):
        circle_polyline_area_metrics(1.0, [])
    with pytest.raises(ValueError):
        circle_polyline_area_metrics(1.0, [math.pi, math.pi, 0.1])
    with pytest.raises(ValueError):
        regular_polygon_circle_area(1.0, 2)
    with pytest.raises(ValueError):
        required_regular_polygon_sides_for_area_error(1.0)
