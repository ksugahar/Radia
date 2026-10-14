"""Focused tests for the optional Radia-to-Xsuite beam tracking bridge."""
from __future__ import annotations

import numpy as np
import pytest

xtrack = pytest.importorskip("xtrack")

from radia.xsuite_bridge import (  # noqa: E402
    AxisAlignedBox,
    RadiaBatchFieldMap,
    first_box_exit_events,
    track_magnetic_fieldmap,
)


def test_radia_fieldmap_batches_points_without_unit_conversion():
    calls = []

    def evaluator(handle, field_type, points):
        points = np.asarray(points, dtype=float)
        calls.append((handle, field_type, points.copy()))
        return np.column_stack((points[:, 0] + 1.0, points[:, 1] + 2.0, points[:, 2] + 3.0))

    fieldmap = RadiaBatchFieldMap(17, evaluator=evaluator)
    bx, by, bz = fieldmap(np.array([0.0, 0.1]), 0.2, np.array([0.3, 0.4]))

    assert calls[0][0:2] == (17, "b")
    np.testing.assert_allclose(calls[0][2], [[0.0, 0.2, 0.3], [0.1, 0.2, 0.4]])
    np.testing.assert_allclose(bx, [1.0, 1.1])
    np.testing.assert_allclose(by, [2.2, 2.2])
    np.testing.assert_allclose(bz, [3.3, 3.4])


def test_zero_field_tracks_straight_and_preserves_momentum():
    def zero_field(x, y, z):
        return np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)

    particles = xtrack.Particles(
        "proton",
        p0c=1.0e9,
        x=[0.01],
        px=[1.0e-3],
        y=[-0.02],
        py=[0.0],
    )
    result = track_magnetic_fieldmap(
        particles,
        zero_field,
        s_start_m=0.0,
        s_end_m=1.0,
        n_steps=100,
    )

    assert result["backend"] == "xtrack.BorisSpatialIntegrator"
    assert result["step_count"] == 100
    assert float(particles.x[0]) == pytest.approx(0.0110000005, rel=2.0e-8)
    assert float(particles.y[0]) == pytest.approx(-0.02, abs=1.0e-14)
    np.testing.assert_allclose(
        result["relative_momentum_deviation_after"],
        result["relative_momentum_deviation_before"],
        atol=1.0e-14,
    )


def test_uniform_transverse_field_deflects_without_momentum_drift():
    def uniform_by(x, y, z):
        return np.zeros_like(x), np.full_like(y, 0.25), np.zeros_like(z)

    particles = xtrack.Particles("proton", p0c=1.0e9, x=[0.0], px=[0.0])
    result = track_magnetic_fieldmap(
        particles,
        uniform_by,
        s_start_m=0.0,
        s_end_m=0.5,
        n_steps=500,
    )

    assert abs(float(particles.x[0])) > 1.0e-3
    np.testing.assert_allclose(
        result["relative_momentum_deviation_after"],
        result["relative_momentum_deviation_before"],
        atol=2.0e-12,
    )


def test_first_boundary_exit_is_reported_per_particle():
    x = np.array([[0.0, 0.0], [0.2, 0.0], [0.6, -0.7]])
    y = np.zeros_like(x)
    z = np.array([[0.0, 0.0], [0.2, 0.2], [0.4, 0.4]])
    events = first_box_exit_events(
        x,
        y,
        z,
        AxisAlignedBox(minimum=(-0.5, -0.5, -0.1), maximum=(0.5, 0.5, 1.0)),
    )

    assert events == [
        {"particle_index": 0, "step_index": 2, "position_m": [0.6, 0.0, 0.4]},
        {"particle_index": 1, "step_index": 2, "position_m": [-0.7, 0.0, 0.4]},
    ]
