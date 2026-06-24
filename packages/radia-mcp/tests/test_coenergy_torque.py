"""Coenergy-derived virtual-work torque helpers."""

import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    coenergy_torque_from_angle_samples,
    coenergy_torque_summary,
)


def test_coenergy_torque_periodic_sinusoid_matches_analytic_derivative():
    samples = 1440
    harmonic = 3
    amplitude = 0.75
    angles = [2.0 * math.pi * index / samples for index in range(samples)]
    coenergy = [2.0 - amplitude * math.cos(harmonic * angle) for angle in angles]
    rows = coenergy_torque_from_angle_samples(angles, coenergy, periodic=True)

    exact = [amplitude * harmonic * math.sin(harmonic * angle) for angle in angles]
    max_abs_error = max(abs(row["torque_Nm"] - ref) for row, ref in zip(rows, exact))
    assert rows[0]["stencil"] == "central_periodic"
    assert max_abs_error < 7.0e-5

    summary = coenergy_torque_summary(angles, coenergy, periodic=True)
    assert summary["n_samples"] == samples
    assert summary["torque_peak_abs_Nm"] == pytest.approx(amplitude * harmonic, rel=4.0e-5)
    assert summary["torque_mean_Nm"] == pytest.approx(0.0, abs=1.0e-13)


def test_coenergy_torque_nonperiodic_linear_table_is_exact():
    angles = [0.0, 0.1, 0.25, 0.4]
    coenergy = [1.0 + 4.0 * angle for angle in angles]
    rows = coenergy_torque_from_angle_samples(angles, coenergy)

    assert [row["stencil"] for row in rows] == ["forward", "central", "central", "backward"]
    assert [row["torque_Nm"] for row in rows] == pytest.approx([4.0, 4.0, 4.0, 4.0])


def test_coenergy_torque_rejects_bad_tables():
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0], [0.0, 1.0])
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0, 2.0], [0.0, 1.0])
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0, 1.0], [0.0, 1.0, 2.0])
    with pytest.raises(ValueError):
        coenergy_torque_from_angle_samples([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], periodic=True, period_rad=0.0)
