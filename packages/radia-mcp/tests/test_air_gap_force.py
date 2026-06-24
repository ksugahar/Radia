"""Air-gap Maxwell pressure / holding-force helpers."""

import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    MU0,
    air_gap_force_summary,
    air_gap_holding_force,
    air_gap_maxwell_pressure,
)


def test_air_gap_pressure_matches_maxwell_stress_at_one_tesla():
    expected = 1.0 / (2.0 * MU0)
    assert air_gap_maxwell_pressure(1.0) == pytest.approx(expected)
    assert air_gap_maxwell_pressure(-1.0) == pytest.approx(expected)


def test_air_gap_force_scales_with_b_squared_area_and_faces():
    base = air_gap_holding_force(0.5, area_m2=2.0e-4)
    assert air_gap_holding_force(1.0, area_m2=2.0e-4) == pytest.approx(4.0 * base)
    assert air_gap_holding_force(0.5, area_m2=4.0e-4) == pytest.approx(2.0 * base)
    assert air_gap_holding_force(0.5, area_m2=2.0e-4, faces=2) == pytest.approx(2.0 * base)


def test_air_gap_force_summary_is_json_friendly_and_self_consistent():
    row = air_gap_force_summary(0.8, area_m2=1.5e-4, faces=2)
    pressure = 0.8 * 0.8 / (2.0 * MU0)
    assert row["B_T"] == pytest.approx(0.8)
    assert row["pressure_Pa"] == pytest.approx(pressure)
    assert row["energy_density_J_per_m3"] == pytest.approx(pressure)
    assert row["force_N"] == pytest.approx(pressure * 1.5e-4 * 2)
    assert row["force_per_area_N_per_m2"] == pytest.approx(pressure)


def test_air_gap_force_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        air_gap_maxwell_pressure(1.0, mu=0.0)
    with pytest.raises(ValueError):
        air_gap_holding_force(1.0, area_m2=-1.0)
    with pytest.raises(ValueError):
        air_gap_holding_force(1.0, area_m2=1.0, faces=0)


if __name__ == "__main__":
    test_air_gap_pressure_matches_maxwell_stress_at_one_tesla()
    test_air_gap_force_scales_with_b_squared_area_and_faces()
    test_air_gap_force_summary_is_json_friendly_and_self_consistent()
    test_air_gap_force_rejects_invalid_inputs()
    print("[OK] air-gap Maxwell pressure and holding-force helpers validated.")
