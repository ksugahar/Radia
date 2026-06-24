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
    air_gap_shear_stress,
    air_gap_shear_torque,
    air_gap_shear_torque_summary,
    maxwell_contour_force_2d,
    maxwell_line_segment_force_2d,
    maxwell_stress_tensor_air,
    maxwell_traction_air,
    maxwell_traction_summary,
)


def test_air_gap_pressure_matches_maxwell_stress_at_one_tesla():
    expected = 1.0 / (2.0 * MU0)
    assert air_gap_maxwell_pressure(1.0) == pytest.approx(expected)
    assert air_gap_maxwell_pressure(-1.0) == pytest.approx(expected)


def test_maxwell_tensor_normal_field_reduces_to_air_gap_pressure():
    pressure = air_gap_maxwell_pressure(1.0)
    tensor = maxwell_stress_tensor_air((0.0, 0.0, 1.0))
    traction = maxwell_traction_air((0.0, 0.0, 1.0), (0.0, 0.0, 1.0))
    summary = maxwell_traction_summary((0.0, 0.0, 1.0), (0.0, 0.0, 2.0), area_m2=2.0e-4)

    assert tensor[0][0] == pytest.approx(-pressure)
    assert tensor[1][1] == pytest.approx(-pressure)
    assert tensor[2][2] == pytest.approx(pressure)
    assert traction == pytest.approx([0.0, 0.0, pressure])
    assert summary["normal_traction_Pa"] == pytest.approx(pressure)
    assert summary["normal_traction_identity_Pa"] == pytest.approx(pressure)
    assert summary["force_N"] == pytest.approx([0.0, 0.0, pressure * 2.0e-4])


def test_maxwell_tensor_tangential_field_is_magnetic_tension():
    pressure = air_gap_maxwell_pressure(1.0)
    traction = maxwell_traction_air((1.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    summary = maxwell_traction_summary((1.0, 0.0, 0.0), (0.0, 0.0, 1.0))

    assert traction == pytest.approx([0.0, 0.0, -pressure])
    assert summary["B_normal_T"] == pytest.approx(0.0)
    assert summary["B_tangent_T"] == pytest.approx(1.0)
    assert summary["normal_traction_Pa"] == pytest.approx(-pressure)
    assert summary["tangential_traction_magnitude_Pa"] == pytest.approx(0.0)


def test_maxwell_traction_oblique_field_decomposes_into_normal_and_tangent():
    # B=(3,4,0), n=x gives Bn=3, Bt=4:
    # traction.n=(9-16)/(2mu), tangent traction=(Bn*Bt)/mu in y.
    summary = maxwell_traction_summary((3.0, 4.0, 0.0), (1.0, 0.0, 0.0))

    assert summary["B_normal_T"] == pytest.approx(3.0)
    assert summary["B_tangent_T"] == pytest.approx(4.0)
    assert summary["normal_traction_Pa"] == pytest.approx(-3.5 / MU0)
    assert summary["normal_traction_identity_Pa"] == pytest.approx(-3.5 / MU0)
    assert summary["tangential_traction_Pa"] == pytest.approx([0.0, 12.0 / MU0, 0.0])
    assert summary["tangential_traction_magnitude_Pa"] == pytest.approx(12.0 / MU0)


def test_air_gap_shear_stress_matches_maxwell_tangential_traction():
    Br = 0.8
    Bt = 0.1
    shear = air_gap_shear_stress(Br, Bt)
    traction = maxwell_traction_summary((Br, Bt, 0.0), (1.0, 0.0, 0.0))

    assert shear == pytest.approx(Br * Bt / MU0)
    assert traction["tangential_traction_Pa"] == pytest.approx([0.0, shear, 0.0])
    assert traction["tangential_traction_magnitude_Pa"] == pytest.approx(abs(shear))


def test_air_gap_shear_torque_scales_with_radius_length_angle_and_sign():
    Br = 0.8
    Bt = 0.1
    radius = 0.05
    length = 0.1
    full = air_gap_shear_torque(Br, Bt, radius, axial_length_m=length)
    half = air_gap_shear_torque(Br, Bt, radius, axial_length_m=length, angle_rad=math.pi)
    reverse = air_gap_shear_torque(Br, -Bt, radius, axial_length_m=length)
    summary = air_gap_shear_torque_summary(Br, Bt, radius, axial_length_m=length)

    assert full == pytest.approx(100.0)
    assert half == pytest.approx(50.0)
    assert reverse == pytest.approx(-100.0)
    assert summary["surface_area_m2"] == pytest.approx(2.0 * math.pi * radius * length)
    assert summary["tangential_force_N"] == pytest.approx(full / radius)
    assert summary["torque_Nm"] == pytest.approx(full)
    assert summary["torque_per_axial_length_N"] == pytest.approx(full / length)


def test_maxwell_line_segment_force_2d_matches_air_gap_pressure():
    pressure = air_gap_maxwell_pressure(1.0)
    row = maxwell_line_segment_force_2d(
        (0.0, -0.5),
        (0.0, 0.5),
        (1.0, 0.0),
        normal_side="right",
    )

    assert row["length_m"] == pytest.approx(1.0)
    assert row["unit_normal"] == pytest.approx([1.0, 0.0])
    assert row["traction_N_per_m2"] == pytest.approx([pressure, 0.0])
    assert row["force_per_depth_N_per_m"] == pytest.approx([pressure, 0.0])
    assert row["normal_force_per_depth_N_per_m"] == pytest.approx(pressure)


def test_maxwell_contour_force_2d_closed_uniform_field_cancels():
    contour = [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]
    summary = maxwell_contour_force_2d(contour, (1.0, 0.0), orientation="ccw")

    assert summary["n_segments"] == 4
    assert summary["polygon_signed_area_m2"] == pytest.approx(2.0)
    assert summary["total_force_per_depth_N_per_m"] == pytest.approx([0.0, 0.0], abs=1.0e-9)
    assert summary["total_force_magnitude_per_depth_N_per_m"] == pytest.approx(0.0, abs=1.0e-9)
    assert summary["sum_abs_normal_force_per_depth_N_per_m"] > 0.0


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
    with pytest.raises(ValueError):
        maxwell_stress_tensor_air((1.0,), mu=MU0)
    with pytest.raises(ValueError):
        maxwell_traction_air((1.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        maxwell_traction_summary((1.0, 0.0), (1.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        air_gap_shear_stress(1.0, 0.1, mu=0.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque(1.0, 0.1, radius_m=-1.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque_summary(1.0, 0.1, radius_m=1.0, axial_length_m=-1.0)
    with pytest.raises(ValueError):
        maxwell_line_segment_force_2d((0.0, 0.0), (0.0, 0.0), (1.0, 0.0))
    with pytest.raises(ValueError):
        maxwell_contour_force_2d([(0.0, 0.0), (1.0, 0.0)], (1.0, 0.0))


if __name__ == "__main__":
    test_air_gap_pressure_matches_maxwell_stress_at_one_tesla()
    test_maxwell_tensor_normal_field_reduces_to_air_gap_pressure()
    test_maxwell_tensor_tangential_field_is_magnetic_tension()
    test_maxwell_traction_oblique_field_decomposes_into_normal_and_tangent()
    test_air_gap_shear_stress_matches_maxwell_tangential_traction()
    test_air_gap_shear_torque_scales_with_radius_length_angle_and_sign()
    test_maxwell_line_segment_force_2d_matches_air_gap_pressure()
    test_maxwell_contour_force_2d_closed_uniform_field_cancels()
    test_air_gap_force_scales_with_b_squared_area_and_faces()
    test_air_gap_force_summary_is_json_friendly_and_self_consistent()
    test_air_gap_force_rejects_invalid_inputs()
    print("[OK] air-gap Maxwell pressure and holding-force helpers validated.")
