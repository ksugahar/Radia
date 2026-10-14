import math

import numpy as np
import pytest

from radia.force import (
    MU0,
    air_gap_shear_torque,
    air_gap_shear_torque_from_angle_samples,
    coenergy_torque_from_angle_samples,
    force_torque_result,
    integrate_lorentz_force,
    integrate_lorentz_force_and_torque,
    integrate_maxwell_surface_force,
    integrate_maxwell_surface_force_and_torque,
    integrate_time_average_lorentz_force_and_torque,
    integrate_time_average_maxwell_surface_force,
    lorentz_force_density,
    maxwell_stress_tensor_air,
    maxwell_traction_air,
    time_average_air_gap_shear_torque_from_angle_samples,
    time_average_lorentz_force_density,
    time_average_maxwell_stress_tensor_air,
    virtual_work_force_from_displacement_samples,
)


def _line_current_field(points, center, current_A):
    relative = np.asarray(points, dtype=float) - np.asarray(center, dtype=float)
    radius_squared = relative[:, 0] ** 2 + relative[:, 1] ** 2
    field = np.zeros_like(relative)
    scale = MU0 * current_A / (2.0 * math.pi * radius_squared)
    field[:, 0] = -scale * relative[:, 1]
    field[:, 1] = scale * relative[:, 0]
    return field


def test_lorentz_density_and_volume_integral_have_si_sign_and_units():
    current_density = np.array([0.0, 0.0, 2.0e6])
    magnetic_flux_density = np.array([0.0, 0.3, 0.0])
    density = lorentz_force_density(current_density, magnetic_flux_density)
    np.testing.assert_allclose(density, [-6.0e5, 0.0, 0.0])
    np.testing.assert_allclose(
        integrate_lorentz_force(current_density, magnetic_flux_density, 2.5e-6),
        [-1.5, 0.0, 0.0],
    )


def test_normal_field_has_positive_maxwell_pressure():
    field = np.array([0.0, 0.0, 1.2])
    normal = np.array([0.0, 0.0, 4.0])
    expected_pressure = 1.2**2 / (2.0 * MU0)
    tensor = maxwell_stress_tensor_air(field)
    traction = maxwell_traction_air(field, normal)
    np.testing.assert_allclose(
        tensor,
        np.diag([-expected_pressure, -expected_pressure, expected_pressure]),
    )
    np.testing.assert_allclose(traction, [0.0, 0.0, expected_pressure])


def test_two_wire_lorentz_and_maxwell_routes_match_analytic_force_per_length():
    source_current_A = 80.0
    target_current_A = 35.0
    separation_m = 0.06
    target_radius_m = 0.008
    expected = -MU0 * source_current_A * target_current_A / (
        2.0 * math.pi * separation_m
    )

    target_area_m2 = math.pi * target_radius_m**2
    target_current_density = np.array(
        [0.0, 0.0, target_current_A / target_area_m2]
    )
    source_field_at_target = np.array(
        [0.0, MU0 * source_current_A / (2.0 * math.pi * separation_m), 0.0]
    )
    lorentz = integrate_lorentz_force(
        target_current_density,
        source_field_at_target,
        target_area_m2,
    )

    sample_count = 4096
    angle = 2.0 * math.pi * (np.arange(sample_count) + 0.5) / sample_count
    normals = np.column_stack((np.cos(angle), np.sin(angle), np.zeros(sample_count)))
    target_center = np.array([separation_m, 0.0, 0.0])
    contour_points = target_center + target_radius_m * normals
    total_field = _line_current_field(
        contour_points,
        np.zeros(3),
        source_current_A,
    ) + _line_current_field(
        contour_points,
        target_center,
        target_current_A,
    )
    line_weights_for_unit_depth_m2 = np.full(
        sample_count,
        target_radius_m * 2.0 * math.pi / sample_count,
    )
    maxwell = integrate_maxwell_surface_force(
        total_field,
        normals,
        line_weights_for_unit_depth_m2,
    )

    np.testing.assert_allclose(lorentz, [expected, 0.0, 0.0], rtol=1.0e-14, atol=1.0e-18)
    np.testing.assert_allclose(maxwell, [expected, 0.0, 0.0], rtol=2.0e-13, atol=1.0e-18)


def test_lorentz_force_and_torque_use_physical_points_and_pivot():
    current_density = [[1.0, 0.0, 0.0]]
    magnetic_flux_density = [[0.0, 1.0, 0.0]]
    point = [[0.0, 1.0, 0.0]]
    force, torque = integrate_lorentz_force_and_torque(
        current_density,
        magnetic_flux_density,
        [2.0],
        point,
    )
    shifted_force, shifted_torque = integrate_lorentz_force_and_torque(
        current_density,
        magnetic_flux_density,
        [2.0],
        point,
        pivot_m=[0.0, 0.5, 0.0],
    )

    np.testing.assert_allclose(force, [0.0, 0.0, 2.0])
    np.testing.assert_allclose(torque, [2.0, 0.0, 0.0])
    np.testing.assert_allclose(shifted_force, force)
    np.testing.assert_allclose(shifted_torque, [1.0, 0.0, 0.0])


def test_maxwell_surface_force_and_torque_obey_resultant_moment_identity():
    field = [[0.8, 0.2, 0.0], [0.4, -0.1, 0.0]]
    normals = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    weights = [0.25, 0.5]
    points = [[0.0, 0.1, 0.0], [0.2, 0.0, 0.0]]
    force, torque = integrate_maxwell_surface_force_and_torque(
        field,
        normals,
        weights,
        points,
    )
    expected_force = integrate_maxwell_surface_force(field, normals, weights)
    tractions = maxwell_traction_air(field, normals)
    expected_torque = np.sum(
        np.cross(np.asarray(points), tractions) * np.asarray(weights)[:, None],
        axis=0,
    )

    np.testing.assert_allclose(force, expected_force)
    np.testing.assert_allclose(torque, expected_torque)


def test_peak_phasor_force_and_stress_reduce_to_half_static_real_fields():
    current = np.array([[0.0, 0.0, 3.0]], dtype=complex)
    field = np.array([[0.0, 0.4, 0.0]], dtype=complex)
    density_static = lorentz_force_density(current.real, field.real)
    density_peak = time_average_lorentz_force_density(current, field, amplitude="peak")
    density_rms = time_average_lorentz_force_density(current, field, amplitude="rms")
    np.testing.assert_allclose(density_peak, 0.5 * density_static)
    np.testing.assert_allclose(density_rms, density_static)

    b = np.array([0.3, -0.2, 0.5], dtype=complex)
    np.testing.assert_allclose(
        time_average_maxwell_stress_tensor_air(b, amplitude="peak"),
        0.5 * maxwell_stress_tensor_air(b.real),
    )
    np.testing.assert_allclose(
        integrate_time_average_maxwell_surface_force(
            [b],
            [[1.0, 0.0, 0.0]],
            [2.0],
            amplitude="peak",
        ),
        0.5
        * integrate_maxwell_surface_force(
            [b.real],
            [[1.0, 0.0, 0.0]],
            [2.0],
        ),
    )


def test_phasor_lorentz_phase_and_torque_are_conjugation_consistent():
    current = [[0.0, 0.0, 2.0 + 0.0j]]
    quadrature_field = [[0.0, 0.0 + 0.5j, 0.0]]
    force, torque = integrate_time_average_lorentz_force_and_torque(
        current,
        quadrature_field,
        [3.0],
        [[0.0, 0.0, 1.0]],
        amplitude="peak",
    )
    np.testing.assert_allclose(force, np.zeros(3), atol=1.0e-15)
    np.testing.assert_allclose(torque, np.zeros(3), atol=1.0e-15)


def test_virtual_work_coenergy_torque_and_air_gap_torque_cover_motor_maglev_routes():
    positions = np.linspace(-0.002, 0.002, 5)
    force_expected = 7.5
    coenergy = 0.25 + force_expected * positions
    np.testing.assert_allclose(
        virtual_work_force_from_displacement_samples(positions, coenergy),
        np.full(positions.shape, force_expected),
    )

    angles = np.linspace(0.0, 0.4, 5)
    torque_expected = -2.25
    angular_coenergy = 1.0 + torque_expected * angles
    np.testing.assert_allclose(
        coenergy_torque_from_angle_samples(angles, angular_coenergy),
        np.full(angles.shape, torque_expected),
    )
    assert air_gap_shear_torque(0.8, 0.1, 0.05, axial_length_m=0.1) == pytest.approx(
        100.0
    )


def test_sampled_static_and_phasor_air_gap_torque_cover_a_full_motor_period():
    angles = np.arange(4) * 0.5 * math.pi
    radial = np.full(4, 0.8)
    tangential = np.full(4, 0.1)
    static = air_gap_shear_torque_from_angle_samples(
        angles,
        radial,
        tangential,
        0.05,
        axial_length_m=0.1,
    )
    phasor = time_average_air_gap_shear_torque_from_angle_samples(
        angles,
        radial.astype(complex),
        tangential.astype(complex),
        0.05,
        axial_length_m=0.1,
        amplitude="peak",
    )
    assert static["integrated_angle_rad"] == pytest.approx(2.0 * math.pi)
    assert static["torque_Nm"] == pytest.approx(100.0)
    assert phasor["torque_Nm"] == pytest.approx(50.0)


def test_force_torque_result_records_frame_pivot_units_and_phasor_convention():
    result = force_torque_result(
        [1.0, 2.0, 3.0],
        [0.1, 0.2, 0.3],
        method="time_average_lorentz",
        frame="rotor",
        pivot_m=[0.0, 0.0, 0.1],
        field_convention="time_average_phasor",
        amplitude="peak",
    )
    assert result["schema"] == "radia.force-result/v1"
    assert result["frame"] == "rotor"
    assert result["pivot_m"] == [0.0, 0.0, 0.1]
    assert result["force_N"] == [1.0, 2.0, 3.0]
    assert result["torque_Nm"] == [0.1, 0.2, 0.3]
    assert result["phasor_amplitude"] == "peak"


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: lorentz_force_density([1.0, 0.0], [0.0, 0.0, 1.0]), "shape"),
        (lambda: maxwell_traction_air([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]), "nonzero"),
        (
            lambda: integrate_lorentz_force(
                [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]],
                [0.0, 1.0, 0.0],
                [1.0, -1.0],
            ),
            ">= 0",
        ),
    ],
)
def test_invalid_force_contracts_fail_loudly(call, message):
    with pytest.raises(ValueError, match=message):
        call()
