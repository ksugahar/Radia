"""Common observation-tube contracts for ESRF/FFAG field comparisons."""

from __future__ import annotations

import numpy as np
import pytest

from radia.accelerator_field_validation import (
    CurvilinearObservationTube,
    MagnetFieldEngine,
    MagnetFieldSample,
    circular_transverse_offsets,
    compare_integrated_multipole_rows,
    longitudinal_reversal_symmetry,
    compare_magnetic_flux_density,
    project_straight_quadrupole_symmetry,
    radial_field_index,
    sample_field_engine,
    transverse_multipole_spectrum,
)


def test_straight_quadrupole_projection_removes_symmetry_forbidden_noise():
    points = np.array([
        [0.03, 0.002, 0.001],
        [-0.02, -0.001, 0.003],
        [0.00, 0.004, -0.002],
    ])
    gradient = 12.0

    def noisy_field(query):
        values = np.column_stack((
            np.zeros(query.shape[0]),
            gradient * query[:, 2],
            gradient * query[:, 1],
        ))
        return values + np.array([0.4, -0.2, 0.3])

    projected = project_straight_quadrupole_symmetry(
        points, noisy_field, longitudinal_axis=0)
    expected = np.column_stack((
        np.zeros(points.shape[0]),
        gradient * points[:, 2],
        gradient * points[:, 1],
    ))
    assert projected == pytest.approx(expected, abs=2.0e-14)


def test_straight_quadrupole_projection_validates_contract():
    with pytest.raises(ValueError, match="longitudinal_axis"):
        project_straight_quadrupole_symmetry(
            np.zeros((1, 3)), lambda points: points, longitudinal_axis=3)
    with pytest.raises(ValueError, match="return"):
        project_straight_quadrupole_symmetry(
            np.zeros((2, 3)), lambda points: np.zeros((1, 3)))


def _straight_tube(offsets: np.ndarray) -> CurvilinearObservationTube:
    s = np.linspace(0.0, 0.2, 5)
    return CurvilinearObservationTube(
        station_s=s,
        center=np.column_stack((np.zeros_like(s), np.zeros_like(s), s)),
        tangent=np.tile([0.0, 0.0, 1.0], (s.size, 1)),
        normal=np.tile([1.0, 0.0, 0.0], (s.size, 1)),
        binormal=np.tile([0.0, 1.0, 0.0], (s.size, 1)),
        transverse_offsets=offsets,
    )


def test_common_tube_pairwise_b_residual_is_gauge_independent():
    tube = _straight_tube(circular_transverse_offsets(0.02, 16))

    def exact(points):
        return np.column_stack((points[:, 1], points[:, 0], 1.5 + 0 * points[:, 2]))

    def shifted(points):
        value = exact(points)
        value[:, 1] += 2.0e-4
        return value

    left = sample_field_engine(MagnetFieldEngine("native_radia", exact), tube)
    right = sample_field_engine(
        MagnetFieldEngine("hdiv_mmm", shifted, vector_potential=lambda p: 0 * p), tube)
    result = compare_magnetic_flux_density((left, right))
    assert result["raw_vector_potential_compared"] is False
    assert result["pairs"][0]["maximum_vector_error_t"] == pytest.approx(2.0e-4)


def test_manual_field_sample_validates_shape_finiteness_and_a_pair():
    tube = _straight_tube(circular_transverse_offsets(0.02, 8))
    expected = (tube.station_count, tube.transverse_point_count, 3)
    values = np.zeros(expected)
    sample = MagnetFieldSample("checked", tube, values, values)
    assert not sample.b_global.flags.writeable
    with pytest.raises(ValueError, match="b_local must have shape"):
        MagnetFieldSample("bad-shape", tube, values, values[:, :-1])
    nonfinite = values.copy()
    nonfinite[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="b_global must be finite"):
        MagnetFieldSample("bad-value", tube, nonfinite, values)
    with pytest.raises(ValueError, match="supplied together"):
        MagnetFieldSample("half-a", tube, values, values, a_global=values)


def test_quadrupole_multipole_is_recovered_in_local_frame():
    gradient = 7.5
    tube = _straight_tube(circular_transverse_offsets(0.01, 24))

    def quadrupole(points):
        return np.column_stack((gradient * points[:, 1],
                                gradient * points[:, 0],
                                np.zeros(points.shape[0])))

    sample = sample_field_engine(MagnetFieldEngine("quadrupole", quadrupole), tube)
    coefficients = transverse_multipole_spectrum(sample, maximum_order=4)
    np.testing.assert_allclose(coefficients[:, 1], gradient, atol=2.0e-13)
    np.testing.assert_allclose(coefficients[:, [0, 2, 3]], 0.0, atol=3.0e-12)


def test_multipole_fit_accepts_nonuniform_circle_samples():
    radius = 0.01
    angles = np.array([
        0.00, 0.17, 0.61, 1.04, 1.77, 2.11,
        2.83, 3.42, 4.08, 4.71, 5.29, 5.91,
    ])
    offsets = radius * np.column_stack((np.cos(angles), np.sin(angles)))
    tube = _straight_tube(offsets)
    expected = np.array([0.2 + 0.1j, 7.5 - 0.3j, -12.0 + 2.0j])

    def field(points):
        coordinate = points[:, 0] + 1j * points[:, 1]
        value = sum(
            coefficient * coordinate ** power
            for power, coefficient in enumerate(expected)
        )
        return np.column_stack((value.imag, value.real, np.zeros(value.size)))

    sample = sample_field_engine(MagnetFieldEngine("nonuniform", field), tube)
    coefficients = transverse_multipole_spectrum(sample, maximum_order=3)
    np.testing.assert_allclose(
        coefficients,
        np.tile(expected, (tube.station_count, 1)),
        atol=2.0e-12,
    )


def test_multipole_fit_rejects_rank_deficient_circle_samples():
    offsets = np.tile([[0.01, 0.0]], (8, 1))
    tube = _straight_tube(offsets)
    values = np.zeros((tube.station_count, tube.transverse_point_count, 3))
    sample = MagnetFieldSample("rank-deficient", tube, values, values)
    with pytest.raises(ValueError, match="full-rank"):
        transverse_multipole_spectrum(sample, maximum_order=4)


def test_integrated_multipole_comparison_separates_main_error_and_units():
    reference = [
        {"order": 1, "integrated_real_t_m_per_m_power": 2.0,
         "integrated_imag_t_m_per_m_power": 0.0,
         "normal_units_at_reference_radius": 10000.0,
         "skew_units_at_reference_radius": 0.0},
        {"order": 2, "integrated_real_t_m_per_m_power": 0.1,
         "integrated_imag_t_m_per_m_power": 0.0,
         "normal_units_at_reference_radius": 5.0,
         "skew_units_at_reference_radius": 0.0},
    ]
    candidate = [
        {"order": 1, "integrated_real_t_m_per_m_power": 1.98,
         "integrated_imag_t_m_per_m_power": 0.0,
         "normal_units_at_reference_radius": 10000.0,
         "skew_units_at_reference_radius": 0.0},
        {"order": 2, "integrated_real_t_m_per_m_power": 0.12,
         "integrated_imag_t_m_per_m_power": 0.01,
         "normal_units_at_reference_radius": 6.0,
         "skew_units_at_reference_radius": 0.5},
    ]
    result = compare_integrated_multipole_rows(reference, candidate, main_order=1)
    assert result["main_relative_error"] == pytest.approx(0.01)
    assert result["harmonics_at_reference_radius"][1]["normal_units_difference"] == 1.0
    assert result["harmonics_at_reference_radius"][1]["skew_units_difference"] == 0.5


def test_integrated_multipole_comparison_rejects_missing_main_order():
    row = {"order": 1, "integrated_real_t_m_per_m_power": 1.0,
           "integrated_imag_t_m_per_m_power": 0.0,
           "normal_units_at_reference_radius": 10000.0,
           "skew_units_at_reference_radius": 0.0}
    with pytest.raises(ValueError, match="main_order 2"):
        compare_integrated_multipole_rows([row], [row], main_order=2)


def test_radial_field_index_recovers_power_law():
    radius = np.geomspace(0.01, 0.1, 21)
    sampled_radius, index = radial_field_index(radius, 3.0 * radius ** 2.5)
    np.testing.assert_allclose(sampled_radius, radius[1:-1])
    np.testing.assert_allclose(index, 2.5, atol=2.0e-14)


def test_left_handed_observation_frame_is_rejected():
    offsets = circular_transverse_offsets(0.01, 8)
    s = np.array([0.0])
    with pytest.raises(ValueError, match="tangent x normal"):
        CurvilinearObservationTube(
            station_s=s, center=np.zeros((1, 3)), tangent=np.array([[0., 0., 1.]]),
            normal=np.array([[1., 0., 0.]]), binormal=np.array([[0., -1., 0.]]),
            transverse_offsets=offsets,
        )


def test_longitudinal_reversal_symmetry_accepts_even_complex_profile():
    stations = np.linspace(-0.05, 0.05, 11)
    profile = (2.0 + 3.0j) * (1.0 - 20.0 * stations ** 2)
    result = longitudinal_reversal_symmetry(stations, profile)
    assert result["relative_rms_defect"] < 1.0e-14
    assert result["odd_to_even_l2"] < 1.0e-14


def test_longitudinal_reversal_symmetry_detects_spurious_odd_fringe():
    stations = np.linspace(-0.05, 0.05, 11)
    profile = 1.0 - 20.0 * stations ** 2 + 2.0 * stations
    result = longitudinal_reversal_symmetry(stations, profile)
    assert result["relative_rms_defect"] > 0.1
    with pytest.raises(ValueError, match="symmetric about zero"):
        longitudinal_reversal_symmetry(stations + 0.01, profile)
