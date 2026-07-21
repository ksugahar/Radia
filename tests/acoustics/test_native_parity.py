"""Independent SciPy parity checks for the shared acoustic C++ kernel."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import eval_legendre, hankel1, jv, spherical_jn, spherical_yn

import radia.acoustics as acoustic
from radia import _radia_pybind as native
from radia_mcp.radia_acoustic import cq_grid_gate
from radia.acoustics.cq import (
    bdf_delta,
    soft_sphere_scattering_complex_k,
)


def _jn(order, value):
    return spherical_jn(order, value)


def _jn_d(order, value):
    return spherical_jn(order, value, derivative=True)


def _jn_dd(order, value):
    return -(2.0 / value) * _jn_d(order, value) - (
        1.0 - order * (order + 1) / value**2
    ) * _jn(order, value)


def _h1(order, value):
    return spherical_jn(order, value) + 1j * spherical_yn(order, value)


def _h1_d(order, value):
    return spherical_jn(order, value, derivative=True) + 1j * spherical_yn(
        order, value, derivative=True
    )


def _point_data(points):
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    radius = np.linalg.norm(points, axis=1)
    return points, radius, points[:, 2] / np.maximum(radius, 1e-30)


def _soft_reference(wavenumber, radius, points, terms):
    _, distance, cosine = _point_data(points)
    scattered = np.zeros(distance.shape, dtype=complex)
    for order in range(terms + 1):
        coefficient = (
            -(1j**order)
            * (2 * order + 1)
            * _jn(order, wavenumber * radius)
            / _h1(order, wavenumber * radius)
        )
        scattered += coefficient * _h1(order, wavenumber * distance) * eval_legendre(
            order, cosine
        )
    return scattered


def _rigid_reference(wavenumber, radius, points, terms):
    _, distance, cosine = _point_data(points)
    scattered = np.zeros(distance.shape, dtype=complex)
    for order in range(terms + 1):
        coefficient = (
            -(1j**order)
            * (2 * order + 1)
            * _jn_d(order, wavenumber * radius)
            / _h1_d(order, wavenumber * radius)
        )
        scattered += coefficient * _h1(order, wavenumber * distance) * eval_legendre(
            order, cosine
        )
    return scattered


def _fluid_reference(
    wavenumber, radius, points, interior_wavenumber, density_ratio, terms
):
    _, distance, cosine = _point_data(points)
    inside = distance <= radius * (1.0 + 1e-12)
    total = np.zeros(distance.shape, dtype=complex)
    x0 = wavenumber * radius
    x1 = interior_wavenumber * radius
    for order in range(terms + 1):
        incident_coefficient = (1j**order) * (2 * order + 1)
        beta = (
            (interior_wavenumber / density_ratio)
            * _jn_d(order, x1)
            / _jn(order, x1)
        )
        scattered_coefficient = -incident_coefficient * (
            wavenumber * _jn_d(order, x0) - beta * _jn(order, x0)
        ) / (wavenumber * _h1_d(order, x0) - beta * _h1(order, x0))
        interior_coefficient = (
            incident_coefficient * _jn(order, x0)
            + scattered_coefficient * _h1(order, x0)
        ) / _jn(order, x1)
        legendre = eval_legendre(order, cosine)
        total[inside] += (
            interior_coefficient
            * _jn(order, interior_wavenumber * np.maximum(distance[inside], 1e-30))
            * legendre[inside]
        )
        total[~inside] += (
            incident_coefficient * _jn(order, wavenumber * distance[~inside])
            + scattered_coefficient * _h1(order, wavenumber * distance[~inside])
        ) * legendre[~inside]
    return total


def _elastic_coefficient(
    order, wavenumber, radius, longitudinal_speed, shear_speed, density_ratio
):
    omega = wavenumber
    k_longitudinal = omega / longitudinal_speed
    x = wavenumber * radius
    xl = k_longitudinal * radius
    fluid_factor = wavenumber / omega**2
    mu = density_ratio * shear_speed**2
    lame_lambda = density_ratio * (longitudinal_speed**2 - 2 * shear_speed**2)
    if shear_speed == 0:
        matrix = np.array(
            [
                [fluid_factor * _h1_d(order, x), -k_longitudinal * _jn_d(order, xl)],
                [
                    _h1(order, x),
                    -lame_lambda * k_longitudinal**2 * _jn(order, xl),
                ],
            ],
            dtype=complex,
        )
        rhs = np.array(
            [-fluid_factor * _jn_d(order, x), -_jn(order, x)], dtype=complex
        )
        return np.linalg.solve(matrix, rhs)[0]

    k_transverse = omega / shear_speed
    xt = k_transverse * radius
    angular = order * (order + 1)
    ur_a = k_longitudinal * _jn_d(order, xl)
    ur_b = angular / radius * _jn(order, xt)
    dur_a = k_longitudinal**2 * _jn_dd(order, xl)
    dur_b = angular * (
        k_transverse * _jn_d(order, xt) / radius
        - _jn(order, xt) / radius**2
    )
    srr_a = -lame_lambda * k_longitudinal**2 * _jn(order, xl) + 2 * mu * dur_a
    srr_b = 2 * mu * dur_b
    va_a = _jn(order, xl) / radius
    va_b = _jn(order, xt) / radius + k_transverse * _jn_d(order, xt)
    vp_a = -_jn(order, xl) / radius**2 + k_longitudinal * _jn_d(order, xl) / radius
    vp_b = (
        -_jn(order, xt) / radius**2
        + k_transverse * _jn_d(order, xt) / radius
        + k_transverse**2 * _jn_dd(order, xt)
    )
    srt_a = mu * (ur_a / radius + vp_a - va_a / radius)
    srt_b = mu * (ur_b / radius + vp_b - va_b / radius)
    matrix = np.array(
        [
            [fluid_factor * _h1_d(order, x), -ur_a, -ur_b],
            [_h1(order, x), srr_a, srr_b],
            [0.0, srt_a, srt_b],
        ],
        dtype=complex,
    )
    rhs = np.array(
        [-fluid_factor * _jn_d(order, x), -_jn(order, x), 0.0], dtype=complex
    )
    return np.linalg.solve(matrix, rhs)[0]


def _elastic_reference(
    wavenumber, radius, points, longitudinal_speed, shear_speed, density_ratio, terms
):
    _, distance, cosine = _point_data(points)
    scattered = np.zeros(distance.shape, dtype=complex)
    for order in range(terms + 1):
        coefficient = _elastic_coefficient(
            order,
            wavenumber,
            radius,
            longitudinal_speed,
            shear_speed,
            density_ratio,
        )
        scattered += (
            (1j**order)
            * (2 * order + 1)
            * coefficient
            * _h1(order, wavenumber * distance)
            * eval_legendre(order, cosine)
        )
    return scattered


def _complex_spherical_j(order, value):
    return np.sqrt(np.pi / (2 * value)) * jv(order + 0.5, value)


def _complex_spherical_h(order, value):
    return np.sqrt(np.pi / (2 * value)) * hankel1(order + 0.5, value)


def _complex_soft_reference(wavenumber, radius, points, terms):
    _, distance, cosine = _point_data(points)
    scattered = np.zeros(len(points), dtype=complex)
    for order in range(terms + 1):
        coefficient = (
            -(1j**order)
            * (2 * order + 1)
            * _complex_spherical_j(order, wavenumber * radius)
            / _complex_spherical_h(order, wavenumber * radius)
        )
        scattered += (
            coefficient
            * _complex_spherical_h(order, wavenumber * distance)
            * eval_legendre(order, cosine)
        )
    return scattered


def test_native_real_sphere_kernels_match_independent_scipy_series():
    wavenumber = 2.3
    radius = 0.8
    exterior = np.array([[0.0, 0.0, 1.1], [1.2, 0.0, 0.4], [-0.3, 1.4, 0.2]])

    soft = acoustic.soft_sphere_scattering(wavenumber, radius, exterior, terms=17)
    np.testing.assert_allclose(
        soft["scattered"], _soft_reference(wavenumber, radius, exterior, 17),
        rtol=2e-12, atol=2e-13,
    )
    assert soft["backend"] == "native-cpp-pybind11"

    rigid = acoustic.rigid_sphere_scattering(wavenumber, radius, exterior, terms=17)
    np.testing.assert_allclose(
        rigid["scattered"], _rigid_reference(wavenumber, radius, exterior, 17),
        rtol=2e-12, atol=2e-13,
    )

    mixed = np.vstack(([[0.0, 0.0, 0.3]], exterior))
    fluid = acoustic.fluid_sphere_scattering(
        wavenumber, radius, mixed, interior_wavenumber=1.35, density_ratio=1.7,
        terms=17,
    )
    np.testing.assert_allclose(
        fluid["total"],
        _fluid_reference(wavenumber, radius, mixed, 1.35, 1.7, 17),
        rtol=8e-12, atol=8e-13,
    )

    elastic = acoustic.elastic_sphere_scattering(
        wavenumber, radius, exterior, longitudinal_speed=2.4, shear_speed=1.1,
        density_ratio=1.8, terms=14,
    )
    np.testing.assert_allclose(
        elastic["scattered"],
        _elastic_reference(wavenumber, radius, exterior, 2.4, 1.1, 1.8, 14),
        rtol=2e-11, atol=2e-12,
    )


def test_native_complex_wavenumber_and_cq_grid_match_independent_formulas():
    wavenumber = 0.75 + 0.45j
    radius = 0.9
    points = np.array([[0.0, 0.0, 1.4], [1.3, 0.0, -0.2]])
    expected = _complex_soft_reference(wavenumber, radius, points, 18)
    actual = soft_sphere_scattering_complex_k(
        wavenumber, radius, points, terms=18
    )
    np.testing.assert_allclose(actual, expected, rtol=3e-11, atol=3e-12)

    zeta = np.array([0.1 + 0.2j, -0.4 + 0.3j])
    np.testing.assert_allclose(
        bdf_delta(zeta, "BDF2"), 1.5 - 2 * zeta + 0.5 * zeta**2,
        rtol=0, atol=2e-15,
    )
    assert np.isscalar(bdf_delta(0.25 + 0.1j, "BDF1"))

    cq_points = np.array([[0.0, 0.0, -3.0], [3.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    grid = native._AcousticCQGrid(16, 0.28, 1.0, "BDF2")
    for cq_wavenumber in grid["cq_wavenumbers"]:
        actual = soft_sphere_scattering_complex_k(
            complex(cq_wavenumber), 1.0, cq_points, terms=28
        )
        expected = _complex_soft_reference(
            complex(cq_wavenumber), 1.0, cq_points, 28
        )
        np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)

    with pytest.raises(ValueError, match="terms must be -1"):
        acoustic.soft_sphere_scattering(2.0, 1.0, cq_points, terms=-2)


@pytest.mark.parametrize(("method", "num_time"), (("BDF1", 15), ("BDF2", 16)))
def test_minimal_mcp_cq_gate_matches_native_python_grid(method, num_time):
    time_step = 0.08
    sound_speed = 1.2
    gate = cq_grid_gate(
        num_time=num_time,
        time_step=time_step,
        sound_speed=sound_speed,
        method=method,
    )
    grid = native._AcousticCQGrid(num_time, time_step, sound_speed, method)
    nodes = np.asarray(grid["cq_nodes"])
    wavenumbers = np.asarray(grid["cq_wavenumbers"])

    assert gate["ok"]
    assert gate["method"] == method
    assert gate["cq_radius"] == pytest.approx(grid["cq_radius"], abs=2e-15)
    assert gate["min_real_s"] == pytest.approx(np.min(nodes.real), abs=2e-14)
    assert gate["max_abs_kappa"] == pytest.approx(
        np.max(np.abs(wavenumbers)), abs=2e-14
    )
