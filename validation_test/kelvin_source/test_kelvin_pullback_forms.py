"""Validation locks for Kelvin A/B pullbacks and source evaluation.

The former A-formulation Kelvin examples contained full research runners for
Radia-vs-FEM and energy checks.  The maintained public surface is now the
``radia.kelvin_source`` API.  These tests keep the mathematical contract close
to that API without depending on a heavy FEM solve.
"""

from __future__ import annotations

import math

import numpy as np

from radia.kelvin_source import (
    A_s_at_obs_with_kelvin,
    biot_savart_A_at_points,
    kelvin_factor_scalar,
    kelvin_map_3d,
    kelvin_nu_factor_3d,
    kelvin_pullback_B_pseudovector,
    kelvin_pullback_vector,
)


def _sample_kelvin_points(center, radius, n=16):
    rng = np.random.default_rng(20260628)
    directions = rng.normal(size=(n, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    rhos = np.linspace(0.25 * radius, 0.90 * radius, n)
    return np.asarray(center, dtype=float) + rhos[:, None] * directions


def _build_loop(radius=0.030, z0=0.0, segments=96, current=1.0):
    theta = np.linspace(0.0, 2.0 * math.pi, segments + 1)
    pts = np.column_stack([
        radius * np.cos(theta),
        radius * np.sin(theta),
        np.full_like(theta, z0),
    ])
    return [[(pts[i].tolist(), pts[i + 1].tolist()) for i in range(segments)]], [current]


def test_1form_pullback_splits_tangential_and_radial_parts():
    center = np.array([0.15, -0.02, 0.03])
    radius = 0.06
    pts = _sample_kelvin_points(center, radius)
    d = pts - center
    n = d / np.linalg.norm(d, axis=1, keepdims=True)

    rng = np.random.default_rng(7)
    v = rng.normal(size=pts.shape)
    tangent = v - np.sum(v * n, axis=1, keepdims=True) * n
    radial = rng.normal(size=(pts.shape[0], 1)) * n
    scale = kelvin_factor_scalar(pts, center, radius)[:, None]

    np.testing.assert_allclose(
        kelvin_pullback_vector(tangent, pts, center, radius),
        scale * tangent,
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        kelvin_pullback_vector(radial, pts, center, radius),
        -scale * radial,
        rtol=1e-12,
        atol=1e-12,
    )


def test_1form_pullback_is_involutive_with_kelvin_map():
    center = np.array([0.0, 0.1, -0.2])
    radius = 1.0
    rng = np.random.default_rng(11)
    r_phys = center + rng.normal(size=(12, 3))
    r_phys += 2.5 * (r_phys - center) / np.linalg.norm(
        r_phys - center, axis=1, keepdims=True)
    a_phys = rng.normal(size=(12, 3))

    r_prime = kelvin_map_3d(r_phys, center, radius)
    a_comp = kelvin_pullback_vector(a_phys, r_prime, center, radius)
    a_back = kelvin_pullback_vector(a_comp, r_phys, center, radius)

    np.testing.assert_allclose(kelvin_map_3d(r_prime, center, radius), r_phys)
    np.testing.assert_allclose(a_back, a_phys, rtol=1e-12, atol=1e-12)


def test_2form_pullback_axis_signs_and_magnitudes():
    center = np.zeros(3)
    radius = 1.0
    b_phys = np.array([0.0, 0.0, 2.5])

    for rho in (0.2, 0.4, 0.8):
        rp_x = np.array([rho, 0.0, 0.0])
        expected_x = -((radius / rho) ** 4) * b_phys
        np.testing.assert_allclose(
            kelvin_pullback_B_pseudovector(b_phys, rp_x, center, radius),
            expected_x,
            rtol=1e-12,
            atol=1e-12,
        )

        rp_z = np.array([0.0, 0.0, rho])
        expected_z = ((radius / rho) ** 4) * b_phys
        np.testing.assert_allclose(
            kelvin_pullback_B_pseudovector(b_phys, rp_z, center, radius),
            expected_z,
            rtol=1e-12,
            atol=1e-12,
        )


def test_kelvin_energy_density_identity_is_local():
    """Nagamine/Sugahara factor makes the B-pullback energy invariant locally."""
    center = np.array([0.05, -0.04, 0.02])
    radius = 0.2
    pts = _sample_kelvin_points(center, radius, n=10)
    rng = np.random.default_rng(19)
    b_phys = rng.normal(size=pts.shape)
    nu0 = 1.0 / (4.0 * math.pi * 1e-7)

    rho_prime = np.linalg.norm(pts - center, axis=1)
    jac_abs = (radius / rho_prime) ** 6
    b_comp = kelvin_pullback_B_pseudovector(b_phys, pts, center, radius)
    nu_comp = nu0 * kelvin_nu_factor_3d(pts, center, radius)

    physical_density_times_jac = 0.5 * nu0 * np.sum(b_phys * b_phys, axis=1) * jac_abs
    computational_density = 0.5 * nu_comp * np.sum(b_comp * b_comp, axis=1)
    np.testing.assert_allclose(
        computational_density,
        physical_density_times_jac,
        rtol=1e-12,
        atol=1e-5,
    )


def test_kelvin_aware_biot_savart_uses_direct_inner_and_pullback_outer():
    center = np.array([0.15, 0.0, 0.0])
    radius = 0.06
    paths, currents = _build_loop()

    inner = np.array([0.010, 0.0, 0.005])
    outer_comp = center + np.array([0.020, 0.005, 0.010])
    obs = np.vstack([inner, outer_comp])

    got = A_s_at_obs_with_kelvin(
        obs,
        paths,
        currents,
        kelvin_center=center,
        R_kelvin=radius,
        factor_mode="pullback",
        n_quad=6,
    )

    expected_inner = biot_savart_A_at_points(inner.reshape(1, 3), paths, currents,
                                             n_quad=6)[0]
    outer_phys = kelvin_map_3d(outer_comp, center, radius)
    outer_a = biot_savart_A_at_points(outer_phys.reshape(1, 3), paths, currents,
                                      n_quad=6)[0]
    expected_outer = kelvin_pullback_vector(outer_a, outer_comp, center, radius)

    np.testing.assert_allclose(got[0], expected_inner, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(got[1], expected_outer, rtol=1e-12, atol=1e-14)
