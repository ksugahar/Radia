"""Validation locks for centralized Kelvin material-factor helpers.

The former Kelvin examples carried many inline variants of the same factors.
This validation lane pins the maintained ``radia.kelvin_source`` API so future
code and MCP knowledge can point here instead of historical docs archives.
"""

from __future__ import annotations

import numpy as np

from radia.kelvin_source import (
    kelvin_factor_2d_inplane,
    kelvin_mu_factor_2d_axial,
    kelvin_mu_factor_3d,
    kelvin_mu_factor_axisym,
    kelvin_nu_factor_2d_axial,
    kelvin_nu_factor_3d,
    kelvin_nu_factor_axisym,
)


def test_3d_nu_and_mu_factors_are_reciprocal_and_match_at_boundary():
    R = 0.2
    center = np.array([0.1, -0.2, 0.3])
    pts = np.array([
        center + [R, 0.0, 0.0],
        center + [0.0, 2 * R, 0.0],
        center + [0.0, 0.0, 0.5 * R],
    ])

    nu = kelvin_nu_factor_3d(pts, center, R)
    mu = kelvin_mu_factor_3d(pts, center, R)

    np.testing.assert_allclose(nu * mu, 1.0)
    assert np.isclose(nu[0], 1.0)
    assert np.isclose(mu[0], 1.0)
    assert np.isclose(nu[1], 4.0)
    assert np.isclose(mu[2], 4.0)


def test_axisymmetric_uses_3d_spherical_distance_and_reciprocal_pair():
    R = 0.05
    z_offset = 5 * R
    r = np.array([R, 0.0, 2 * R])
    z = np.array([z_offset, z_offset + R, z_offset])

    nu = kelvin_nu_factor_axisym(r, z, z_offset, R)
    mu = kelvin_mu_factor_axisym(r, z, z_offset, R)

    np.testing.assert_allclose(nu, [1.0, 1.0, 4.0])
    np.testing.assert_allclose(nu * mu, 1.0)


def test_2d_inplane_is_identity_and_axial_slot_is_fourth_power():
    R = 0.4
    pts = np.array([[R, 0.0], [2 * R, 0.0], [0.0, 0.5 * R]])

    assert kelvin_factor_2d_inplane() == 1.0
    np.testing.assert_allclose(kelvin_nu_factor_2d_axial(pts, (0.0, 0.0), R),
                               [1.0, 16.0, 1.0 / 16.0])
    np.testing.assert_allclose(kelvin_mu_factor_2d_axial(pts, (0.0, 0.0), R),
                               [1.0, 1.0 / 16.0, 16.0])


def test_kelvin_center_regularization_is_explicit_and_finite():
    center = np.zeros(3)
    R = 1.0
    pt = np.zeros(3)

    assert np.isfinite(kelvin_nu_factor_3d(pt, center, R, rho_min=0.25))
    assert np.isfinite(kelvin_mu_factor_3d(pt, center, R, rho_min=0.25))
    assert np.isclose(kelvin_nu_factor_3d(pt, center, R, rho_min=0.25), 0.25 ** 2)
    assert np.isclose(kelvin_mu_factor_3d(pt, center, R, rho_min=0.25), 4.0 ** 2)
