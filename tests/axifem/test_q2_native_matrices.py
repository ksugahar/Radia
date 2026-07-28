"""Independent gates for the shared-native Q2 Henrotte element matrices."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


sys.path.insert(0, str(Path(__file__).parent / "_reference_python"))
import axifem_quad_q2 as q2_reference  # noqa: E402


axifem = pytest.importorskip("radia.axifem")


def test_q2_interior_matches_independent_gauss_reference():
    mu0 = 4.0 * np.pi * 1.0e-7
    sigma = 5.8e7
    values = (1.0e-3, 2.0e-3, -0.5e-3, 0.5e-3, mu0, sigma)

    actual = axifem.q2_magnetic_element_matrices(*values)
    # The radial 1/s terms are rational, so the historical 8-point teaching
    # reference is not a machine-precision oracle. Use a denser independent
    # quadrature here when checking the generated closed form.
    q2_reference._GL8_x, q2_reference._GL8_w = (  # noqa: SLF001
        np.polynomial.legendre.leggauss(64)
    )
    expected_k, _ = q2_reference.element_matrices_q2_numerical(
        *values[:4], values[4], values[4]
    )
    expected_m = q2_reference.element_sigma_mass_q2_numerical(
        *values[:4], values[5]
    )

    np.testing.assert_allclose(
        actual["stiffness"], expected_k, rtol=3.0e-10, atol=3.0e-9
    )
    np.testing.assert_allclose(
        actual["sigma_mass"], expected_m, rtol=1.0e-11, atol=2.0e-13
    )
    assert actual["axis_touching"] is False


def test_q2_axis_rows_are_explicitly_eliminated():
    mu0 = 4.0 * np.pi * 1.0e-7
    actual = axifem.q2_magnetic_element_matrices(
        0.0, 1.0e-3, -0.5e-3, 0.5e-3, mu0, 5.8e7
    )

    stiffness = np.asarray(actual["stiffness"])
    sigma_mass = np.asarray(actual["sigma_mass"])
    assert actual["axis_touching"] is True
    assert np.isfinite(stiffness).all()
    assert np.isfinite(sigma_mass).all()
    np.testing.assert_allclose(stiffness, stiffness.T, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(sigma_mass, sigma_mass.T, rtol=0.0, atol=0.0)
    for index in (0, 3, 7):
        np.testing.assert_array_equal(stiffness[index, :], 0.0)
        np.testing.assert_array_equal(stiffness[:, index], 0.0)
        np.testing.assert_array_equal(sigma_mass[index, :], 0.0)
        np.testing.assert_array_equal(sigma_mass[:, index], 0.0)


@pytest.mark.parametrize(
    "values",
    [
        (2.0e-3, 1.0e-3, 0.0, 1.0e-3, 1.0, 1.0),
        (1.0e-3, 2.0e-3, 1.0e-3, 0.0, 1.0, 1.0),
        (1.0e-3, 2.0e-3, 0.0, 1.0e-3, 0.0, 1.0),
        (1.0e-3, 2.0e-3, 0.0, 1.0e-3, 1.0, -1.0),
    ],
)
def test_q2_invalid_inputs_fail_loudly(values):
    with pytest.raises(ValueError):
        axifem.q2_magnetic_element_matrices(*values)
