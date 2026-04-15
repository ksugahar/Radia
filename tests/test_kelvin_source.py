"""Sanity tests for kelvin_source.biot_savart_A_at_points + Kelvin map.

Compares Biot-Savart A_z at the axis of a single circular loop
with the analytical formula (on axis only, since this is where the
closed-form is simple):

    A_theta(rho=R, z=z0) analytical via elliptic integrals -- we
    simplify by testing A at the CENTER of a single small loop which
    is zero by symmetry except for the self-induced component.

More useful sanity check: compare two different discretizations of
the same coil for self-consistency.
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src', 'radia'))

MU_0 = 4e-7 * math.pi


def build_loop(R, z0, n, current):
    """Discretize a circular loop in xy-plane at z=z0 into n segments."""
    thetas = np.linspace(0, 2 * math.pi, n + 1)
    pts = np.column_stack([
        R * np.cos(thetas), R * np.sin(thetas),
        np.full_like(thetas, z0)])
    segs = [(pts[i].tolist(), pts[i + 1].tolist()) for i in range(n)]
    return [segs], [current]


def test_on_axis_B_circular_loop():
    """Analytical B_z on axis of a single circular loop:

        B_z(0, 0, z) = mu_0 * I * R^2 / (2 * (R^2 + z^2)^{3/2})

    We compute A at (eps, 0, z) and take finite-difference curl to
    get B_z ~ (A_theta_at_z+dz - A_theta_at_z-dz) / (2 dz) -- too
    fragile. Instead: just check that A_phi scales ~ R^2/(R^2+z^2)^{1/2}
    magnitude.
    """
    from kelvin_source import biot_savart_A_at_points

    R = 0.010
    I = 1.0
    paths, currents = build_loop(R, 0.0, 200, I)
    # Evaluate A on axis at several z
    for z in [0.0, 0.005, 0.010, 0.020]:
        obs = np.array([[1e-6, 0.0, z]])
        A = biot_savart_A_at_points(obs, paths, currents, n_quad=6)
        # On axis, dipole expansion for small rho:
        #   A_phi(rho, z) ~ mu_0 m * rho / (4 pi (R^2 + z^2)^{3/2})
        # where m = I * pi * R^2 (magnetic moment), so
        #   A_phi ~ (mu_0 I / 4) * rho * R^2 / (R^2 + z^2)^{3/2}
        rho = 1e-6
        expected_A_y = (MU_0 * I / 4.0) * rho * R ** 2 \
            / (R ** 2 + z ** 2) ** 1.5
        A_y = A[0, 1]
        rel_err = abs(A_y - expected_A_y) / abs(expected_A_y)
        print(f"  z={z*1e3:5.1f}mm   A_y={A_y:.4e}   "
              f"expected={expected_A_y:.4e}   rel_err={rel_err*100:.3f}%")
        assert rel_err < 0.02, f"Loop A at z={z} off by {rel_err*100:.2f}%"


def test_kelvin_map_involution():
    from kelvin_source import kelvin_map_3d

    center = np.array([0.0, 0.0, 0.0])
    R = 1.0
    rng = np.random.default_rng(42)
    pts = rng.normal(size=(10, 3))
    pts[:, :] *= 2 + np.linalg.norm(pts, axis=1)[:, None]
    mapped = kelvin_map_3d(pts, center, R)
    double = kelvin_map_3d(mapped, center, R)
    err = np.max(np.abs(double - pts))
    print(f"  max involution error: {err:.2e}")
    assert err < 1e-10


def test_kelvin_map_on_sphere():
    """Points on the Kelvin sphere map to themselves."""
    from kelvin_source import kelvin_map_3d

    center = np.array([0.1, -0.2, 0.3])
    R = 0.5
    thetas = np.linspace(0, 2 * math.pi, 12, endpoint=False)
    pts = center + R * np.column_stack([np.cos(thetas),
                                         np.sin(thetas),
                                         np.zeros_like(thetas)])
    mapped = kelvin_map_3d(pts, center, R)
    err = np.max(np.abs(mapped - pts))
    print(f"  on-sphere invariance error: {err:.2e}")
    assert err < 1e-12


def test_scalar_factor_at_sphere_boundary():
    """Kelvin factor (R/rho')^2 = 1 at rho' = R."""
    from kelvin_source import kelvin_factor_scalar

    center = np.array([0.0, 0.0, 0.0])
    R = 0.5
    pts = np.array([[R, 0, 0], [0, R, 0], [0, 0, R], [-R, 0, 0]])
    f = kelvin_factor_scalar(pts, center, R)
    print(f"  factors at rho'=R: {f}")
    assert np.allclose(f, 1.0)


def test_A_s_kelvin_pathway():
    """A_s_at_obs_with_kelvin: inner points unchanged, shell points mapped."""
    from kelvin_source import A_s_at_obs_with_kelvin, biot_savart_A_at_points

    R_loop = 0.030
    R_kelvin = 0.060
    offset = np.array([0.150, 0.0, 0.0])   # Kelvin sphere center (Sugahara 2022)
    I = 1.0
    paths, currents = build_loop(R_loop, 0.0, 200, I)
    # Inner point (physical inner domain, far from Kelvin sphere)
    p_in = np.array([0.010, 0.0, 0.005])
    # Kelvin shell point (computational frame): offset +/- dx within radius
    p_out_shell = offset + np.array([0.020, 0.005, 0.010])
    # Its Kelvin-inverse physical position:
    from kelvin_source import kelvin_map_3d
    p_out_phys = kelvin_map_3d(p_out_shell, offset, R_kelvin)

    obs = np.vstack([p_in, p_out_shell])
    A_full = A_s_at_obs_with_kelvin(
        obs, paths, currents,
        kelvin_center=offset.tolist(), R_kelvin=R_kelvin,
        factor_mode='scalar')

    # Inner vs direct
    A_in_direct = biot_savart_A_at_points(p_in.reshape(1, 3),
                                           paths, currents, n_quad=8)[0]
    err_in = np.max(np.abs(A_full[0] - A_in_direct))
    print(f"  inner A match: err = {err_in:.2e}")
    assert err_in < 1e-12

    # Shell: expected = A(physical_equivalent) * (R/rho')^2
    A_at_phys = biot_savart_A_at_points(p_out_phys.reshape(1, 3),
                                         paths, currents, n_quad=8)[0]
    rho_prime = np.linalg.norm(p_out_shell - offset)
    factor = (R_kelvin / rho_prime) ** 2
    expected = A_at_phys * factor
    err_sh = np.max(np.abs(A_full[1] - expected))
    rel_sh = err_sh / np.max(np.abs(expected))
    print(f"  shell factor match: err = {err_sh:.2e}  rel = {rel_sh:.2e}")
    # Quadrature-level tolerance (loop discretization + n_quad=8).
    assert rel_sh < 1e-6


if __name__ == "__main__":
    print("test_on_axis_B_circular_loop"); test_on_axis_B_circular_loop()
    print("test_kelvin_map_involution"); test_kelvin_map_involution()
    print("test_kelvin_map_on_sphere"); test_kelvin_map_on_sphere()
    print("test_scalar_factor_at_sphere_boundary");
    test_scalar_factor_at_sphere_boundary()
    print("test_A_s_kelvin_pathway"); test_A_s_kelvin_pathway()
    print("\nall tests passed.")
