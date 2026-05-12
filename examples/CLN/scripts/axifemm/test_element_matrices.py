"""
Unit tests for axifemm_core element matrices.

Validates that:
1. Mr and Mz are symmetric.
2. For an interior triangle (no axis), formulas give sensible scale.
3. Element matrices satisfy "constant solution -> zero residual" in appropriate sense.
4. Axis-touching cases (flag=1, flag=2) execute without NaN/Inf.
"""

from __future__ import annotations

import numpy as np

from axifemm_core import element_matrices, _element_geom, MU0


def test_symmetry_interior():
    rn = np.array([1.0e-3, 3.0e-3, 2.0e-3])
    zn = np.array([0.0, 0.0, 1.5e-3])
    M_e, Mr, Mz, *_ = element_matrices(rn, zn, MU0, MU0)
    assert np.allclose(Mr, Mr.T), "Mr not symmetric"
    assert np.allclose(Mz, Mz.T), "Mz not symmetric"
    print("[OK] interior triangle Mr, Mz are symmetric")


def test_symmetry_axis_2nodes():
    rn = np.array([0.0, 1.0e-3, 0.0])
    zn = np.array([0.0, 0.5e-3, 1.0e-3])
    M_e, Mr, Mz, *_ = element_matrices(rn, zn, MU0, MU0)
    assert np.allclose(Mr, Mr.T)
    assert np.allclose(Mz, Mz.T)
    assert np.all(np.isfinite(M_e))
    print("[OK] flag=2 (2 nodes on axis) executes without NaN/Inf")


def test_symmetry_axis_1node():
    rn = np.array([0.0, 1.0e-3, 1.5e-3])
    zn = np.array([0.0, 0.0, 1.0e-3])
    M_e, Mr, Mz, *_ = element_matrices(rn, zn, MU0, MU0)
    assert np.allclose(Mr, Mr.T)
    assert np.allclose(Mz, Mz.T)
    assert np.all(np.isfinite(M_e))
    print("[OK] flag=1 (1 node on axis) executes without NaN/Inf")


def test_interpolation_consistency():
    """
    For a triangle with phi = c0 + c1*r^2 + c2*z (the {1,r^2,z} basis itself),
    the resulting B_z should be c1/pi (constant) and B_r should be -c2/(2*pi*r).

    We verify the postprocess routine matches.
    """
    from axifemm_core import field_in_element

    rn = np.array([1.0e-3, 3.0e-3, 2.0e-3])
    zn = np.array([0.0, 0.0, 1.5e-3])
    c0, c1, c2 = 1.5, 2.7, 3.1
    phi_e = c0 + c1 * rn ** 2 + c2 * zn

    # evaluate at centroid
    r_c = rn.mean()
    z_c = zn.mean()
    A_phi, B_r, B_z = field_in_element(rn, zn, phi_e, r_c, z_c)

    expected_Bz = c1 / np.pi
    expected_Br = -c2 / (2 * np.pi * r_c)
    assert abs(B_z - expected_Bz) < 1e-10, f"B_z mismatch: {B_z} vs {expected_Bz}"
    assert abs(B_r - expected_Br) < 1e-10, f"B_r mismatch: {B_r} vs {expected_Br}"
    print(f"[OK] interpolation: B_z={B_z:.6f} B_r={B_r:.6f}")


def test_const_phi_zero_energy():
    """
    Constant flux phi field => B_z = 0 and B_r = 0 => zero energy.
    Note: prob3big.cpp solves in V (= A_phi/mu_0)-form; constant phi corresponds
    to V_j ~ 1/r_j (since phi = 2*pi*r*mu_0*V). So we test that the V-form
    matrices annihilate the (1/r_0, 1/r_1, 1/r_2) vector.
    """
    rn = np.array([1.0e-3, 3.0e-3, 2.0e-3])
    zn = np.array([0.0, 0.0, 1.5e-3])
    M_e, Mr, Mz, *_ = element_matrices(rn, zn, MU0, MU0)
    v_const_phi = 1.0 / rn  # V-form of constant phi
    res_Mr = Mr @ v_const_phi
    res_Mz = Mz @ v_const_phi
    print(f"  Mr * (1/r_j) = {res_Mr}")
    print(f"  Mz * (1/r_j) = {res_Mz}")
    assert np.allclose(res_Mr, 0, atol=1e-10), "Mr * (1/r) != 0"
    assert np.allclose(res_Mz, 0, atol=1e-10), "Mz * (1/r) != 0"
    print("[OK] const phi (V-form: 1/r) annihilates Mr, Mz")


def test_geom_consistency():
    """
    Verify a_hat = d_meeker / (4*R) where d_meeker = sum_j s_j * p_j
    and the standard P1 area is a = (p[0]*q[1] - p[1]*q[0])/2.
    """
    rn = np.array([1.0e-3, 3.0e-3, 2.0e-3])
    zn = np.array([0.0, 0.0, 1.5e-3])
    p, q, g, a, R, a_hat, R_hat = _element_geom(rn, zn)
    s = rn ** 2
    d_meeker = sum(s[j] * p[j] for j in range(3))
    a_hat_check = d_meeker / (4 * R)
    print(f"  d_meeker = {d_meeker:.6e}, a_hat = {a_hat:.6e}, check = {a_hat_check:.6e}")
    assert abs(a_hat - a_hat_check) < 1e-15
    print("[OK] a_hat = d_meeker / (4*R)")


if __name__ == "__main__":
    test_geom_consistency()
    test_symmetry_interior()
    test_symmetry_axis_2nodes()
    test_symmetry_axis_1node()
    test_interpolation_consistency()
    test_const_phi_zero_energy()
    print("\nAll tests passed.")
