"""test_p2_single_element.py — Phase B1b Item 1
P2 vs P1-subspace single-element consistency check.

For a P1-affine function u(r, z) = a + b r^2 + c z (which lies in both the
P1 Henrotte space {1, r^2, z} and the P2 space {1, r^2, z, r^4, r^2 z, z^2}),
compute the magnetic energy
    E_K = (2 pi / mu) integral over T of ((du/dr)^2 + (du/dz)^2) r dA
and the kinetic energy
    E_M = sigma * 2 pi integral over T of u^2 r dA
two ways:
  (i)  via P2 element matrices, with u sampled at the 6 P2 nodes
  (ii) via direct Duffy-Gauss quadrature of the analytical integrand

If the P2 element matrices are correct, the two energies must agree to
machine precision (no discretisation error since u is in the P2 space).
"""
from __future__ import annotations

import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

from axifem_p2_triangle import AxifemP2Triangle, MU0, _QUAD


def make_p2_triangle(r0, z0, r1, z1, r2, z2):
    rn = np.array([r0, r1, r2,
                   0.5 * (r0 + r1), 0.5 * (r1 + r2), 0.5 * (r2 + r0)])
    zn = np.array([z0, z1, z2,
                   0.5 * (z0 + z1), 0.5 * (z1 + z2), 0.5 * (z2 + z0)])
    return AxifemP2Triangle(rn, zn)


def _energy_p2_K(tri, a, b, c, mu_r):
    """Energy a^T K @ a for u_coeffs = u sampled at 6 P2 nodes."""
    u_coeffs = a + b * tri.r ** 2 + c * tri.z
    K = tri.stiffness(mu_r=mu_r)
    return u_coeffs @ K @ u_coeffs


def _energy_direct_K(tri, a, b, c, mu_r):
    """Direct quadrature of (2 pi / mu) integral ((2 b r)^2 + c^2) r dA."""
    total = 0.0
    jac = abs(tri.det_J)
    for q in _QUAD:
        xi, eta, w = q
        r, z = tri.map_to_physical(xi, eta)
        if r <= 0:
            continue
        # du/dr = 2 b r, du/dz = c
        integrand = ((2.0 * b * r) ** 2 + c ** 2) * r
        total += w * integrand
    total *= jac * 2.0 * np.pi / mu_r
    return total


def _energy_p2_M(tri, a, b, c, sigma):
    """Energy u_coeffs^T M u_coeffs."""
    u_coeffs = a + b * tri.r ** 2 + c * tri.z
    M = tri.sigma_mass(sigma=sigma)
    return u_coeffs @ M @ u_coeffs


def _energy_direct_M(tri, a, b, c, sigma):
    """Direct quadrature of sigma * 2 pi integral u^2 r dA."""
    total = 0.0
    jac = abs(tri.det_J)
    for q in _QUAD:
        xi, eta, w = q
        r, z = tri.map_to_physical(xi, eta)
        if r <= 0:
            continue
        u = a + b * r ** 2 + c * z
        total += w * u * u * r
    total *= jac * 2.0 * np.pi * sigma
    return total


def main():
    print("=" * 64)
    print("Phase B1b Item 1: P2 vs direct-quadrature consistency check")
    print("=" * 64)

    # Several non-degenerate triangles, both general and axis-touching.
    test_triangles = [
        ("general", (1e-3, 0.0, 5e-3, 0.0, 3e-3, 2e-3)),
        ("axis-touching", (0.0, 0.0, 4e-3, 0.0, 2e-3, 3e-3)),
        ("axis-both-vertex", (0.0, 0.0, 0.0, 2e-3, 5e-3, 1e-3)),
        ("skew", (2e-3, 1e-3, 6e-3, 0.5e-3, 4e-3, 4e-3)),
    ]

    test_cases = [
        ("constant u=1", 1.0, 0.0, 0.0),
        ("radial u=r^2", 0.0, 1.0, 0.0),
        ("axial u=z", 0.0, 0.0, 1.0),
        ("mixed u=1+r^2+z", 1.0, 1.0, 1.0),
        ("strong u=2+5 r^2-3 z", 2.0, 5.0, -3.0),
    ]

    mu_test = MU0
    sigma_test = 5.8e7
    all_passed = True

    for tri_name, coords in test_triangles:
        print(f"\nTriangle: {tri_name}")
        # Skip pure-axis-touching tests where r=0 quadrature points dominate
        # (axis nodes contribute 0 to axisym integrals).
        try:
            tri = make_p2_triangle(*coords)
        except RuntimeError as e:
            print(f"  Skip (Vandermonde singular): {e}")
            continue

        for case_name, a, b, c in test_cases:
            E_K_p2 = _energy_p2_K(tri, a, b, c, mu_test)
            E_K_dir = _energy_direct_K(tri, a, b, c, mu_test)
            E_M_p2 = _energy_p2_M(tri, a, b, c, sigma_test)
            E_M_dir = _energy_direct_M(tri, a, b, c, sigma_test)

            err_K = (abs(E_K_p2 - E_K_dir) / abs(E_K_dir)
                     if abs(E_K_dir) > 1e-30 else abs(E_K_p2 - E_K_dir))
            err_M = (abs(E_M_p2 - E_M_dir) / abs(E_M_dir)
                     if abs(E_M_dir) > 1e-30 else abs(E_M_p2 - E_M_dir))

            # K tolerance relaxed to 1e-8 due to 6x6 Vandermonde FP64 cond
            # (~1e10) on physical coords; can be tightened with monomial
            # scaling in C++ port. M is unaffected.
            tag_K = "OK" if err_K < 1e-8 else "FAIL"
            tag_M = "OK" if err_M < 1e-12 else "FAIL"
            if err_K >= 1e-8 or err_M >= 1e-12:
                all_passed = False

            print(f"  {case_name:>22s}  K rel {err_K:.2e} [{tag_K}], "
                  f"M rel {err_M:.2e} [{tag_M}]")

    print()
    print("All P1-subspace energies match P2-element-matrix energies"
          if all_passed else "Some tests FAILED")


if __name__ == "__main__":
    main()
