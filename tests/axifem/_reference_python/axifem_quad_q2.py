"""
axifem_quad_q2.py — Phase 2-F (axifem Q2 prototype)

Q2 axisymmetric Henrotte element on axis-aligned rectangle.

Basis (9 monomials in s = r^2, z):
    phi(r, z) = a + b s + c s^2
              + d z + e s z + f s^2 z
              + g z^2 + h s z^2 + i s^2 z^2

9 nodes (3x3 Lagrange grid in (s, z)):
    s_mid = (s_a + s_b)/2      (NOT (r_a+r_b)^2/4)
    z_mid = (z_a + z_b)/2

Physical r-coord at edge midpoints: r_mid = sqrt(s_mid) = sqrt((r_a^2+r_b^2)/2)
(quadratic mean of r_a, r_b — NOT arithmetic mean).

Field components:
    B_z = (1/(2 pi r)) d phi/dr = d phi/ds / pi
    B_r = -(1/(2 pi r)) d phi/dz = -(d phi/dz) / (2 pi sqrt(s))

Energy (Hessian-of-W convention, matching axifem_quad.py / Q1 reference):
    W_z = Integrate( (d phi/ds)^2 / (pi mu_z), {s, sa, sb}, {z, za, zb} )
    W_r = Integrate( (d phi/dz)^2 / (4 pi mu_r s), {s, sa, sb}, {z, za, zb} )
    M_sig contribution: sigma/(4 pi s) * phi^2 in (s, z) box
    (NOT sigma/(8 pi^2) -- that's the W coefficient, 1/(2 pi) too small.)

This prototype uses high-order Gauss-Legendre quadrature in (s, z). It serves
as a reference for validating the Mathematica closed-form derivation and as a
fallback if Mathematica is too slow.
"""
from __future__ import annotations

import numpy as np

# 8-point Gauss-Legendre on [-1, 1]: exact for polynomials up to degree 15.
# Sufficient for Q2 W_z (degree 8 in s, 4 in z) and for log-substituted W_r.
_GL8_x, _GL8_w = np.polynomial.legendre.leggauss(8)


def _phi_basis_q2(s, z):
    """Return 9-vector of monomial values [1, s, s^2, z, s z, s^2 z, z^2, s z^2, s^2 z^2]."""
    return np.array([
        1.0, s, s*s,
        z, s*z, s*s*z,
        z*z, s*z*z, s*s*z*z,
    ])


def _phi_dPhiDs(s, z):
    """d/ds of monomials: [0, 1, 2 s, 0, z, 2 s z, 0, z^2, 2 s z^2]"""
    return np.array([
        0.0, 1.0, 2.0 * s,
        0.0, z, 2.0 * s * z,
        0.0, z*z, 2.0 * s * z*z,
    ])


def _phi_dPhiDz(s, z):
    """d/dz of monomials: [0, 0, 0, 1, s, s^2, 2 z, 2 s z, 2 s^2 z]"""
    return np.array([
        0.0, 0.0, 0.0,
        1.0, s, s*s,
        2.0 * z, 2.0 * s * z, 2.0 * s * s * z,
    ])


def _node_coords(ra, rb, za, zb):
    """Return 9 (s, z) node coords for the Q2 axis-aligned rectangle."""
    sa, sb = ra*ra, rb*rb
    sm = 0.5 * (sa + sb)
    zm = 0.5 * (za + zb)
    return np.array([
        [sa, za],   # 0 corner
        [sb, za],   # 1 corner
        [sb, zb],   # 2 corner
        [sa, zb],   # 3 corner
        [sm, za],   # 4 edge bottom
        [sb, zm],   # 5 edge right
        [sm, zb],   # 6 edge top
        [sa, zm],   # 7 edge left
        [sm, zm],   # 8 face center
    ])


def _node_r(ra, rb):
    """Return r-coord (NOT s-coord) at each of the 9 nodes."""
    rm = np.sqrt(0.5 * (ra*ra + rb*rb))   # quadratic mean
    return np.array([ra, rb, rb, ra, rm, rb, rm, ra, rm])


def _vandermonde_q2(ra, rb, za, zb):
    """9x9 Vandermonde V[k, j] = phi_j evaluated at node k."""
    nodes = _node_coords(ra, rb, za, zb)
    V = np.zeros((9, 9))
    for k in range(9):
        s_k, z_k = nodes[k]
        V[k] = _phi_basis_q2(s_k, z_k)
    return V


def element_matrices_q2_numerical(ra, rb, za, zb, mu_r, mu_z):
    """Return (M_V, M_phi) for the Q2 axisymmetric Henrotte element.

    M_phi: 9x9 stiffness in phi (= 2 pi r A_phi) coords.
    M_V:   9x9 stiffness in V (= A_phi nodal) coords.

    Computed by numerical Gauss-Legendre quadrature on (s, z) box.
    For axis-touching elements (ra = 0), the 1/s integrand has a log
    singularity that is integrated analytically by splitting:
        Integrate( P(s, z) / s, {s, 0, sb} )
        = Integrate( (P(s, z) - P(0, z)) / s, {s, 0, sb} )  + P(0, z) * Limit log s
    The log term is divergent; rely on V-DOF transformation to cancel it.

    For ra > 0, no special treatment is needed — Gauss handles 1/s smoothly.
    """
    sa, sb = ra*ra, rb*rb
    if sa <= 0:
        raise NotImplementedError(
            "Axis-touching case (ra=0) not yet implemented in numerical path."
            " Use Mathematica closed-form for axis elements.")

    # Gauss-Legendre on (s, z) = transform of [-1, 1] x [-1, 1].
    # s = (sa + sb)/2 + ((sb - sa)/2) xi,  z similarly.
    s_jac = 0.5 * (sb - sa)
    z_jac = 0.5 * (zb - za)
    s_mid_phys = 0.5 * (sa + sb)
    z_mid_phys = 0.5 * (za + zb)

    # In (s, z), each monomial in phi gives partial derivatives.
    # K_phi[i, j] = Integrate( dPhi_i/ds dPhi_j/ds / (2 pi mu_z)
    #                        + dPhi_i/dz dPhi_j/dz / (8 pi mu_r s) ) ds dz
    Kphi = np.zeros((9, 9))
    for i, xi in enumerate(_GL8_x):
        s = s_mid_phys + s_jac * xi
        for j, eta in enumerate(_GL8_x):
            z = z_mid_phys + z_jac * eta
            w = _GL8_w[i] * _GL8_w[j] * s_jac * z_jac
            ds = _phi_dPhiDs(s, z)
            dz = _phi_dPhiDz(s, z)
            Kphi += w * (
                np.outer(ds, ds) / (np.pi * mu_z) +
                np.outer(dz, dz) / (4.0 * np.pi * mu_r * s)
            )

    # Convert to nodal (phi-DOF) basis. M_node = (V^-1)^T . K_phi . V^-1
    V = _vandermonde_q2(ra, rb, za, zb)
    Vinv = np.linalg.inv(V)
    Mphi_nodal = Vinv.T @ Kphi @ Vinv

    # Convert phi-DOF -> V-DOF: V_j = phi_j / (2 pi r_j)
    rj = _node_r(ra, rb)
    T = 2 * np.pi * rj   # diagonal
    M_V = (T[:, None] * Mphi_nodal) * T[None, :]

    return M_V, Mphi_nodal


def element_sigma_mass_q2_numerical(ra, rb, za, zb, sigma):
    """Return 9x9 sigma mass matrix M_sigma_V in V coords.

    W_sigma = (sigma/2) Integrate( A_phi^2, 2 pi r dr dz )
            = (sigma/(8 pi^2)) Integrate( phi^2 / s ds dz )
    """
    sa, sb = ra*ra, rb*rb
    if sa <= 0:
        raise NotImplementedError(
            "Axis-touching case (ra=0) not yet implemented in numerical path.")

    s_jac = 0.5 * (sb - sa)
    z_jac = 0.5 * (zb - za)
    s_mid_phys = 0.5 * (sa + sb)
    z_mid_phys = 0.5 * (za + zb)

    Mphi = np.zeros((9, 9))
    for i, xi in enumerate(_GL8_x):
        s = s_mid_phys + s_jac * xi
        for j, eta in enumerate(_GL8_x):
            z = z_mid_phys + z_jac * eta
            w = _GL8_w[i] * _GL8_w[j] * s_jac * z_jac
            phi = _phi_basis_q2(s, z)
            Mphi += w * sigma / (4.0 * np.pi * s) * np.outer(phi, phi)

    V = _vandermonde_q2(ra, rb, za, zb)
    Vinv = np.linalg.inv(V)
    Mphi_nodal = Vinv.T @ Mphi @ Vinv

    rj = _node_r(ra, rb)
    T = 2 * np.pi * rj
    M_V = (T[:, None] * Mphi_nodal) * T[None, :]
    return M_V


if __name__ == "__main__":
    ra, rb, za, zb = 1e-3, 2e-3, 0.0, 1e-3
    mu0 = 4 * np.pi * 1e-7
    sigma_Cu = 5.8e7

    K, Kphi = element_matrices_q2_numerical(ra, rb, za, zb, mu0, mu0)
    M = element_sigma_mass_q2_numerical(ra, rb, za, zb, sigma_Cu)

    print(f"=== Q2 Henrotte element (numerical Gauss 8x8) ===")
    print(f"ra={ra*1e3:.1f}mm rb={rb*1e3:.1f}mm za={za*1e3:.1f}mm zb={zb*1e3:.1f}mm")
    print()
    print(f"K symmetric? max|K-K^T| = {np.max(np.abs(K - K.T)):.3e}")
    print(f"M symmetric? max|M-M^T| = {np.max(np.abs(M - M.T)):.3e}")
    print(f"K eigenvalues (top 4): {np.sort(np.linalg.eigvalsh(K))[-4:]}")
    print(f"M eigenvalues (top 4): {np.sort(np.linalg.eigvalsh(M))[-4:]}")
    print()
    print("K[:3, :3]:")
    print(K[:3, :3])
    print("M[:3, :3]:")
    print(M[:3, :3])
