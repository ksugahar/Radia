"""
sigma_mass.py — Sigma mass matrix for axifem (V-form, FEMM-style basis).

For Hiruma 3-term CLN we need:
  M_sigma_ij = sigma * Integrate( A_phi(r,z) * v_phi(r,z) * 2*pi*r dr dz )  over conductor

In our convention:
- We solve for V_j at each node where V_j = A_phi at node j (= phi_j / (2*pi*r_j))
- phi(r,z) = c0 + c1*r^2 + c2*z is linear in (s=r^2, z)
- A_phi(r,z) = phi(r,z) / (2*pi*r) is NON-linear in (r,z), but well-behaved (no singularity since phi ~ r at axis)

So per element with nodal values (V_0, V_1, V_2):
- phi_j = 2*pi*r_j*V_j
- Solve [1, r_j^2, z_j] coefficients to get (c0, c1, c2) in terms of (phi_0, phi_1, phi_2)
- A_phi(r,z) = (c0 + c1*r^2 + c2*z) / (2*pi*r)
- M_e[i,j] V_i V_j = sigma * Integrate( A_phi^2 * 2*pi*r dr dz )

Implementation: numerical quadrature on each triangle. Use Hammer 7-point rule (exact for P5).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


# Hammer 7-point Gauss quadrature on the reference triangle (vertices (0,0),(1,0),(0,1))
# Exact for polynomials up to degree 5.
_S15 = np.sqrt(15.0)
_HAMMER7 = np.array([
    [1.0 / 3.0,                  1.0 / 3.0,                  9.0 / 80.0],
    [(6.0 + _S15) / 21.0,        (6.0 + _S15) / 21.0,        (155.0 - _S15) / 2400.0],
    [(9.0 - 2.0 * _S15) / 21.0,  (6.0 + _S15) / 21.0,        (155.0 - _S15) / 2400.0],
    [(6.0 + _S15) / 21.0,        (9.0 - 2.0 * _S15) / 21.0,  (155.0 - _S15) / 2400.0],
    [(6.0 - _S15) / 21.0,        (6.0 - _S15) / 21.0,        (155.0 + _S15) / 2400.0],
    [(9.0 + 2.0 * _S15) / 21.0,  (6.0 - _S15) / 21.0,        (155.0 + _S15) / 2400.0],
    [(6.0 - _S15) / 21.0,        (9.0 + 2.0 * _S15) / 21.0,  (155.0 + _S15) / 2400.0],
])


def element_sigma_mass(rn: np.ndarray, zn: np.ndarray, sigma: float) -> np.ndarray:
    """
    Per-element 3x3 sigma mass matrix in V coordinates.

    M_e[i,j] = sigma * Integrate_T( N_A_i(r,z) * N_A_j(r,z) * 2*pi*r dr dz )

    where N_A_j is the A_phi-basis function for node j:
      N_A_j(r,z) = (alpha_j + beta_j*r^2 + gamma_j*z) * 2*pi*r_j / (2*pi*r)
                 = r_j * (alpha_j + beta_j*r^2 + gamma_j*z) / r

    The (alpha, beta, gamma) coefficients come from inverting the Vandermonde
    [1, r_j^2, z_j] -> Kronecker delta system.

    Uses 7-point Hammer quadrature in physical (r,z) space.
    """
    # Vandermonde for {1, r^2, z} basis on the triangle
    M_vand = np.column_stack([np.ones(3), rn ** 2, zn])
    # Solve M_vand @ coeffs = I  =>  coeffs[k, j] = k-th coefficient of phi-basis for node j
    # Then phi-basis: psi_j(r,z) = coeffs[0,j] + coeffs[1,j]*r^2 + coeffs[2,j]*z
    # Such that psi_j(r_k, z_k) = delta_jk
    try:
        inv_vand = np.linalg.inv(M_vand)
    except np.linalg.LinAlgError:
        # Degenerate triangle
        return np.zeros((3, 3))
    # inv_vand has rows = (alpha, beta, gamma) basis coefficients
    # column j gives the (alpha_j, beta_j, gamma_j) for the j-th node basis function
    alpha = inv_vand[0, :]  # (3,)
    beta = inv_vand[1, :]
    gamma = inv_vand[2, :]

    # Affine map ref triangle (xi, eta) -> physical (r, z)
    # r = r_0 + (r_1 - r_0)*xi + (r_2 - r_0)*eta
    # z = z_0 + (z_1 - z_0)*xi + (z_2 - z_0)*eta
    # |Jacobian| = |(r_1 - r_0)*(z_2 - z_0) - (r_2 - r_0)*(z_1 - z_0)| = 2 * |area|
    drxi = rn[1] - rn[0]
    dreta = rn[2] - rn[0]
    dzxi = zn[1] - zn[0]
    dzeta = zn[2] - zn[0]
    detJ = abs(drxi * dzeta - dreta * dzxi)  # = 2 * triangle area in (r,z) space

    M_e = np.zeros((3, 3))
    for xi, eta, w in _HAMMER7:
        r_q = rn[0] + drxi * xi + dreta * eta
        z_q = zn[0] + dzxi * xi + dzeta * eta
        if r_q <= 0:
            # Should not happen with interior Gauss points, but skip if so
            continue
        # phi basis values: psi_j(r_q, z_q) = alpha_j + beta_j*r_q^2 + gamma_j*z_q
        psi = alpha + beta * (r_q ** 2) + gamma * z_q  # (3,)
        # A_phi-basis: N_A_j = psi_j / (2*pi*r_q) * (2*pi*r_j) = r_j * psi_j / r_q
        # In V coordinates: V_j corresponds to A_phi at node j
        N_A = rn * psi / r_q  # (3,)
        # Integrand at quad point: sigma * N_A_i * N_A_j * 2*pi*r_q
        integrand_contrib = sigma * np.outer(N_A, N_A) * 2.0 * np.pi * r_q
        M_e += integrand_contrib * w * detJ
    return M_e


def assemble_sigma_mass(
    nodes: np.ndarray,
    triangles: np.ndarray,
    sigma_per_elem: np.ndarray,
) -> sp.csr_matrix:
    """Assemble global sigma mass matrix; only elements with sigma>0 contribute."""
    Nnode = nodes.shape[0]
    Nelem = triangles.shape[0]

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []

    for ie in range(Nelem):
        sigma = sigma_per_elem[ie]
        if sigma == 0.0:
            continue
        idx = triangles[ie]
        rn = nodes[idx, 0]
        zn = nodes[idx, 1]
        M_e = element_sigma_mass(rn, zn, sigma)
        for j in range(3):
            for k in range(3):
                rows.append(idx[j])
                cols.append(idx[k])
                vals.append(M_e[j, k])

    return sp.csr_matrix((vals, (rows, cols)), shape=(Nnode, Nnode))


def _self_test():
    """Sanity test."""
    # Single triangle, all interior, centroid r ~ 2mm
    rn = np.array([1.0e-3, 3.0e-3, 2.0e-3])
    zn = np.array([0.0, 0.0, 1.5e-3])
    sigma = 5.8e7  # Cu
    M_e = element_sigma_mass(rn, zn, sigma)
    print("element_sigma_mass:")
    print(M_e)
    # Symmetry
    assert np.allclose(M_e, M_e.T), "M_e not symmetric"
    # Positive definite (all eigvals > 0)
    eigs = np.linalg.eigvalsh(M_e)
    print(f"Eigenvalues: {eigs}")
    assert np.all(eigs > -1e-15), "M_e not PSD"
    print("[OK] sigma mass matrix is symmetric and PSD")


if __name__ == "__main__":
    _self_test()
