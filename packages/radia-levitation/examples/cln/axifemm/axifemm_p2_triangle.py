"""axifemm_p2_triangle.py — P2 Henrotte triangle reference (Python).

P2 quadratic triangle in physical (r, z) using the FEMM-Henrotte
convention (matches axifemm_core.py P1 implementation):

The DOF v_i represents A_phi sampled at the i-th node (r_i, z_i).
The vector basis is
    N_A_i(r, z) = (r_i / r) * psi_i(r, z),
where psi_i is the scalar 6-dim Lagrange interpolant in
    span{1, r^2, z, r^4, r^2 z, z^2}
with nodal DOFs at the 3 vertices + 3 edge midpoints (standard P2 layout).
This choice eliminates the A_phi/r^2 vector-Laplacian terms and produces a
clean axisymmetric bilinear form:

  Stiffness  K[i, j] = (2 pi / mu) * r_i * r_j *
              integral over T of
              ( (∂psi_i/∂r)(∂psi_j/∂r) + (∂psi_i/∂z)(∂psi_j/∂z) ) / r dA
  Sigma-mass M[i, j] = 2 pi sigma * r_i * r_j *
              integral over T of psi_i * psi_j / r dA
  Excitation b[i]    = pi B_0 * r_i *
              integral over T of psi_i * r dA           (A_s = B_0 r/2)

For axis-touching cells (some node has r_i = 0), the corresponding row /
column of K and M is automatically zero because of the r_i factor — the
axis DOF is decoupled and the regularity A_phi = 0 at r = 0 is enforced
naturally.

Axis vertices give an integrable 1/r singularity in the integrand; Duffy-
Gauss 8x8 = 64-pt quadrature handles it with finite (though degraded)
accuracy. For the C++ port a specialized axis-touching quadrature should
be added.

History:
  Phase B1a — first attempt using SCALAR Laplacian formulation (r weight,
              no r_i*r_j factor). Disk test matched scalar-Laplacian
              analytical (28.83 us, J_0 zeros) but sphere disagreed with
              Stoll vector-Laplacian by factor 3.4x.
  Phase B1c (this version) — corrected FEMM Henrotte vector formulation.
"""
from __future__ import annotations

import sys

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

MU0 = 4.0e-7 * np.pi
EPS_AXIS = 1.0e-12


# ---------------------------------------------------------------------------
# Duffy-Gauss quadrature on the reference triangle (xi, eta in T_0)
# T_0 = { (xi, eta) : xi >= 0, eta >= 0, xi + eta <= 1 }, area = 1/2.
# Duffy transform: (u, v) in [0, 1]^2 -> (xi, eta) = (u, v(1-u)), Jacobian (1-u).
# n-point tensor-product Gauss-Legendre gives degree (2n-1) in each variable.
# For P2 Henrotte basis, integrand has polynomial degree <= 10, so n=8 (degree
# 15 exact in each variable) is comfortably sufficient.
# ---------------------------------------------------------------------------


def _duffy_gauss_points(n=8):
    """Tensor-product Gauss-Legendre on [0,1]^2 mapped to reference triangle
    via the Duffy transform. Returns array of shape (n^2, 3) of (xi, eta, w).
    Weights sum to area of T_0 = 1/2 (verified in test).
    """
    nodes, weights = np.polynomial.legendre.leggauss(n)
    # Map from [-1, 1] to [0, 1]:
    u_nodes = 0.5 * (nodes + 1.0)
    u_w = 0.5 * weights
    v_nodes = u_nodes.copy()
    v_w = u_w.copy()
    pts = []
    for i, u in enumerate(u_nodes):
        for j, v in enumerate(v_nodes):
            xi = u
            eta = v * (1.0 - u)
            jacobian_duffy = (1.0 - u)  # |d(xi, eta) / d(u, v)|
            w_combined = u_w[i] * v_w[j] * jacobian_duffy
            pts.append((xi, eta, w_combined))
    return np.array(pts)


# Module-level cached quadrature (Duffy-Gauss 8x8 = 64-pt reference triangle).
_QUAD = _duffy_gauss_points(n=8)


# ---------------------------------------------------------------------------
# P2 Henrotte triangle in physical (r, z)
# ---------------------------------------------------------------------------


class AxifemmP2Triangle:
    """P2 quadratic Lagrange triangle in physical (r, z) coords with
    Henrotte basis space span{1, r^2, z, r^4, r^2 z, z^2}.

    Args:
        rn: (6,) array of radial coords of 6 nodes. Order:
                [v0, v1, v2, m01, m12, m20]
            where v_i are vertices and m_ij is midpoint of edge (v_i, v_j).
        zn: (6,) array of axial coords (same ordering).

    The mid-edge nodes are at the *physical* midpoint of the straight edge
    between two vertices (no curved-edge support in Phase B1a).
    """

    NUM_DOFS = 6
    MONOMIALS = 6  # {1, r^2, z, r^4, r^2 z, z^2}

    def __init__(self, rn, zn):
        rn = np.asarray(rn, dtype=float)
        zn = np.asarray(zn, dtype=float)
        if rn.shape != (6,) or zn.shape != (6,):
            raise ValueError("rn, zn must each be length 6")
        self.r = rn.copy()
        self.z = zn.copy()
        self._setup_vandermonde()
        self._setup_physical_geometry()

    def _monomials(self, r, z):
        """Return (6,) vector of monomial values m_j(r, z)."""
        s = r * r
        return np.array([1.0, s, z, s * s, s * z, z * z])

    def _monomial_grads(self, r, z):
        """Return (6, 2) array of gradients [dm_j/dr, dm_j/dz]."""
        s = r * r
        return np.array([
            [0.0,        0.0],
            [2.0 * r,    0.0],
            [0.0,        1.0],
            [4.0 * r * s, 0.0],            # d/dr (r^4) = 4 r^3 = 4 r * s
            [2.0 * r * z, s],              # d/dr (r^2 z) = 2 r z
            [0.0,        2.0 * z],
        ])

    def _setup_vandermonde(self):
        """V[i, j] = m_j(r_i, z_i). Invert -> nodal-basis coeffs.

        Shape functions phi_i(r, z) = sum_j C[i, j] * m_j(r, z)
        with C = V^{-T} (so phi_i(r_k, z_k) = delta_ik).
        Storing C as self.coeffs of shape (6, 6).
        """
        V = np.zeros((6, 6))
        for i in range(6):
            V[i, :] = self._monomials(self.r[i], self.z[i])
        # Solve V^T C^T = I -> C = (V^{-1})^T. Equivalently solve V x = e_i for
        # each i to get coefficients of phi_i in monomial basis.
        try:
            self.coeffs = np.linalg.solve(V, np.eye(6))
            # self.coeffs has shape (6, 6): row j gives monomial coeffs of
            # nodal-basis function for node index given by the column.
            # i.e., phi_i(r, z) = sum_j self.coeffs[j, i] * m_j(r, z).
        except np.linalg.LinAlgError as e:
            raise RuntimeError(
                f"P2 Vandermonde singular for this triangle: rn={self.r} zn={self.z}"
            ) from e

    def _setup_physical_geometry(self):
        """Linear (P1) affine map (xi, eta) -> (r, z) using the 3 vertex
        nodes (indices 0, 1, 2). For straight-edge P2 the mid-edge nodes
        sit exactly on the affine image of the reference-triangle edge
        midpoints, so the Jacobian is constant = 2 * |T|.
        """
        r0, r1, r2 = self.r[0], self.r[1], self.r[2]
        z0, z1, z2 = self.z[0], self.z[1], self.z[2]
        # |T| via the cross-product of edges.
        self.det_J = (r1 - r0) * (z2 - z0) - (z1 - z0) * (r2 - r0)
        self.area = 0.5 * abs(self.det_J)
        self.v0 = np.array([r0, z0])
        self.e1 = np.array([r1 - r0, z1 - z0])
        self.e2 = np.array([r2 - r0, z2 - z0])

    def map_to_physical(self, xi, eta):
        """Reference (xi, eta) -> physical (r, z) using affine map of the
        3 vertex nodes (P1 mapping; straight-edge P2 is consistent)."""
        rp = self.v0[0] + self.e1[0] * xi + self.e2[0] * eta
        zp = self.v0[1] + self.e1[1] * xi + self.e2[1] * eta
        return rp, zp

    def shape(self, r, z):
        """Return (6,) vector phi_i(r, z) for i = 0..5."""
        m = self._monomials(r, z)
        return self.coeffs.T @ m  # phi_i = sum_j C[j, i] m_j

    def shape_grad(self, r, z):
        """Return (6, 2) array d phi_i / d(r, z)."""
        dm = self._monomial_grads(r, z)  # (6, 2)
        return self.coeffs.T @ dm

    # -----------------------------------------------------------------------
    # Element matrices via Hammer 25-point quadrature on the physical triangle
    # -----------------------------------------------------------------------

    def stiffness(self, mu_r=MU0, mu_z=None):
        """FEMM Henrotte stiffness for A_phi:
              K[i, j] = (2 pi r_i r_j) * integral over T of
                        ( d psi_i/dr d psi_j/dr / mu_r
                        + d psi_i/dz d psi_j/dz / mu_z ) / r dA
        Decouples axis DOFs (r_i = 0) automatically.
        """
        if mu_z is None:
            mu_z = mu_r
        K = np.zeros((6, 6))
        jac = abs(self.det_J)
        two_pi = 2.0 * np.pi
        for q in range(_QUAD.shape[0]):
            xi, eta, w = _QUAD[q]
            r, z = self.map_to_physical(xi, eta)
            if r <= EPS_AXIS:
                continue  # integrand has integrable 1/r singularity; skip
            grad = self.shape_grad(r, z)  # (6, 2) = d psi_i / d(r, z)
            d_r = grad[:, 0]
            d_z = grad[:, 1]
            K += w * (np.outer(d_r, d_r) / mu_r
                      + np.outer(d_z, d_z) / mu_z) / r
        K *= jac * two_pi
        K *= np.outer(self.r, self.r)
        return K

    def sigma_mass(self, sigma):
        """FEMM Henrotte sigma-mass for A_phi:
              M[i, j] = (2 pi sigma r_i r_j) * integral over T of
                        psi_i psi_j / r dA
        Decouples axis DOFs (r_i = 0) automatically.
        """
        M = np.zeros((6, 6))
        jac = abs(self.det_J)
        two_pi = 2.0 * np.pi
        for q in range(_QUAD.shape[0]):
            xi, eta, w = _QUAD[q]
            r, z = self.map_to_physical(xi, eta)
            if r <= EPS_AXIS:
                continue
            phi = self.shape(r, z)  # (6,)
            M += w * np.outer(phi, phi) / r
        M *= jac * two_pi * sigma
        M *= np.outer(self.r, self.r)
        return M

    def excitation(self, B0):
        """b[i] = integral over T of A_s . N_A_i dV
                = integral over T of (B_0 r/2) (r_i/r) psi_i 2 pi r dr dz
                = pi B_0 r_i * integral over T of psi_i r dA.
        Axis DOFs (r_i = 0) get zero excitation automatically.
        """
        b = np.zeros(6)
        jac = abs(self.det_J)
        for q in range(_QUAD.shape[0]):
            xi, eta, w = _QUAD[q]
            r, z = self.map_to_physical(xi, eta)
            if r <= EPS_AXIS:
                continue
            phi = self.shape(r, z)
            b += w * phi * r
        b *= jac * np.pi * B0
        b *= self.r
        return b


# ---------------------------------------------------------------------------
# Smoke tests
# ---------------------------------------------------------------------------


def _make_test_triangle(r0, z0, r1, z1, r2, z2):
    """Return AxifemmP2Triangle with straight-edge mid-side nodes."""
    rn = np.array([r0, r1, r2,
                   0.5 * (r0 + r1), 0.5 * (r1 + r2), 0.5 * (r2 + r0)])
    zn = np.array([z0, z1, z2,
                   0.5 * (z0 + z1), 0.5 * (z1 + z2), 0.5 * (z2 + z0)])
    return AxifemmP2Triangle(rn, zn)


def _test_partition_of_unity():
    """sum_i phi_i(r, z) == 1 for any (r, z) in the triangle."""
    tri = _make_test_triangle(1e-3, 0.0, 5e-3, 0.0, 3e-3, 2e-3)
    pts = [(2e-3, 0.5e-3), (3e-3, 1e-3), (4e-3, 1.5e-3)]
    for (r, z) in pts:
        s = tri.shape(r, z).sum()
        assert abs(s - 1.0) < 1e-12, f"Partition-of-unity failed at ({r}, {z}): sum = {s}"
    print("[OK] partition of unity at 3 interior points")


def _test_nodal_interpolation():
    """phi_i(r_j, z_j) = delta_ij."""
    tri = _make_test_triangle(1e-3, 0.0, 5e-3, 0.0, 3e-3, 2e-3)
    for j in range(6):
        phi_at_node = tri.shape(tri.r[j], tri.z[j])
        expected = np.zeros(6)
        expected[j] = 1.0
        err = np.linalg.norm(phi_at_node - expected)
        assert err < 1e-10, f"Nodal interp failed at node {j}: |phi - e_j| = {err}"
    print("[OK] nodal interpolation phi_i(r_j, z_j) = delta_ij")


def _test_axis_decoupling():
    """For a triangle with one axis vertex (r_0 = 0), the corresponding
    row/column of K and M must be zero (axis DOF decoupled in FEMM
    Henrotte convention via the r_i factor)."""
    rn = np.array([0.0, 4e-3, 2e-3,
                   2e-3, 3e-3, 1e-3])
    zn = np.array([0.0, 0.0, 3e-3,
                   0.0, 1.5e-3, 1.5e-3])
    tri = AxifemmP2Triangle(rn, zn)
    K = tri.stiffness(mu_r=MU0)
    M = tri.sigma_mass(sigma=5.8e7)
    # Axis vertex is node index 0 (r_0 = 0). Its entire row and column
    # should be zero.
    K_row0 = np.max(np.abs(K[0, :]))
    M_row0 = np.max(np.abs(M[0, :]))
    assert K_row0 < 1e-30, f"Axis K row not zero: max={K_row0}"
    assert M_row0 < 1e-30, f"Axis M row not zero: max={M_row0}"
    print(f"[OK] axis-vertex DOF decoupled (K row, M row exactly zero)")


def _test_mass_positive_definite():
    """sigma-mass matrix should be SPD (positive definite)."""
    tri = _make_test_triangle(1e-3, 0.0, 5e-3, 0.0, 3e-3, 2e-3)
    M = tri.sigma_mass(sigma=5.8e7)
    eigvals = np.linalg.eigvalsh(M)
    assert np.all(eigvals > 0), f"M not SPD: min eigenvalue = {eigvals.min()}"
    print(f"[OK] sigma-mass SPD, eigvals range [{eigvals.min():.2e}, "
          f"{eigvals.max():.2e}]")


def _test_symmetry():
    """K and M are symmetric."""
    tri = _make_test_triangle(1e-3, 0.0, 5e-3, 0.0, 3e-3, 2e-3)
    K = tri.stiffness()
    M = tri.sigma_mass(sigma=5.8e7)
    err_K = np.linalg.norm(K - K.T) / max(np.linalg.norm(K), 1e-30)
    err_M = np.linalg.norm(M - M.T) / max(np.linalg.norm(M), 1e-30)
    assert err_K < 1e-10, f"K asymmetry: {err_K}"
    assert err_M < 1e-10, f"M asymmetry: {err_M}"
    print(f"[OK] K, M symmetric to {max(err_K, err_M):.2e}")


def _test_area_via_quadrature():
    """Integrate "1" over the triangle: should equal |T|."""
    rn = np.array([1e-3, 5e-3, 3e-3, 3e-3, 4e-3, 2e-3])
    zn = np.array([0.0, 0.0, 2e-3, 0.0, 1e-3, 1e-3])
    tri = AxifemmP2Triangle(rn, zn)
    total = 0.0
    jac = abs(tri.det_J)
    for q in range(_QUAD.shape[0]):
        xi, eta, w = _QUAD[q]
        total += w
    total *= jac
    expected = tri.area
    err = abs(total - expected) / expected
    assert err < 1e-12, f"Triangle area quadrature: {total} vs {expected}, rel {err}"
    print(f"[OK] triangle area via Hammer 25 = {total:.6e} (expected {expected:.6e})")


def main():
    print("=" * 60)
    print("axifemm P2 Henrotte triangle Python reference — smoke tests")
    print("=" * 60)
    _test_area_via_quadrature()
    _test_partition_of_unity()
    _test_nodal_interpolation()
    _test_symmetry()
    _test_axis_decoupling()
    _test_mass_positive_definite()
    print()
    print("All smoke tests passed.")


if __name__ == "__main__":
    main()
