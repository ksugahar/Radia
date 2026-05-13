"""axifemm_quad_q2_curved.py — Phase B5 curved Q2 Henrotte element (Python prototype).

9-node general quadrilateral in physical (r, z) for the axisymmetric FEMM-
Henrotte formulation, with biquadratic isoparametric geometric mapping.

DOF basis (matching axifemm_quad_q2.py axis-aligned case):
    monomials in (s = r^2, z):
        {1, s, s^2, z, s z, s^2 z, z^2, s z^2, s^2 z^2}
    Vandermonde over physical node positions (r_i, z_i).
    psi_i(r, z) = sum_j C[j, i] m_j(r, z),  with C = V^{-1},  psi_i(r_k, z_k) = delta_ik.

Even r-powers ensure the Henrotte form ((d psi/dr)^2 + (d psi/dz)^2) / r
integrand is finite at the axis (r → 0) without Duffy clustering.

Geometric mapping (xi, eta) ∈ [-1, 1]^2 → (r, z):
    standard biquadratic Lagrange interpolation through 9 nodes:
        r(xi, eta) = sum_i N_i(xi, eta) r_i,  z analogous.
    where N_i are biquadratic Lagrange on the reference square.

Node ordering (matches axifemm_quad_q2.py and standard FEMM Q2 layout):
    0: ( -1, -1 ) corner
    1: ( +1, -1 ) corner
    2: ( +1, +1 ) corner
    3: ( -1, +1 ) corner
    4: (  0, -1 ) edge bottom
    5: ( +1,  0 ) edge right
    6: (  0, +1 ) edge top
    7: ( -1,  0 ) edge left
    8: (  0,  0 ) face center

FEMM-Henrotte bilinear forms (V_i = A_phi at node i):
    K[i, j] = (2 pi r_i r_j) ∫_T ((∂psi_i/∂r)(∂psi_j/∂r)/mu_r
                                 + (∂psi_i/∂z)(∂psi_j/∂z)/mu_z) / r dr dz
    M[i, j] = (2 pi sigma r_i r_j) ∫_T psi_i psi_j / r dr dz
    b[i]    = pi B_0 r_i ∫_T psi_i r dr dz       (uniform A_s = B_0 r / 2)

Axis-touching elements (some r_i = 0) decouple automatically via the r_i factor.

Phase B5 status: Python prototype + element-level tests. C++ port (NGSolve
AxiHenrotteQ2FESpace) deferred to Phase B6.
"""
from __future__ import annotations

import sys

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

MU0 = 4.0e-7 * np.pi
EPS_AXIS = 1.0e-12


# ---------------------------------------------------------------------------
# Reference-element data (FEMM Q2 node ordering)
# ---------------------------------------------------------------------------

# (xi_i, eta_i) reference coords for the 9 nodes.
_XI_REF = np.array([-1.0,  1.0,  1.0, -1.0,  0.0,  1.0,  0.0, -1.0,  0.0])
_ETA_REF = np.array([-1.0, -1.0,  1.0,  1.0, -1.0,  0.0,  1.0,  0.0,  0.0])

# 1D Lagrange index in {0, 1, 2} for nodes on {-1, 0, +1}:
#   index 0 ↔ t = -1,  index 1 ↔ t = 0,  index 2 ↔ t = +1.
_LIDX_XI = np.array([0, 2, 2, 0, 1, 2, 1, 0, 1], dtype=int)
_LIDX_ETA = np.array([0, 0, 2, 2, 0, 1, 2, 1, 1], dtype=int)


def _lagrange_1d(t):
    """Return (3,) of 1D Lagrange basis values at t for nodes {-1, 0, +1}."""
    return np.array([0.5 * t * (t - 1.0), 1.0 - t * t, 0.5 * t * (t + 1.0)])


def _lagrange_1d_deriv(t):
    """Return (3,) of 1D Lagrange basis derivatives at t for nodes {-1, 0, +1}."""
    return np.array([t - 0.5, -2.0 * t, t + 0.5])


# ---------------------------------------------------------------------------
# AxifemmQuadQ2Curved
# ---------------------------------------------------------------------------


class AxifemmQuadQ2Curved:
    """9-node curved biquadratic FEMM-Henrotte axisymmetric element.

    Args:
        rn: (9,) radial coords of 9 nodes (FEMM ordering above).
        zn: (9,) axial coords of 9 nodes.
    """

    NUM_DOFS = 9
    MONOMIALS = 9

    def __init__(self, rn, zn):
        rn = np.asarray(rn, dtype=float)
        zn = np.asarray(zn, dtype=float)
        if rn.shape != (9,) or zn.shape != (9,):
            raise ValueError("rn, zn must each be length 9")
        if np.any(rn < -EPS_AXIS):
            raise ValueError(f"r-coords must be non-negative: {rn}")
        self.r = rn.copy()
        self.z = zn.copy()
        # Conditioning: scale r by L_r = max |r_i|, scale z by L_z = z half-range.
        # Scaling r preserves even-r-power structure (s' = (r/L_r)^2 = s/L_r^2).
        # Shifting z by element midpoint puts z' ∈ [-1, +1] for box element.
        self.L_r = max(float(np.max(np.abs(self.r))), 1e-30)
        z_min, z_max = float(np.min(self.z)), float(np.max(self.z))
        self.z_c = 0.5 * (z_min + z_max)
        self.L_z = max(0.5 * (z_max - z_min), 1e-30)
        self._setup_vandermonde()

    # ---- monomial basis in scaled (s' = (r/L_r)^2, z' = (z-z_c)/L_z) ----

    def _scaled(self, r, z):
        return r / self.L_r, (z - self.z_c) / self.L_z

    def _monomials(self, r, z):
        """Return (9,) monomial values m_j in scaled coords (r', z')."""
        rp, zp = self._scaled(r, z)
        sp = rp * rp
        return np.array([
            1.0,                # 0
            sp,                 # 1: s'
            sp * sp,            # 2: s'^2
            zp,                 # 3: z'
            sp * zp,            # 4: s' z'
            sp * sp * zp,       # 5: s'^2 z'
            zp * zp,            # 6: z'^2
            sp * zp * zp,       # 7: s' z'^2
            sp * sp * zp * zp,  # 8: s'^2 z'^2
        ])

    def _monomial_grads_rz(self, r, z):
        """Return (9, 2): [d m_j / dr, d m_j / dz] in *physical* coords.
        Uses chain rule via scaled coords (r' = r/L_r, z' = (z-z_c)/L_z):
            d m_j / dr = (∂m_j/∂r') / L_r,
            d m_j / dz = (∂m_j/∂z') / L_z.
        """
        rp, zp = self._scaled(r, z)
        sp = rp * rp
        dr_prime = np.array([
            0.0,
            2.0 * rp,
            4.0 * rp * sp,
            0.0,
            2.0 * rp * zp,
            4.0 * rp * sp * zp,
            0.0,
            2.0 * rp * zp * zp,
            4.0 * rp * sp * zp * zp,
        ])
        dz_prime = np.array([
            0.0,
            0.0,
            0.0,
            1.0,
            sp,
            sp * sp,
            2.0 * zp,
            2.0 * sp * zp,
            2.0 * sp * sp * zp,
        ])
        return np.stack([dr_prime / self.L_r, dz_prime / self.L_z], axis=1)

    def _setup_vandermonde(self):
        """Compute monomial-to-nodal coeffs C = V^{-1}.

        psi_i(r, z) = sum_j C[j, i] m_j(r, z),  psi_i(r_k, z_k) = delta_ik.
        """
        V = np.zeros((9, 9))
        for i in range(9):
            V[i, :] = self._monomials(self.r[i], self.z[i])
        try:
            self.coeffs = np.linalg.solve(V, np.eye(9))
        except np.linalg.LinAlgError as e:
            raise RuntimeError(
                f"Q2 curved Vandermonde singular for this element:\n"
                f"  rn = {self.r}\n  zn = {self.z}"
            ) from e
        self._vandermonde_cond = np.linalg.cond(V)

    # ---- nodal basis evaluation (in physical (r, z)) ----

    def shape(self, r, z):
        """Return (9,) psi_i(r, z) for i = 0..8."""
        return self.coeffs.T @ self._monomials(r, z)

    def shape_grad_rz(self, r, z):
        """Return (9, 2) [d psi_i / dr, d psi_i / dz] for i = 0..8."""
        return self.coeffs.T @ self._monomial_grads_rz(r, z)

    # ---- biquadratic isoparametric geometric mapping ----

    @staticmethod
    def _geom_shape_ref(xi, eta):
        """Return (9,) N_i(xi, eta), the biquadratic Lagrange shape functions
        on the reference square [-1, 1]^2 with the FEMM Q2 node ordering."""
        Lx = _lagrange_1d(xi)
        Ly = _lagrange_1d(eta)
        return Lx[_LIDX_XI] * Ly[_LIDX_ETA]

    @staticmethod
    def _geom_shape_grad_ref(xi, eta):
        """Return (9, 2) [dN_i/dxi, dN_i/deta]."""
        Lx = _lagrange_1d(xi)
        dLx = _lagrange_1d_deriv(xi)
        Ly = _lagrange_1d(eta)
        dLy = _lagrange_1d_deriv(eta)
        dN = np.zeros((9, 2))
        dN[:, 0] = dLx[_LIDX_XI] * Ly[_LIDX_ETA]
        dN[:, 1] = Lx[_LIDX_XI] * dLy[_LIDX_ETA]
        return dN

    def map_to_physical(self, xi, eta):
        N = self._geom_shape_ref(xi, eta)
        return float(self.r @ N), float(self.z @ N)

    def jacobian(self, xi, eta):
        """Return 2x2 Jacobian [[dr/dxi, dr/deta], [dz/dxi, dz/deta]]."""
        dN = self._geom_shape_grad_ref(xi, eta)
        return np.array([
            [self.r @ dN[:, 0], self.r @ dN[:, 1]],
            [self.z @ dN[:, 0], self.z @ dN[:, 1]],
        ])

    # ---- element matrices ----

    def stiffness(self, mu_r=MU0, mu_z=None, n_quad=10):
        """FEMM Henrotte stiffness:
              K[i, j] = (2 pi r_i r_j) ∫_T
                  ((∂psi_i/∂r)(∂psi_j/∂r)/mu_r
                 + (∂psi_i/∂z)(∂psi_j/∂z)/mu_z) / r dr dz
        Tensor-product Gauss-Legendre n_quad x n_quad on [-1, 1]^2.
        """
        if mu_z is None:
            mu_z = mu_r
        pts, wts = np.polynomial.legendre.leggauss(n_quad)
        K = np.zeros((9, 9))
        for ix, xi in enumerate(pts):
            for iy, eta in enumerate(pts):
                r, z = self.map_to_physical(xi, eta)
                if r <= EPS_AXIS:
                    continue  # axis-edge sample point: r_i factor will zero out anyway
                J = self.jacobian(xi, eta)
                detJ = abs(J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0])
                grad = self.shape_grad_rz(r, z)
                d_r = grad[:, 0]
                d_z = grad[:, 1]
                w = wts[ix] * wts[iy] * detJ
                K += w * (np.outer(d_r, d_r) / mu_r
                          + np.outer(d_z, d_z) / mu_z) / r
        K *= 2.0 * np.pi
        K *= np.outer(self.r, self.r)
        return K

    def sigma_mass(self, sigma, n_quad=10):
        """FEMM Henrotte sigma-mass:
              M[i, j] = (2 pi sigma r_i r_j) ∫_T psi_i psi_j / r dr dz
        """
        pts, wts = np.polynomial.legendre.leggauss(n_quad)
        M = np.zeros((9, 9))
        for ix, xi in enumerate(pts):
            for iy, eta in enumerate(pts):
                r, z = self.map_to_physical(xi, eta)
                if r <= EPS_AXIS:
                    continue
                J = self.jacobian(xi, eta)
                detJ = abs(J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0])
                phi = self.shape(r, z)
                w = wts[ix] * wts[iy] * detJ
                M += w * np.outer(phi, phi) / r
        M *= 2.0 * np.pi * sigma
        M *= np.outer(self.r, self.r)
        return M

    def excitation(self, B0, n_quad=10):
        """Uniform B_z = B_0 source (A_s = B_0 r / 2):
              b[i] = pi B_0 r_i ∫_T psi_i r dr dz
        """
        pts, wts = np.polynomial.legendre.leggauss(n_quad)
        b = np.zeros(9)
        for ix, xi in enumerate(pts):
            for iy, eta in enumerate(pts):
                r, z = self.map_to_physical(xi, eta)
                J = self.jacobian(xi, eta)
                detJ = abs(J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0])
                phi = self.shape(r, z)
                w = wts[ix] * wts[iy] * detJ
                b += w * phi * r
        b *= np.pi * B0
        b *= self.r
        return b


# ---------------------------------------------------------------------------
# Helper constructors
# ---------------------------------------------------------------------------


def make_axis_aligned_q2_curved(ra, rb, za, zb, mid_r="quad"):
    """Build an axis-aligned curved Q2 element.

    mid_r ∈ {"quad", "arith"} chooses the placement of the edge mid-nodes
    along r:
        "quad":  rm = sqrt((ra^2 + rb^2) / 2)   (matches axifemm_quad_q2.py)
        "arith": rm = (ra + rb) / 2             (natural for unmapped (r, z) Q2)

    Mid-nodes along z always use the arithmetic mean (z is the natural coord).
    """
    if mid_r == "quad":
        rm = float(np.sqrt(0.5 * (ra * ra + rb * rb)))
    elif mid_r == "arith":
        rm = 0.5 * (ra + rb)
    else:
        raise ValueError(f"mid_r must be 'quad' or 'arith', got {mid_r}")
    zm = 0.5 * (za + zb)
    rn = np.array([ra, rb, rb, ra, rm, rb, rm, ra, rm])
    zn = np.array([za, za, zb, zb, za, zm, zb, zm, zm])
    return AxifemmQuadQ2Curved(rn, zn)


def make_straight_q2_curved(rs, zs):
    """Build a (general 4-sided) Q2 element with straight edges + midpoint mids.

    Args:
        rs, zs: (4,) arrays of corner (r, z) coords in CCW order.
    """
    rs = np.asarray(rs, dtype=float)
    zs = np.asarray(zs, dtype=float)
    if rs.shape != (4,) or zs.shape != (4,):
        raise ValueError("rs, zs must each be length 4 (4 corners)")
    rn = np.zeros(9)
    zn = np.zeros(9)
    rn[0:4] = rs
    zn[0:4] = zs
    # edge mids: 4 = mid(0,1), 5 = mid(1,2), 6 = mid(2,3), 7 = mid(3,0)
    rn[4] = 0.5 * (rs[0] + rs[1]);  zn[4] = 0.5 * (zs[0] + zs[1])
    rn[5] = 0.5 * (rs[1] + rs[2]);  zn[5] = 0.5 * (zs[1] + zs[2])
    rn[6] = 0.5 * (rs[2] + rs[3]);  zn[6] = 0.5 * (zs[2] + zs[3])
    rn[7] = 0.5 * (rs[3] + rs[0]);  zn[7] = 0.5 * (zs[3] + zs[0])
    # center: average of 4 corners
    rn[8] = 0.25 * rs.sum();  zn[8] = 0.25 * zs.sum()
    return AxifemmQuadQ2Curved(rn, zn)
