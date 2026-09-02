"""09_wedge_basis.py

Non-separable wedge basis for cube Mixed Galerkin -- attempt to break
the 0.33% wall-band floor.

The separable f(x)f(y)f(z) basis is asymptotically rank-deficient
(diagnosed in 08_edge_corner_basis.py: condition number explodes,
solver gives garbage).  Here we try a basis that has intrinsic 2D
structure on the (x, y) plane normal to a z-axis edge:

    psi_e(r, s) = [e^{-t r_xy} - e^{-tx} - e^{-ty} + 1]
                  * exp(-c (x^2 + y^2) / L^2)
                  * f(z)

with r_xy = sqrt(x^2 + y^2), c chosen so the bump kills the far-face
contribution (c ~ 10 / (tL/2)^2 ensures decay before x = L/2).  The
bracket vanishes on x = 0 and y = 0 (adjacent face Dirichlet); f(z)
vanishes on z = 0 and z = L; the bump handles x = L, y = L far-decay.

We add ONE z-edge basis at corner (x=0, y=0) and 11 symmetry images
to cover all 12 cube edges.  Numerical Gauss-Legendre quadrature for
all integrals.

The fundamental question: does Schur composition with non-separable
basis break the 0.33% floor?  Schur is just a generic projection;
the answer is determined by whether the basis is asymptotically
linearly INDEPENDENT from the corner envelope at high f.
"""

from __future__ import annotations

import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.cube3d_foster import K_SIBC_cube3d, Y_DC_cube3d, Y_mellin_cube3d

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cube3d(L, SIGMA)
K_SIBC = K_SIBC_cube3d(L, SIGMA, MU)
V_CUBE = L ** 3

# Gauss-Legendre quadrature with skin-graded nodes near the boundary.
# Use refinement: more points near x = 0, L (where wedge basis lives).
N_QUAD_PER_AXIS = 24


def _glnodes(n):
    """Gauss-Legendre nodes and weights on [-1, 1]."""
    return np.polynomial.legendre.leggauss(n)


def _gradient_xy(x, y, t):
    """Compute |grad psi_xy|^2 contribution at one point.

    psi_xy = e^{-t r} - e^{-tx} - e^{-ty} + 1, r = sqrt(x^2+y^2)
    d/dx (psi_xy) = -t (x/r) e^{-t r} + t e^{-tx}
    d/dy (psi_xy) = -t (y/r) e^{-t r} + t e^{-ty}

    Returns the (du/dx)^2 + (du/dy)^2 for u = e^{-tr} - e^{-tx} - e^{-ty} + 1.
    """
    if abs(x) + abs(y) < 1e-15:
        return complex(0.0)  # singular but supports measure zero
    r = cmath.sqrt(x * x + y * y)
    dpsi_dx = -t * (x / r) * cmath.exp(-t * r) + t * cmath.exp(-t * x)
    dpsi_dy = -t * (y / r) * cmath.exp(-t * r) + t * cmath.exp(-t * y)
    return dpsi_dx * dpsi_dx + dpsi_dy * dpsi_dy


def _grad_xy_field_dot_xyz(x, y, t):
    """(d/dx, d/dy, 0) of u = e^{-tr} - e^{-tx} - e^{-ty} + 1, returned as tuple."""
    r = cmath.sqrt(x * x + y * y) if (x * x + y * y) > 0 else complex(1e-15)
    dpsi_dx = -t * (x / r) * cmath.exp(-t * r) + t * cmath.exp(-t * x)
    dpsi_dy = -t * (y / r) * cmath.exp(-t * r) + t * cmath.exp(-t * y)
    return dpsi_dx, dpsi_dy


def _u_xy(x, y, t):
    """u(x, y) = e^{-t r_xy} - e^{-tx} - e^{-ty} + 1."""
    r = cmath.sqrt(x * x + y * y)
    return cmath.exp(-t * r) - cmath.exp(-t * x) - cmath.exp(-t * y) + 1.0


def _bump(x, y, c_param):
    """Smooth far-decay bump: exp(-c (x^2 + y^2) / L^2)."""
    return math.exp(-c_param * (x * x + y * y) / (L * L))


def _bump_grad(x, y, c_param):
    """Gradient of bump in (x, y) plane."""
    g = _bump(x, y, c_param)
    return -2 * c_param * x / (L * L) * g, -2 * c_param * y / (L * L) * g


# ---------------------------------------------------------------------------
# Numerical integration of K_ss, K_bs, b for wedge basis
# ---------------------------------------------------------------------------


def assemble_wedge_z(s, c_param=10.0):
    """Build the contribution of 1 z-axis edge (at x=0, y=0) wedge basis to
    the surface block K_ss, K_bs (vs bulk phi_111), and b vector.

    Wedge basis: psi_e(r, s) = u(x, y) * bump(x, y) * f(z)
    where u(x, y) = e^{-t r_xy} - e^{-tx} - e^{-ty} + 1
          f(z)    = cosh((z-L/2) t) / cosh(L t / 2) - 1
    """
    t = cmath.sqrt(s * MS)
    nodes, weights = _glnodes(N_QUAD_PER_AXIS)
    # map nodes from [-1, 1] to [0, L]
    pts_x = (nodes + 1) * L / 2
    wts_x = weights * L / 2

    # 1D integrals over z (analytic for f(z))
    Lt2 = L * t / 2
    if abs(Lt2.real) > 50:
        cosh_Lt2 = cmath.exp(abs(Lt2)) / 2
        tanh_Lt2 = complex(1.0, 0.0) if Lt2.real > 0 else complex(-1.0, 0.0)
        sech2_Lt2 = complex(0.0, 0.0)
    else:
        cosh_Lt2 = cmath.cosh(Lt2)
        tanh_Lt2 = cmath.tanh(Lt2)
        sech2_Lt2 = 1.0 / cosh_Lt2 ** 2

    # F0 = ∫_0^L f dz
    F0_z = (2 / t) * tanh_Lt2 - L
    # J3 = ∫_0^L f^2 dz
    J3_z = (L / 2) * sech2_Lt2 - 3 * tanh_Lt2 / t + L
    # J2 = ∫_0^L (df/dz)^2 dz
    J2_z = t * tanh_Lt2 - (t ** 2 * L / 2) * sech2_Lt2
    # phi1_z = ∫_0^L f(z) sin(pi z / L) dz
    phi1_z = 2 * L * math.pi / (math.pi ** 2 + L ** 2 * t ** 2) - 2 * L / math.pi
    # int sin(pi z / L) dz = 2 L / pi
    int_sin_z = 2 * L / math.pi
    # int sin^2(pi z / L) dz = L / 2
    int_sin2_z = L / 2

    # Numerical 2D quadrature over (x, y)
    # I_mass_xy = ∫∫ |u * bump|^2 dx dy
    # I_stiff_xy = ∫∫ |grad_xy (u * bump)|^2 dx dy
    # I_b_xy    = ∫∫ u * bump dx dy
    # I_overlap_xy_phi = ∫∫ u * bump * sin(pi x/L) sin(pi y/L) dx dy

    I_mass = complex(0.0)
    I_stiff = complex(0.0)
    I_b = complex(0.0)
    I_overlap = complex(0.0)
    for i, x in enumerate(pts_x):
        for j, y in enumerate(pts_x):
            w = wts_x[i] * wts_x[j]
            u = _u_xy(x, y, t)
            b = _bump(x, y, c_param)
            ub = u * b
            I_mass += w * ub * ub  # bilinear, no conjugation
            # Gradient of u*b = (du/dx * b + u * db/dx, du/dy * b + u * db/dy)
            du_dx, du_dy = _grad_xy_field_dot_xyz(x, y, t)
            db_dx, db_dy = _bump_grad(x, y, c_param)
            grad_x = du_dx * b + u * db_dx
            grad_y = du_dy * b + u * db_dy
            I_stiff += w * (grad_x ** 2 + grad_y ** 2)
            I_b += w * ub
            I_overlap += w * ub * math.sin(math.pi * x / L) * math.sin(math.pi * y / L)

    # Combine 2D xy integrals with 1D z integrals
    # K_ss (mass+stiff for psi_e):
    #   ∫(grad psi_e)^2 dV = I_stiff_xy * J3_z (grad in xy) + I_mass_xy * J2_z (grad in z)
    #   ∫|psi_e|^2 dV = I_mass_xy * J3_z
    K_ee_grad = I_stiff * J3_z + I_mass * J2_z
    K_ee_mass = I_mass * J3_z
    K_ee = K_ee_grad + s * MS * K_ee_mass

    # K_bs (overlap of psi_e with bulk phi_{1,1,1} = sin(pi x/L) sin(pi y/L) sin(pi z/L)):
    #   ∫(grad psi_e . grad phi + sMS psi_e * phi) dV
    #   = (lambda_111 + sMS) * ∫ psi_e * phi_111 dV   (since phi is eigenfunction of Lap)
    #   ∫ psi_e * phi_111 dV = I_overlap_xy_phi * phi1_z
    lambda_111 = 3 * math.pi ** 2 / L ** 2
    K_be = (lambda_111 + s * MS) * I_overlap * phi1_z

    # b_e = ∫ psi_e dV = I_b_xy * F0_z
    b_e = I_b * F0_z

    return K_ee, K_be, b_e


# ---------------------------------------------------------------------------
# Existing baseline (rank-1 bulk + 1-DOF corner envelope)
# ---------------------------------------------------------------------------


def _tanh_safe(z):
    if abs(z.real) > 50:
        return complex(1.0, 0.0) if z.real > 0 else complex(-1.0, 0.0)
    return cmath.tanh(z)


def _sech2_safe(z):
    if abs(z.real) > 50:
        return complex(0.0, 0.0)
    return 1.0 / cmath.cosh(z) ** 2


def J_2_exact(t):
    Lt2 = L * t / 2
    return t * _tanh_safe(Lt2) - (t ** 2 * L / 2) * _sech2_safe(Lt2)


def J_3_exact(t):
    Lt2 = L * t / 2
    return (L / 2) * _sech2_safe(Lt2) - 3 * _tanh_safe(Lt2) / t + L


def F_0_exact(t):
    return (2 / t) * _tanh_safe(L * t / 2) - L


def phi1_overlap(t):
    return 2 * L * math.pi / (math.pi ** 2 + L ** 2 * t ** 2) - 2 * L / math.pi


def Y_mixed_baseline(s):
    """rank-1 bulk + 1-DOF psi_xyz tensor envelope (baseline)."""
    t = cmath.sqrt(s * MS)
    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)
    phi1 = phi1_overlap(t)
    lambda_111 = 3 * math.pi ** 2 / L ** 2
    K_b = (lambda_111 + s * MS) * (L / 2) ** 3
    b_b = (2 * L / math.pi) ** 3
    K_ss = 3 * J_2 * J_3 ** 2 + s * MS * J_3 ** 3
    b_psi = F_0 ** 3
    K_bs = (lambda_111 + s * MS) * phi1 ** 3
    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


def Y_mixed_with_wedge(s, c_param=10.0):
    """rank-1 bulk + 1-DOF corner envelope + 1 z-edge wedge (at x=0, y=0)."""
    t = cmath.sqrt(s * MS)
    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)
    phi1 = phi1_overlap(t)
    lambda_111 = 3 * math.pi ** 2 / L ** 2

    # Bulk
    K_b = (lambda_111 + s * MS) * (L / 2) ** 3
    b_b = (2 * L / math.pi) ** 3

    # Corner envelope
    K_ss_corn = 3 * J_2 * J_3 ** 2 + s * MS * J_3 ** 3
    b_psi_corn = F_0 ** 3
    K_bs_corn = (lambda_111 + s * MS) * phi1 ** 3

    # Wedge (1 edge at z-axis, corner (0, 0))
    K_ee, K_be, b_e = assemble_wedge_z(s, c_param=c_param)

    # Cross between corner envelope and wedge basis:
    #   psi_corn = f(x) f(y) f(z), psi_wedge = u(x,y) bump f(z)
    #   ∫(grad psi_corn . grad psi_wedge + sMS psi_corn psi_wedge) dV
    # This requires another 2D xy quadrature.  Skip for now and treat the
    # two surface bases as independent (off-diagonal K_corn_wedge = 0).
    # This is approximate but tests the principle.
    K_corn_wedge = 0.0

    # Assemble 3x3 (bulk + corner + wedge)
    K_full = np.array(
        [
            [K_b,         K_bs_corn,    K_be],
            [K_bs_corn,   K_ss_corn,    K_corn_wedge],
            [K_be,        K_corn_wedge, K_ee],
        ],
        dtype=complex,
    )
    b_full = np.array([b_b, b_psi_corn, b_e], dtype=complex)
    xi = np.linalg.solve(K_full, -s * MS * b_full)
    v_avg = (xi @ b_full) / V_CUBE
    return Y_DC * (1 + v_avg)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------


def main():
    print("=== cube Mixed Galerkin: non-separable wedge basis test ===")
    print(f"L = {L * 1e3} mm, K_SIBC = {K_SIBC:.4e}, Quad pts per axis = {N_QUAD_PER_AXIS}")
    print("Wedge: psi_e(r) = [e^{-tr} - e^{-tx} - e^{-ty} + 1] * exp(-10 r^2/L^2) * f(z)")
    print("Symmetry would give 12 edges; this is 1 z-axis edge test.")
    print()
    print(
        f"{'f (Hz)':>10}  {'|Y_Mellin|':>14}  "
        f"{'err_baseline':>12}  {'err_wedge':>12}  {'reduction':>10}"
    )

    for f in [1e3, 1e4, 1e5, 1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        Y_M = Y_mellin_cube3d(s, L, SIGMA, MU)
        Y_b = Y_mixed_baseline(s)
        try:
            Y_w = Y_mixed_with_wedge(s)
            err_b = abs(Y_b - Y_M) / abs(Y_M) * 100
            err_w = abs(Y_w - Y_M) / abs(Y_M) * 100
            red = err_b / err_w if err_w > 1e-12 else float("inf")
            print(
                f"{f:10.2e}  {abs(Y_M):14.6e}  "
                f"{err_b:11.6f}%  {err_w:11.6f}%  {red:9.2f}x"
            )
        except Exception as e:
            print(f"{f:10.2e}  failed: {e}")


if __name__ == "__main__":
    main()
