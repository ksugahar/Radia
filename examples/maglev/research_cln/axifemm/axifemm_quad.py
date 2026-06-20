"""
axifemm_quad.py — Phase 1b: FEMM/Henrotte Q1 axisymmetric on axis-aligned rectangles.

Basis: phi(r, z) = a + b*r^2 + c*z + d*r^2*z   (Q1 in (s=r^2, z) coords)

Field components:
  B_z = (1/(2 pi r)) d phi/dr = (b + d*z)/pi               (linear in z, const in r)
  B_r = -(1/(2 pi r)) d phi/dz = -(c + d*s)/(2 pi sqrt(s))  (1/r structure, no singularity)

For an axis-aligned rectangle [r_a, r_b] x [z_a, z_b] with nodes:
  node 0: (r_a, z_a)        (sa, za)
  node 1: (r_b, z_a)        (sb, za)
  node 2: (r_b, z_b)        (sb, zb)
  node 3: (r_a, z_b)        (sa, zb)

We solve for V_j = A_phi at node j (same convention as Triangle P1 axifemm).

For axis-touching rectangle (r_a = 0 => sa = 0): the 1/s singularity in
B_r^2 integral produces log(sb/sa) -> infinity. Need limit case (Limit[..., sa -> 0]).

This file provides:
  - element_matrices_quad(ra, rb, za, zb, mu_r, mu_z): (M_e, Mr, Mz)
  - element_sigma_mass_quad(ra, rb, za, zb, sigma): M_sigma
  - assemble_quad(...): global K and M_sigma
  - structured_mesh_disk(R_disk, t_disk, R_air, Z_air, NR_disk, Nz_disk, NR_air, Nz_air)

Closed-form element matrices: filled in from Mathematica derivation
(quad_q1_henrotte_matrices.json). Until then, numerical quadrature is used.
"""

from __future__ import annotations

import os
import json
import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))


# 4x4 Gauss-Legendre quadrature on reference square (-1, 1) x (-1, 1).
_GL4_x = np.array([-0.8611363115940526, -0.3399810435848563,
                    0.3399810435848563,  0.8611363115940526])
_GL4_w = np.array([ 0.3478548451374538,  0.6521451548625461,
                    0.6521451548625461,  0.3478548451374538])


def _phi_basis_inverse(ra, rb, za, zb):
    """Compute the (a, b, c, d) -> (phi_0, phi_1, phi_2, phi_3) inverse.

    phi(r, z) = a + b*r^2 + c*z + d*r^2*z

    Nodal values:
        phi_j = a + b*s_j + c*z_j + d*s_j*z_j
    where (s_j, z_j) are the (r^2, z) coords of node j.

    Returns the 4x4 matrix M such that (a, b, c, d) = M . (phi_0, phi_1, phi_2, phi_3).
    """
    sa, sb = ra ** 2, rb ** 2
    # node order: (sa, za), (sb, za), (sb, zb), (sa, zb)
    V = np.array([
        [1, sa, za, sa * za],
        [1, sb, za, sb * za],
        [1, sb, zb, sb * zb],
        [1, sa, zb, sa * zb],
    ])
    return np.linalg.inv(V)


def element_matrices_quad(ra, rb, za, zb, mu_r, mu_z):
    """
    Return (M_e, Mr, Mz) for an axis-aligned rectangle in V (= A_phi) coordinates.

    Energy: W = (1/2) Integrate( B_r^2/mu_r + B_z^2/mu_z, 2 pi r dr dz )
    With phi = a + b*s + c*z + d*s*z, s = r^2:
      B_z = (b + d*z)/pi
      B_r = -(c + d*s)/(2 pi sqrt(s))

    W_z = Integrate( (b + d z)^2 / (2 pi mu_z), {s, sa, sb}, {z, za, zb} )
        = (1/(2 pi mu_z)) (sb - sa) ((b za)^2 + (b zb)^2 ... ) -- closed form
    W_r involves Integrate( 1/s, sa, sb ) = log(sb/sa).

    For axis-touching (sa = 0): W_r has log singularity -> use Limit treatment.

    Implementation: closed-form when available (Mathematica-derived), otherwise
    numerical Gauss quadrature.
    """
    return _element_matrices_quad_closed_form(ra, rb, za, zb, mu_r, mu_z)


def _element_matrices_quad_closed_form(ra, rb, za, zb, mu_r, mu_z):
    """
    Closed-form 4x4 element matrices in V coords.

    Derived in derive_quad_q1_henrotte.wls:
      W_z = (sb - sa)/(6 pi mu_z) * [b^2 (za^2 + za zb + zb^2) + 3 b d ... + d^2 ...]
      W_r = (zb - za)/(8 pi mu_r) * [c^2 log(sb/sa) + 2 c d (sb - sa) + d^2 (sb^2 - sa^2)/2]

    Will be replaced once Mathematica derivation completes; for now use numerical.
    """
    # Coefficient transformation matrix: (a,b,c,d) = inv_V . (phi_0, phi_1, phi_2, phi_3)
    inv_V = _phi_basis_inverse(ra, rb, za, zb)

    # Closed-form for W_z and W_r in (b, c, d) coefficients.
    # (W_z does not depend on a; W_r does not depend on a either.)
    sa, sb = ra ** 2, rb ** 2

    # W_z = Integrate[ (b + d z)^2 / (2 pi mu_z), {s, sa, sb}, {z, za, zb} ]
    # Inner integral over z: Integrate[(b + d z)^2, z, za, zb]
    #     = b^2 (zb - za) + b d (zb^2 - za^2) + d^2 (zb^3 - za^3)/3
    # Outer over s: multiply by (sb - sa)/(2 pi mu_z)
    Iz_b2 = (sb - sa) * (zb - za) / (2 * np.pi * mu_z)
    Iz_bd = (sb - sa) * (zb ** 2 - za ** 2) / (2 * np.pi * mu_z)
    Iz_d2 = (sb - sa) * (zb ** 3 - za ** 3) / (3 * 2 * np.pi * mu_z)

    # Hessian in (b, c, d): b derivative twice, c twice, etc.
    # 2W_z = (sb - sa)/(pi mu_z) [ b^2 (zb - za) + b d (zb^2 - za^2) + d^2 (zb^3 - za^3)/3 ]
    # M_bd_z[b, b] = (sb - sa)*(zb - za)/(pi mu_z)
    # M_bd_z[b, d] = (sb - sa)(zb^2 - za^2)/(2 pi mu_z)  (1/2 from 2*coeff)
    # M_bd_z[d, d] = (sb - sa)(zb^3 - za^3)/(3 pi mu_z)
    # All others zero.
    M_z_coef = np.zeros((4, 4))  # rows/cols: (a, b, c, d)
    Hz_bb = (sb - sa) * (zb - za) / (np.pi * mu_z)
    Hz_bd = (sb - sa) * (zb ** 2 - za ** 2) / (2 * np.pi * mu_z)
    Hz_dd = (sb - sa) * (zb ** 3 - za ** 3) / (3 * np.pi * mu_z)
    M_z_coef[1, 1] = Hz_bb  # b-b
    M_z_coef[1, 3] = Hz_bd  # b-d
    M_z_coef[3, 1] = Hz_bd
    M_z_coef[3, 3] = Hz_dd

    # W_r = Integrate[ (c + d s)^2 / (8 pi mu_r s), {s, sa, sb}, {z, za, zb} ]
    # Inner over s: Integrate[(c + d s)^2 / s, sa, sb]
    #     = c^2 log(sb/sa) + 2 c d (sb - sa) + d^2 (sb^2 - sa^2)/2
    # Outer over z: multiply by (zb - za)/(8 pi mu_r)
    eps = 1e-30
    if sa < eps:
        # axis-touching: log(sb/sa) -> infinity, but coefficient c may be zero
        # for axis-symmetric solutions. Apply Limit treatment.
        # Actually: for axis-touching rectangle, the proper limit needs to keep
        # only finite pieces. We'll handle this via a special-case path below.
        log_sbsa = np.inf  # placeholder; should not be used
    else:
        log_sbsa = np.log(sb / sa)

    if sa < eps:
        # Axis-touching limit: c (linear in z direction) coefficient must be zero
        # for finite energy. The remaining contributions are from d (s*z mixed term).
        # Limit calculation: as sa -> 0, c^2 log(sb/sa) -> infinity unless c = 0.
        # Physically, on axis r=0, the basis function values (psi at axis nodes) must be zero
        # for finite A_phi (since A_phi = phi/(2*pi*r) and phi must be -> 0 at axis to keep A finite).
        # In V coordinates, V_j at axis nodes is constrained to 0 by Dirichlet BC.
        # The remaining d^2 contribution: limit_{sa->0} d^2 (sb^2 - sa^2)/2 = d^2 sb^2/2
        Hr_cc = 0.0  # excluded from active dofs by axis Dirichlet
        Hr_cd = 0.0
        # For d^2 term: 2 c d (sb - sa) goes to 2 c d sb as sa->0. If c is forced to 0, this is 0.
        # d^2 (sb^2 - sa^2)/2 -> d^2 sb^2 / 2
        Hr_dd = (zb - za) * (sb ** 2) / (8 * np.pi * mu_r)
    else:
        Hr_cc = (zb - za) * log_sbsa / (4 * np.pi * mu_r)
        Hr_cd = (zb - za) * (sb - sa) / (4 * np.pi * mu_r)
        Hr_dd = (zb - za) * (sb ** 2 - sa ** 2) / (8 * np.pi * mu_r)

    M_r_coef = np.zeros((4, 4))
    M_r_coef[2, 2] = Hr_cc  # c-c
    M_r_coef[2, 3] = Hr_cd  # c-d
    M_r_coef[3, 2] = Hr_cd
    M_r_coef[3, 3] = Hr_dd

    M_coef = M_z_coef + M_r_coef

    # Transform from (a, b, c, d) coords to (phi_0, phi_1, phi_2, phi_3) coords:
    # M_phi = inv_V^T . M_coef . inv_V   (since coeff = inv_V . phi)
    M_phi = inv_V.T @ M_coef @ inv_V

    # Transform from phi to V: phi_j = 2*pi*r_j*V_j, so V_j = phi_j/(2*pi*r_j)
    # M_V = T . M_phi . T   where T = diag(2*pi*r_j)
    r_nodes = np.array([ra, rb, rb, ra])
    T = np.diag(2 * np.pi * r_nodes)
    M_V = T @ M_phi @ T
    M_V = (M_V + M_V.T) / 2.0

    # Sign conventions: M_V is the FEMM-style (negative of physical PD stiffness).
    # In disk_hiruma we negate; here we return positive-definite directly so users
    # can use it for eigenvalue problems without thinking about sign.
    Mr_only = inv_V.T @ M_r_coef @ inv_V
    Mz_only = inv_V.T @ M_z_coef @ inv_V
    Mr_only = T @ Mr_only @ T
    Mz_only = T @ Mz_only @ T

    return M_V, Mr_only, Mz_only


def element_sigma_mass_quad(ra, rb, za, zb, sigma):
    """
    Sigma mass matrix for axis-aligned rectangle, V coords.

    M_sigma_ij V_i V_j = sigma * Integrate(A_phi^2 * 2 pi r dr dz)
                       = sigma/(8 pi^2) * Integrate(phi^2 / s ds dz)

    Closed-form derived from Mathematica.
    """
    sa, sb = ra ** 2, rb ** 2
    inv_V = _phi_basis_inverse(ra, rb, za, zb)

    # phi^2 = (a + b*s + c*z + d*s*z)^2
    # Expand and integrate (1/s)*phi^2 over [sa,sb] x [za,zb]
    # Coefficients of phi^2 grouped:
    #   a^2 / s              => log(sb/sa)
    #   2 a b                => (sb - sa)
    #   2 a c z / s          => z * log(sb/sa)
    #   2 a d z              => z * (sb - sa)
    #   b^2 s                => (sb^2 - sa^2)/2
    #   2 b c z              => z * (sb - sa)
    #   2 b d s z            => z * (sb^2 - sa^2)/2
    #   c^2 z^2 / s          => z^2 * log(sb/sa)
    #   2 c d z^2            => z^2 * (sb - sa)
    #   d^2 s z^2            => z^2 * (sb^2 - sa^2)/2

    eps = 1e-30
    log_sb_sa = np.log(sb / sa) if sa > eps else None
    Sb_Sa_1 = sb - sa
    Sb_Sa_2 = (sb ** 2 - sa ** 2) / 2

    # z-integrals: int_za^zb 1 dz = (zb - za); z dz = (zb^2-za^2)/2; z^2 dz = (zb^3-za^3)/3
    Iz0 = zb - za
    Iz1 = (zb ** 2 - za ** 2) / 2
    Iz2 = (zb ** 3 - za ** 3) / 3

    if sa > eps:
        # General case: all integrals finite
        # Hessian in (a, b, c, d): from 2*W_sigma differentiation
        # W_sigma = sigma/(8 pi^2) * sum_terms
        # Hessian H_ij = d^2 W / d(coef_i) d(coef_j) = (sigma/(4 pi^2)) * (2*coeff of coef_i*coef_j)
        # i.e., H = (sigma/(4 pi^2)) * (matrix of 2*[i==j coefficient] or [i != j coefficient])
        # Actually: W = sum_j sum_k coef_j * coef_k * I_jk where I is symmetric and coefficient
        # of coef_j*coef_k is summed (j != k), so M = 2 * I
        # Or: just derive H as second deriv

        # Coefficient matrix of phi^2: phi^2 = sum_{i,j} c_i c_j * (basis_i * basis_j)
        # where basis = (1, s, z, s*z). So phi^2 / s = sum_{i,j} c_i c_j * (b_i b_j / s)
        # And we need int/8pi^2 of phi^2/s over the rectangle.
        # Let I_ij = Integrate(b_i(s,z) b_j(s,z) / s, ds dz)
        # Then Hessian H_ij = sigma/(8 pi^2) * 2 * I_ij = sigma/(4 pi^2) I_ij  (with i,j on (a,b,c,d))

        I = np.zeros((4, 4))  # I[a,a]=int 1/s; I[a,b]=int s/s=1; etc.
        # b_0 = 1, b_1 = s, b_2 = z, b_3 = s*z
        # I[a,a] = integ_dz (zb-za) * integ_ds 1/s = (zb-za) * log_sb_sa
        I[0, 0] = Iz0 * log_sb_sa
        # I[a,b] = (zb-za) * Sb_Sa_1
        I[0, 1] = Iz0 * Sb_Sa_1
        I[1, 0] = I[0, 1]
        # I[a,c] = (z dz) * log_sb_sa = Iz1 * log_sb_sa
        I[0, 2] = Iz1 * log_sb_sa
        I[2, 0] = I[0, 2]
        # I[a,d] = z dz * Sb_Sa_1 = Iz1 * Sb_Sa_1
        I[0, 3] = Iz1 * Sb_Sa_1
        I[3, 0] = I[0, 3]
        # I[b,b] = (zb-za) * Sb_Sa_2
        I[1, 1] = Iz0 * Sb_Sa_2
        # I[b,c] = z dz * Sb_Sa_1 = Iz1 * Sb_Sa_1 (equals I[a,d]!)
        I[1, 2] = Iz1 * Sb_Sa_1
        I[2, 1] = I[1, 2]
        # I[b,d] = z dz * Sb_Sa_2 = Iz1 * Sb_Sa_2
        I[1, 3] = Iz1 * Sb_Sa_2
        I[3, 1] = I[1, 3]
        # I[c,c] = z^2 dz * log_sb_sa = Iz2 * log_sb_sa
        I[2, 2] = Iz2 * log_sb_sa
        # I[c,d] = z^2 dz * Sb_Sa_1 = Iz2 * Sb_Sa_1
        I[2, 3] = Iz2 * Sb_Sa_1
        I[3, 2] = I[2, 3]
        # I[d,d] = z^2 dz * Sb_Sa_2 = Iz2 * Sb_Sa_2
        I[3, 3] = Iz2 * Sb_Sa_2

        H_coef = sigma * I / (4 * np.pi)
    else:
        # Axis-touching: a, c coefficients are constrained to zero (from Dirichlet phi=0 at axis).
        # We zero those rows/cols and use limit for the rest.
        # Limit sa -> 0:
        #   log(sb/sa) terms: filtered out (a, c -> 0)
        #   (sb - sa) -> sb
        #   (sb^2 - sa^2)/2 -> sb^2/2
        I = np.zeros((4, 4))
        # Only b, d contribute (since a, c forced to 0)
        # I[b,b] = Iz0 * sb^2 / 2
        I[1, 1] = Iz0 * (sb ** 2) / 2
        # I[b,d] = Iz1 * sb^2 / 2
        I[1, 3] = Iz1 * (sb ** 2) / 2
        I[3, 1] = I[1, 3]
        # I[d,d] = Iz2 * sb^2 / 2
        I[3, 3] = Iz2 * (sb ** 2) / 2
        # b-d terms also need (sb - sa)*z components — but those go through (b,d) to (b,d) via
        # multiplication, all already accounted for above.
        # Note: zero the (a, c) rows/cols ensures Dirichlet-axis values not contribute
        H_coef = sigma * I / (4 * np.pi)

    # Transform to phi nodal values, then to V
    M_phi = inv_V.T @ H_coef @ inv_V
    r_nodes = np.array([ra, rb, rb, ra])
    T = np.diag(2 * np.pi * r_nodes)
    M_V = T @ M_phi @ T
    M_V = (M_V + M_V.T) / 2.0
    return M_V


def structured_disk_mesh(R_disk, t_disk, R_air, Z_air, NR_disk, Nz_disk, NR_air_extra, Nz_air_extra):
    """
    Build a structured rectangular mesh for axisymmetric Cu disk in vacuum.

    Returns (nodes, quads, materials):
      nodes:     (Nnode, 2) array of (r, z) coordinates
      quads:     (Nelem, 4) integer indices, ordered (sa,za),(sb,za),(sb,zb),(sa,zb)
      materials: (Nelem,) array of strings ("conductor", "air")

    Disk: r in [0, R_disk], z in [-t_disk/2, t_disk/2]   NR_disk x Nz_disk cells
    Air:  the rest, extending to R_air, |z| < Z_air
          NR_air_extra columns added to the right of disk; Nz_air_extra rows above/below.
    """
    # r grid: NR_disk cells inside disk, NR_air_extra cells in air
    r_disk = np.linspace(0, R_disk, NR_disk + 1)
    # geometric expansion for air (denser near disk)
    r_air = np.geomspace(R_disk, R_air, NR_air_extra + 1)[1:]  # skip duplicate R_disk
    r = np.concatenate([r_disk, r_air])

    # z grid: split into bottom-air, disk, top-air
    z_disk = np.linspace(-t_disk / 2, t_disk / 2, Nz_disk + 1)
    z_air_above = np.geomspace(t_disk / 2 + 1e-6, Z_air, Nz_air_extra + 1)[1:]
    z_air_below = -np.geomspace(t_disk / 2 + 1e-6, Z_air, Nz_air_extra + 1)[1:][::-1]
    z = np.concatenate([z_air_below, z_disk, z_air_above])

    NR = len(r) - 1
    Nz = len(z) - 1

    # nodes
    nodes = np.array([(rr, zz) for zz in z for rr in r])
    Nnode = nodes.shape[0]

    def nid(i, j):  # i in r, j in z
        return j * (NR + 1) + i

    quads = []
    materials = []
    for j in range(Nz):
        for i in range(NR):
            n0 = nid(i, j)
            n1 = nid(i + 1, j)
            n2 = nid(i + 1, j + 1)
            n3 = nid(i, j + 1)
            quads.append([n0, n1, n2, n3])
            # Material: conductor if cell is fully in disk
            r_a = r[i]
            r_b = r[i + 1]
            z_a = z[j]
            z_b = z[j + 1]
            r_center = (r_a + r_b) / 2
            z_center = (z_a + z_b) / 2
            if r_center < R_disk and abs(z_center) < t_disk / 2:
                materials.append("conductor")
            else:
                materials.append("air")
    return nodes, np.array(quads), np.array(materials)


def assemble_quad(nodes, quads, mat_mu_r, mat_mu_z):
    """Assemble global K matrix from Q1 element matrices."""
    Nnode = nodes.shape[0]
    rows, cols, vals = [], [], []
    for ie in range(quads.shape[0]):
        idx = quads[ie]
        r0, z0 = nodes[idx[0]]
        r1, z1 = nodes[idx[1]]
        r2, z2 = nodes[idx[2]]
        r3, z3 = nodes[idx[3]]
        ra, rb = r0, r1
        za, zb = z0, z2
        if abs(r0 - r3) > 1e-12 or abs(r1 - r2) > 1e-12:
            raise ValueError("Quad not axis-aligned (r-edges differ)")
        if abs(z0 - z1) > 1e-12 or abs(z2 - z3) > 1e-12:
            raise ValueError("Quad not axis-aligned (z-edges differ)")
        M_e, _, _ = element_matrices_quad(ra, rb, za, zb, mat_mu_r[ie], mat_mu_z[ie])
        for j in range(4):
            for k in range(4):
                rows.append(idx[j])
                cols.append(idx[k])
                vals.append(M_e[j, k])
    return sp.csr_matrix((vals, (rows, cols)), shape=(Nnode, Nnode))


def assemble_sigma_mass_quad(nodes, quads, sigma_per_elem):
    """Assemble sigma mass matrix from Q1 element sigma mass."""
    Nnode = nodes.shape[0]
    rows, cols, vals = [], [], []
    for ie in range(quads.shape[0]):
        if sigma_per_elem[ie] == 0.0:
            continue
        idx = quads[ie]
        r0, z0 = nodes[idx[0]]
        r1, z1 = nodes[idx[1]]
        r2, z2 = nodes[idx[2]]
        ra, rb = r0, r1
        za, zb = z0, z2
        M_e = element_sigma_mass_quad(ra, rb, za, zb, sigma_per_elem[ie])
        for j in range(4):
            for k in range(4):
                rows.append(idx[j])
                cols.append(idx[k])
                vals.append(M_e[j, k])
    return sp.csr_matrix((vals, (rows, cols)), shape=(Nnode, Nnode))


def _self_test():
    """Sanity test on a single non-axis quadrilateral."""
    print("=== axifemm_quad self-test ===")
    ra, rb, za, zb = 1e-3, 2e-3, 0.0, 1e-3
    mu0 = 4 * np.pi * 1e-7
    M_e, Mr, Mz = element_matrices_quad(ra, rb, za, zb, mu0, mu0)
    print("M_e (4x4):")
    print(M_e)
    eigs = np.linalg.eigvalsh((M_e + M_e.T) / 2.0)
    print("Eigenvalues:", eigs)
    # Test: constant-phi (V_j = const/(2*pi*r_j)) should give zero energy
    # phi_j = 2*pi*r_j*V_j = const. So V_j = const/(2*pi*r_j) = const/(2*pi)*[1/ra, 1/rb, 1/rb, 1/ra]
    test_v = np.array([1/ra, 1/rb, 1/rb, 1/ra])
    energy = test_v @ M_e @ test_v
    print(f"V=(1/r) energy = {energy:.6e} (should be ~0)")

    # Sigma mass
    M_sig = element_sigma_mass_quad(ra, rb, za, zb, 5.8e7)
    print("M_sigma (4x4):")
    print(M_sig)
    eigs_sig = np.linalg.eigvalsh((M_sig + M_sig.T) / 2.0)
    print("Sigma eigvals:", eigs_sig)


if __name__ == "__main__":
    _self_test()
