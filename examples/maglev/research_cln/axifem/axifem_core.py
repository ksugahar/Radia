"""
axifem_core.py — Phase 1 prototype: FEMM-style axisymmetric magnetostatic solver.

Translates the StaticAxisymmetric() function of FEMM's prob3big.cpp
(David Meeker) directly to NumPy/scipy.

Key idea (Henrotte/Meeker):
- 3 nodes per triangle (P1-equivalent DOF count)
- Shape functions are linear in {1, r^2, z}, NOT linear in {1, r, z}
- Yields B_z constant per element, B_r ~ 1/r per element
- Removes 1/r singularities at axis r=0

UNKNOWN VARIABLE CONVENTION (matches prob3big.cpp):
We solve for V_j at each node, where V is FEMM's internal flux variable.
The flux phi = 2*pi*r*A_phi is interpolated linearly in (s=r^2, z) coordinates,
and at each node phi_j = 2*pi*r_j*mu_0*V_j (FEMM convention; V is roughly H_phi).
Equivalently: A_phi_j = mu_0 * V_j.

Units in this prototype: pure SI (m, Tesla, H/m, A/m).
Pass mu_r, mu_z per element as ABSOLUTE permeability (H/m), not relative.
The element matrices follow prob3big.cpp formulas literally; the FEMM
"c = 4*pi*1e-5" cm-unit factor is NOT applied here.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EPS_AXIS = 1.0e-12
MU0 = 4.0e-7 * np.pi


# ---------------------------------------------------------------------------
# Element-level routines (translation of prob3big.cpp lines 173-273)
# ---------------------------------------------------------------------------


def _element_geom(rn: np.ndarray, zn: np.ndarray):
    """Geometric coefficients for a triangle with vertex coords (rn[i], zn[i]).

    Returns (p, q, g, a, R, a_hat, R_hat) following prob3big.cpp notation:
      p[i]  = z_{i+1} - z_{i+2}     (b's in Allaire/Meeker)
      q[i]  = r_{i+2} - r_{i+1}     (c's in Allaire/Meeker)
      g[i]  = midpoint r of edge opposite to node i
      a     = standard P1 triangle area
      R     = centroid r
      a_hat = "modified" area for {1, r^2, z} basis (= d_meeker / (4 R))
      R_hat = effective radial average for the z-component matrix
              (with axis-touching cases handled).
    """
    p = np.array([zn[1] - zn[2], zn[2] - zn[0], zn[0] - zn[1]])
    q = np.array([rn[2] - rn[1], rn[0] - rn[2], rn[1] - rn[0]])
    g = np.array([(rn[1] + rn[2]) / 2, (rn[0] + rn[2]) / 2, (rn[0] + rn[1]) / 2])

    a = (p[0] * q[1] - p[1] * q[0]) / 2.0
    R = (rn[0] + rn[1] + rn[2]) / 3.0
    a_hat = sum(rn[j] ** 2 * p[j] / (4.0 * R) for j in range(3))

    flag = sum(1 for j in range(3) if abs(rn[j]) < EPS_AXIS)

    if flag == 2:
        R_hat = R
    elif flag == 1:
        # which node is on the axis?
        if abs(rn[0]) < EPS_AXIS:
            r1, r2 = rn[1], rn[2]
        elif abs(rn[1]) < EPS_AXIS:
            r1, r2 = rn[2], rn[0]
        else:
            r1, r2 = rn[0], rn[1]
        if abs(r1 - r2) < EPS_AXIS:
            R_hat = r2 / 2.0
        else:
            R_hat = (r1 - r2) / (2.0 * np.log(r1) - 2.0 * np.log(r2))
    else:
        # interior; q[i] == 0 sub-cases first
        if abs(q[0]) < EPS_AXIS:
            R_hat = q[1] ** 2 / (2.0 * (-q[1] + rn[0] * np.log(rn[0] / rn[2])))
        elif abs(q[1]) < EPS_AXIS:
            R_hat = q[2] ** 2 / (2.0 * (-q[2] + rn[1] * np.log(rn[1] / rn[0])))
        elif abs(q[2]) < EPS_AXIS:
            R_hat = q[0] ** 2 / (2.0 * (-q[0] + rn[2] * np.log(rn[2] / rn[1])))
        else:
            R_hat = -(q[0] * q[1] * q[2]) / (
                2.0 * (
                    q[0] * rn[0] * np.log(rn[0])
                    + q[1] * rn[1] * np.log(rn[1])
                    + q[2] * rn[2] * np.log(rn[2])
                )
            )
    return p, q, g, a, R, a_hat, R_hat


def element_matrices(rn: np.ndarray, zn: np.ndarray, mu_r: float, mu_z: float):
    """Return (M_e, Mr, Mz, p, q, g, a, R, a_hat, R_hat).

    M_e = Mr/mu_z + Mz/mu_r  (FEMM convention from prob3big.cpp line 539)
    """
    p, q, g, a, R, a_hat, R_hat = _element_geom(rn, zn)

    # --- Mr (Mx in FEMM): contribution from B_r component ---
    K_r = -1.0 / (2.0 * a_hat * R)
    Mr = np.zeros((3, 3))
    for j in range(3):
        for k in range(j, 3):
            Mr[j, k] = K_r * p[j] * rn[j] * p[k] * rn[k]
            Mr[k, j] = Mr[j, k]
    # axis-singular diagonal boost (prob3big.cpp 250-251)
    diag_sum = Mr[0, 0] + Mr[1, 1] + Mr[2, 2]
    for j in range(3):
        if abs(rn[j]) < EPS_AXIS:
            Mr[j, j] += diag_sum

    # --- Mz (My in FEMM): contribution from B_z component ---
    K_z = -1.0 / (2.0 * a_hat * R_hat)
    Mz = np.zeros((3, 3))
    for j in range(3):
        for k in range(j, 3):
            Mz[j, k] = K_z * (q[j] * rn[j]) * (q[k] * rn[k]) * (g[j] / R) * (g[k] / R)
            Mz[k, j] = Mz[j, k]

    # FEMM: Me = Mx/mu2 + My/mu1 (line 539);  mu1=mu_r, mu2=mu_z
    M_e = Mr / mu_z + Mz / mu_r

    return M_e, Mr, Mz, p, q, g, a, R, a_hat, R_hat


# ---------------------------------------------------------------------------
# Source contributions
# ---------------------------------------------------------------------------


def magnetization_rhs(rn, zn, M_dir_deg: float, H_c: float):
    """Per-element RHS contribution from magnetization (constant H_c, direction MagDir).

    Translates prob3big.cpp lines 314-352 (volume divergence -> edge integral).

    Parameters
    ----------
    rn, zn : node r, z coordinates (length 3)
    M_dir_deg : magnetization direction angle in degrees (FEMM convention,
                0 = +r, 90 = +z)
    H_c : magnitude of H_c in A/m (= |M| for permanent magnet using B = mu0(H+M))
    """
    t = M_dir_deg
    be = np.zeros(3)
    cos_t = np.cos(t * np.pi / 180.0)
    sin_t = np.sin(t * np.pi / 180.0)
    for j in range(3):
        k = (j + 1) % 3
        r_mid = (rn[j] + rn[k]) / 2.0
        # In SI: K = -mu0 * r_mid * H_c * (cos*(r_k - r_j) + sin*(z_k - z_j))
        # FEMM uses 0.0001 (cm units); we drop that and absorb mu0 elsewhere.
        K = -r_mid * H_c * (cos_t * (rn[k] - rn[j]) + sin_t * (zn[k] - zn[j]))
        be[j] += K
        be[k] += K
    return be


def current_rhs(rn, zn, J_phi: float):
    """Per-element RHS contribution from prescribed J_phi (A/m^2).

    Translates prob3big.cpp lines 296-308 (volume integral of J).
    For phi = 2*pi*r*A formulation, the contribution is -(2*pi)*R*J_phi*a/3 per node.
    """
    p = np.array([zn[1] - zn[2], zn[2] - zn[0], zn[0] - zn[1]])
    q = np.array([rn[2] - rn[1], rn[0] - rn[2], rn[1] - rn[0]])
    a = (p[0] * q[1] - p[1] * q[0]) / 2.0
    R = (rn[0] + rn[1] + rn[2]) / 3.0
    K = -2.0 * np.pi * R * J_phi * a / 3.0
    return np.array([K, K, K])


# ---------------------------------------------------------------------------
# Global assembly
# ---------------------------------------------------------------------------


def assemble(
    nodes: np.ndarray,         # (Nnode, 2): each row (r, z)
    triangles: np.ndarray,     # (Nelem, 3): node indices
    mat_mu_r: np.ndarray,      # (Nelem,) radial permeability per element [H/m]
    mat_mu_z: np.ndarray,      # (Nelem,) axial permeability per element  [H/m]
    mat_H_c: np.ndarray | None = None,
    mat_M_dir_deg: np.ndarray | None = None,
    mat_Jphi: np.ndarray | None = None,
):
    """Assemble global stiffness K and RHS b. Returns (K_csr, b)."""
    Nnode = nodes.shape[0]
    Nelem = triangles.shape[0]

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    b = np.zeros(Nnode)

    for ie in range(Nelem):
        idx = triangles[ie]
        rn = nodes[idx, 0]
        zn = nodes[idx, 1]
        M_e, *_ = element_matrices(rn, zn, mat_mu_r[ie], mat_mu_z[ie])
        # scatter (note FEMM accumulates with a sign flip; keep positive here)
        for j in range(3):
            for k in range(3):
                rows.append(idx[j])
                cols.append(idx[k])
                vals.append(M_e[j, k])

        # source contributions
        if mat_H_c is not None and mat_H_c[ie] != 0.0:
            be = magnetization_rhs(rn, zn, mat_M_dir_deg[ie], mat_H_c[ie])
            for j in range(3):
                b[idx[j]] += be[j]
        if mat_Jphi is not None and mat_Jphi[ie] != 0.0:
            bj = current_rhs(rn, zn, mat_Jphi[ie])
            for j in range(3):
                b[idx[j]] += bj[j]

    K = sp.csr_matrix((vals, (rows, cols)), shape=(Nnode, Nnode))
    return K, b


def apply_dirichlet(K: sp.csr_matrix, b: np.ndarray, dirichlet_nodes: np.ndarray, values: np.ndarray):
    """Apply Dirichlet BC: set V[i] = values[i] for i in dirichlet_nodes.

    Modifies K, b in-place (K becomes a copy), returns (K_bc, b_bc).
    """
    K = K.tolil()
    b = b.copy()
    for nid, val in zip(dirichlet_nodes, values):
        b -= np.asarray(K[:, nid].toarray()).flatten() * val
        K[nid, :] = 0
        K[:, nid] = 0
        K[nid, nid] = 1.0
        b[nid] = val
    return K.tocsr(), b


def solve(K: sp.csr_matrix, b: np.ndarray):
    """Direct sparse solve. For NMR-size problems this is fine."""
    return spla.spsolve(K, b)


# ---------------------------------------------------------------------------
# Postprocessing
# ---------------------------------------------------------------------------


def field_in_element(rn, zn, phi_e, r_pt, z_pt):
    """Evaluate (A_phi, B_r, B_z) at a point inside the triangle.

    phi_e: nodal values of phi = 2*pi*r*A_phi at the 3 vertices.
    Uses the {1, r^2, z} interpolation:
        phi = c0 + c1*s + c2*z, with s = r^2.
    """
    s_n = rn ** 2
    M = np.column_stack([np.ones(3), s_n, zn])  # 3x3
    coeffs = np.linalg.solve(M, phi_e)          # [c0, c1, c2]
    c0, c1, c2 = coeffs

    s_pt = r_pt ** 2
    phi = c0 + c1 * s_pt + c2 * z_pt
    A_phi = phi / (2.0 * np.pi * r_pt) if r_pt > EPS_AXIS else 0.0
    B_z = c1 / np.pi
    B_r = -c2 / (2.0 * np.pi * r_pt) if r_pt > EPS_AXIS else 0.0
    return A_phi, B_r, B_z


def find_element(nodes, triangles, r_pt, z_pt):
    """Find the index of the triangle containing (r_pt, z_pt). Brute force."""
    for ie in range(triangles.shape[0]):
        idx = triangles[ie]
        r0, z0 = nodes[idx[0]]
        r1, z1 = nodes[idx[1]]
        r2, z2 = nodes[idx[2]]
        # barycentric test in (r, z)
        denom = (z1 - z2) * (r0 - r2) + (r2 - r1) * (z0 - z2)
        if abs(denom) < 1e-30:
            continue
        a = ((z1 - z2) * (r_pt - r2) + (r2 - r1) * (z_pt - z2)) / denom
        b = ((z2 - z0) * (r_pt - r2) + (r0 - r2) * (z_pt - z2)) / denom
        c = 1.0 - a - b
        if a >= -1e-9 and b >= -1e-9 and c >= -1e-9:
            return ie
    return -1


# ---------------------------------------------------------------------------
# Self-test (single-triangle sanity)
# ---------------------------------------------------------------------------


def _self_test():
    """Sanity test on a single triangle."""
    rn = np.array([0.0, 1.0e-3, 0.0])
    zn = np.array([0.0, 0.0, 1.0e-3])
    M_e, Mr, Mz, p, q, g, a, R, a_hat, R_hat = element_matrices(rn, zn, MU0, MU0)
    print("p =", p)
    print("q =", q)
    print("g =", g)
    print(f"a={a:.4e}  R={R:.4e}  a_hat={a_hat:.4e}  R_hat={R_hat:.4e}")
    print("Mr =\n", Mr)
    print("Mz =\n", Mz)
    print("M_e =\n", M_e)


if __name__ == "__main__":
    _self_test()
