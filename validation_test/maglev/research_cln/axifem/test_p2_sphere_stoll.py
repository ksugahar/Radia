"""test_p2_sphere_stoll.py — Phase B1b Item 3
P2 sphere eigenvalue compared with Stoll analytical (finite air box).

Setup: Cu sphere R=10 mm in vacuum, sigma=5.8e7 S/m. Modelled in
axisymmetric (rho, theta) -> (r, z) coordinates with structured
triangulation. Air box of large radius R_air=50 mm (5x sphere radius)
serves as a finite approximation to infinity (proper Kelvin
transformation deferred to Phase B3).

Stoll analytical Foster spectrum for Cu sphere:
    tau_n = mu_0 sigma a^2 / (n pi)^2     (n = 1, 2, 3, ...)
For Cu sphere R=10 mm:  tau_1 = 4 pi * 1e-7 * 5.8e7 * (0.01)^2 / pi^2
                              = 738.5 us  (leading)

Reference baseline:
- axifem Q1 + 500 mm air box: tau_pair[0] = leading 0.21% off Stoll.
- VIM DD axisym 270 cells: leading tau_pair 0.2% off Stoll.

We solve the generalized eigenvalue problem with Schur-complement
reduction onto conductor DOFs:
    S u_c = lambda M_cc u_c,   S = K_cc - K_ca K_aa^{-1} K_ac

(M is non-zero only on conductor DOFs, hence semi-definite; the
Schur reduction eliminates the air block which has zero M.)

The leading eigenvalue gives the slowest-decaying eddy-current mode.
"""
from __future__ import annotations

import sys
sys.stdout.reconfigure(encoding="utf-8")

import time

import numpy as np
import scipy.linalg as sla

sys.path.insert(0, str(__file__.rsplit("\\", 1)[0] if "\\" in __file__
                       else __file__.rsplit("/", 1)[0]))

from axifem_p2_triangle import AxifemP2Triangle, MU0


R_SPHERE = 10.0e-3   # 10 mm
R_AIR = 50.0e-3      # 50 mm = 5x sphere (finite air-box approximation)
SIGMA = 5.8e7
EPS_AXIS = 1.0e-12


def stoll_tau_n(n, a=R_SPHERE, sigma=SIGMA):
    """Stoll Foster spectrum: tau_n = mu_0 sigma a^2 / (n pi)^2."""
    return MU0 * sigma * a ** 2 / (n * np.pi) ** 2


def generate_sphere_mesh(N_rho_in, N_rho_out, N_theta):
    """Structured (rho, theta) mesh of sphere + air, mapped to (r, z).

    Args:
        N_rho_in: radial cells inside conductor (rho in [0, R_sphere])
        N_rho_out: radial cells outside conductor (rho in [R_sphere, R_air])
        N_theta: meridional cells (theta in [0, pi])

    Returns:
        nodes_p2: (N, 2) array of (r, z) coords
        triangles_p2: (M, 6) array of indices [v0, v1, v2, m01, m12, m20]
        is_conductor: (M,) array of bool — True for cells inside sphere
    """
    rho_in = np.linspace(0.0, R_SPHERE, N_rho_in + 1)
    rho_out = np.linspace(R_SPHERE, R_AIR, N_rho_out + 1)
    # Combine (avoid duplicate at R_SPHERE)
    rho_grid = np.concatenate([rho_in, rho_out[1:]])
    theta_grid = np.linspace(0.0, np.pi, N_theta + 1)
    N_rho = len(rho_grid) - 1  # number of radial cells

    # Vertex coords on (rho, theta) grid
    nodes_p1 = []
    for j, theta in enumerate(theta_grid):
        for i, rho in enumerate(rho_grid):
            r = rho * np.sin(theta)
            z = rho * np.cos(theta)
            # Clean up axis-touching: snap small r to 0
            if abs(r) < EPS_AXIS:
                r = 0.0
            nodes_p1.append([r, z])
    nodes_p1 = np.array(nodes_p1)
    n_per_row = len(rho_grid)

    triangles_p1 = []
    is_cond_p1 = []
    for j in range(N_theta):
        for i in range(N_rho):
            n00 = j * n_per_row + i
            n10 = j * n_per_row + i + 1
            n01 = (j + 1) * n_per_row + i
            n11 = (j + 1) * n_per_row + i + 1
            # Skip degenerate cells where all 3 vertices on axis r=0
            # (i.e., both rho_a=0 and theta=0 or pi). For our grid only
            # rho=0 + any theta. Such pencils don't form triangles meaningfully.
            # The triangulation below (anti-diagonal split) creates one tri
            # that has 2 axis-touching vertices and 1 off-axis — OK in
            # Henrotte basis as long as Vandermonde is non-degenerate.
            t1 = (n00, n10, n11)
            t2 = (n00, n11, n01)
            for tri in (t1, t2):
                rho_center = np.mean([rho_grid[k]
                                      for k in [
                                          i if tri[0] in [n00, n01] else i + 1,
                                          i if tri[1] in [n00, n01] else i + 1,
                                          i if tri[2] in [n00, n01] else i + 1,
                                      ]])
                is_cond = rho_center < R_SPHERE
                triangles_p1.append(tri)
                is_cond_p1.append(is_cond)
    triangles_p1 = np.array(triangles_p1)
    is_conductor_per_cell = np.array(is_cond_p1)

    # P2 nodes: append mid-edge nodes
    edge_to_mid = {}
    nodes_p2_list = nodes_p1.tolist()
    triangles_p2 = []
    for tri in triangles_p1:
        v0, v1, v2 = tri

        def get_mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in edge_to_mid:
                mid = 0.5 * (nodes_p1[a] + nodes_p1[b])
                if abs(mid[0]) < EPS_AXIS:
                    mid = np.array([0.0, mid[1]])
                edge_to_mid[key] = len(nodes_p2_list)
                nodes_p2_list.append(mid.tolist())
            return edge_to_mid[key]

        m01 = get_mid(v0, v1)
        m12 = get_mid(v1, v2)
        m20 = get_mid(v2, v0)
        triangles_p2.append([v0, v1, v2, m01, m12, m20])

    nodes_p2 = np.array(nodes_p2_list)
    triangles_p2 = np.array(triangles_p2)
    return nodes_p2, triangles_p2, is_conductor_per_cell


def boundary_outer_nodes(nodes, tol=1e-9):
    """Find nodes on outer air boundary (rho ~ R_AIR) for Dirichlet BC."""
    r = nodes[:, 0]
    z = nodes[:, 1]
    rho = np.sqrt(r * r + z * z)
    return np.where(np.abs(rho - R_AIR) < tol)[0]


def conductor_dofs(nodes, triangles, is_cond_per_cell, tol=EPS_AXIS):
    """Return indices of DOFs that touch any conductor cell."""
    cond_set = set()
    for tri_idx, tri in enumerate(triangles):
        if is_cond_per_cell[tri_idx]:
            for idx in tri:
                cond_set.add(int(idx))
    return np.array(sorted(cond_set))


def assemble_p2_split(nodes, triangles, is_cond_per_cell, sigma):
    """Assemble global K (everywhere) and M (only on conductor cells).

    M is non-zero only on conductor cells (sigma * mass), so it is
    semi-definite over the full DOF set. We will use Schur complement
    later to extract a positive-definite block.
    """
    N = nodes.shape[0]
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    skipped_degenerate = 0
    for tri_idx, tri in enumerate(triangles):
        rn = nodes[tri, 0]
        zn = nodes[tri, 1]
        # Skip pencil cells where all 6 nodes collapse onto axis (rare)
        if np.all(rn < EPS_AXIS):
            skipped_degenerate += 1
            continue
        try:
            t_elem = AxifemP2Triangle(rn, zn)
        except RuntimeError:
            skipped_degenerate += 1
            continue
        # Skip degenerate (collapsed-area) triangles
        if t_elem.area < 1e-20:
            skipped_degenerate += 1
            continue
        Ke = t_elem.stiffness(mu_r=MU0)
        if is_cond_per_cell[tri_idx]:
            Me = t_elem.sigma_mass(sigma=sigma)
        else:
            Me = np.zeros((6, 6))
        for i_l in range(6):
            for j_l in range(6):
                K[tri[i_l], tri[j_l]] += Ke[i_l, j_l]
                M[tri[i_l], tri[j_l]] += Me[i_l, j_l]
    if skipped_degenerate:
        print(f"  (skipped {skipped_degenerate} degenerate cells)")
    return K, M


def solve_schur_eigenvalue(K, M, fixed_idx, cond_idx, nodes, n_eigs=5):
    """Schur-complement Cauer extraction with Dirichlet BC.

    Free DOFs := all DOFs minus (outer-boundary + axis r=0 DOFs).
    Axis DOFs are auto-decoupled by FEMM Henrotte r_i factor (zero rows in
    K, M) and so are pinned to 0 to avoid singular K_aa in the Schur step.

    Within free DOFs:
        cond := conductor DOFs (subset of free)
        air  := free DOFs not on any conductor cell
    Schur:
        S = K_cc - K_ca K_aa^{-1} K_ac
    where indices are in the free subset.

    Solve S u_c = lambda M_cc u_c (lambda = 1/tau_decay).
    """
    N = K.shape[0]
    axis_idx = np.where(nodes[:, 0] < EPS_AXIS)[0]
    pinned = np.union1d(fixed_idx, axis_idx)
    free_idx = np.setdiff1d(np.arange(N), pinned)
    cond_in_free = np.array([i for i, idx in enumerate(free_idx)
                              if idx in set(cond_idx.tolist())])
    air_in_free = np.array([i for i in range(len(free_idx))
                             if i not in set(cond_in_free.tolist())])

    K_free = K[np.ix_(free_idx, free_idx)]
    M_free = M[np.ix_(free_idx, free_idx)]

    if len(air_in_free) > 0:
        K_cc = K_free[np.ix_(cond_in_free, cond_in_free)]
        K_ca = K_free[np.ix_(cond_in_free, air_in_free)]
        K_ac = K_free[np.ix_(air_in_free, cond_in_free)]
        K_aa = K_free[np.ix_(air_in_free, air_in_free)]
        X = sla.solve(K_aa, K_ac, assume_a="sym")
        S = K_cc - K_ca @ X
    else:
        S = K_free[np.ix_(cond_in_free, cond_in_free)]
    M_cc = M_free[np.ix_(cond_in_free, cond_in_free)]
    S = 0.5 * (S + S.T)
    M_cc = 0.5 * (M_cc + M_cc.T)

    # Diagnostic: M_cc rank deficiency check + drop zero-mass rows.
    diag = np.diag(M_cc)
    keep = diag > 1e-30 * diag.max()
    n_dropped = (~keep).sum()
    if n_dropped > 0:
        S = S[np.ix_(keep, keep)]
        M_cc = M_cc[np.ix_(keep, keep)]
    n_total = S.shape[0]
    n_take = min(n_eigs, n_total)
    eigvals = sla.eigvalsh(S, M_cc, subset_by_index=[0, n_take - 1])
    pos = eigvals[eigvals > 0]
    return pos, n_dropped


def main():
    print("=" * 64)
    print("Phase B1b Item 3: P2 sphere vs Stoll analytical")
    print("=" * 64)
    print(f"\nCu sphere R={R_SPHERE*1e3:.1f} mm, sigma={SIGMA:.2e} S/m")
    print(f"Air box R_air={R_AIR*1e3:.1f} mm (Dirichlet A_phi=0 on outer)")
    tau_stoll_1 = stoll_tau_n(1)
    tau_stoll_2 = stoll_tau_n(2)
    tau_stoll_3 = stoll_tau_n(3)
    print(f"\nStoll analytical Foster spectrum (sphere in infinite vacuum):")
    print(f"  tau_1 = {tau_stoll_1*1e6:.3f} us")
    print(f"  tau_2 = {tau_stoll_2*1e6:.3f} us")
    print(f"  tau_3 = {tau_stoll_3*1e6:.3f} us")

    # Mesh sweep
    levels = [
        (4, 4, 12),
        (6, 6, 16),
        (8, 8, 24),
        (10, 10, 30),
    ]
    print(f"\nMesh sweep (N_rho_in / N_rho_out / N_theta):")
    print(f"  {'(in, out, th)':>15s}  {'N_p2':>6s}  {'tau_1 [us]':>12s}  "
          f"{'gap [%]':>8s}  {'assembly [s]':>12s}")

    for N_in, N_out, N_th in levels:
        t0 = time.time()
        nodes, triangles, is_cond_cell = generate_sphere_mesh(N_in, N_out, N_th)
        bd_idx = boundary_outer_nodes(nodes)
        cond_idx = conductor_dofs(nodes, triangles, is_cond_cell)
        K, M = assemble_p2_split(nodes, triangles, is_cond_cell, SIGMA)
        t_assemble = time.time() - t0

        try:
            evals, n_dropped = solve_schur_eigenvalue(K, M, bd_idx, cond_idx,
                                                       nodes, n_eigs=3)
            if len(evals) > 0:
                lam_1 = evals[0]
                tau_1 = 1.0 / lam_1
                gap = (tau_1 - tau_stoll_1) / tau_stoll_1 * 100
                print(f"  ({N_in:>2d}, {N_out:>2d}, {N_th:>2d})    "
                      f"{nodes.shape[0]:>6d}  {tau_1*1e6:>12.4f}  "
                      f"{gap:>+8.3f}  {t_assemble:>12.2f}  "
                      f"(dropped {n_dropped} zero-mass DOFs)")
            else:
                print(f"  ({N_in:>2d}, {N_out:>2d}, {N_th:>2d})    no positive eigvals")
        except Exception as e:
            print(f"  ({N_in:>2d}, {N_out:>2d}, {N_th:>2d})    FAIL: {e}")

    print()
    print("Reference (for context):")
    print(f"  axifem Q1 (axis-aligned + 500mm air box): 0.21% off Stoll")
    print(f"  VIM DD axisym 270 cells reduced:           0.20% off Stoll")
    print()
    print("Note: This is finite air box, NOT proper Kelvin transformation.")
    print("Phase B3 (Kelvin extension) is the next step.")


if __name__ == "__main__":
    main()
