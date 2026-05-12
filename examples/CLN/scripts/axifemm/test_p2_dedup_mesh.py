"""Bisect 2: build the dedup'd sphere mesh in pure Python (no NGSolve),
assemble via AxifemmP2Triangle, run Python Schur. If result matches the
non-dedup'd Python ref (737 us), then dedup mesh is fine and bug is C++-side
NGSolve assembly. If it gives 50% off like C++ test, dedup mesh itself is
wrong.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

from math import pi
import numpy as np
import scipy.linalg as sla

sys.path.insert(0, ".")
from axifemm_p2_triangle import AxifemmP2Triangle, MU0

R_SPHERE = 10.0e-3
R_AIR = 50.0e-3
SIGMA = 5.8e7


def build_dedup_mesh(N_rho_in, N_rho_out, N_theta):
    """Same algorithm as test_p2_cpp_sphere.make_sphere_axisym_trig_mesh but
    in pure Python: axis-vertex dedup + degenerate-pnum filtering. Returns
    (nodes_p2, triangles_p2, is_cond_cell, outer_boundary_idx)."""
    rho_in = np.linspace(0.0, R_SPHERE, N_rho_in + 1)
    rho_out = np.linspace(R_SPHERE, R_AIR, N_rho_out + 1)
    rho_grid = np.concatenate([rho_in, rho_out[1:]])
    theta_grid = np.linspace(0.0, pi, N_theta + 1)
    n_rho = len(rho_grid)

    # Vertex dedup: rho=0 collapse by z key.
    nodes_p1 = []  # list of [r, z]
    pnums = [None] * (n_rho * (N_theta + 1))
    axis_z_to_pnum = {}

    def add_or_get(r, z):
        if abs(r) < 1e-12:
            key = round(z, 15)
            if key in axis_z_to_pnum:
                return axis_z_to_pnum[key]
            idx = len(nodes_p1)
            nodes_p1.append([0.0, z])
            axis_z_to_pnum[key] = idx
            return idx
        idx = len(nodes_p1)
        nodes_p1.append([r, z])
        return idx

    for j, theta in enumerate(theta_grid):
        for i, rho in enumerate(rho_grid):
            r = rho * np.sin(theta)
            z = rho * np.cos(theta)
            pnums[j * n_rho + i] = add_or_get(r, z)

    def coords_of(pnum_idx):
        j = pnum_idx // n_rho
        i = pnum_idx % n_rho
        theta = theta_grid[j]
        rho = rho_grid[i]
        return rho * np.sin(theta), rho * np.cos(theta)

    def triangle_has_off_axis(node_indices, atol=1e-12):
        for idx in node_indices:
            r, _ = coords_of(idx)
            if abs(r) > atol:
                return True
        return False

    def distinct(triple):
        return len(set(triple)) == 3

    triangles_p1 = []
    is_cond_p1 = []
    for j in range(N_theta):
        for i in range(n_rho - 1):
            i00 = j * n_rho + i
            i10 = j * n_rho + i + 1
            i01 = (j + 1) * n_rho + i
            i11 = (j + 1) * n_rho + i + 1
            n00, n10, n01, n11 = pnums[i00], pnums[i10], pnums[i01], pnums[i11]
            rho_mid = 0.5 * (rho_grid[i] + rho_grid[i + 1])
            is_cond = rho_mid < R_SPHERE
            if triangle_has_off_axis([i00, i10, i11]) and distinct((n00, n10, n11)):
                triangles_p1.append([n00, n10, n11])
                is_cond_p1.append(is_cond)
            if triangle_has_off_axis([i00, i11, i01]) and distinct((n00, n11, n01)):
                triangles_p1.append([n00, n11, n01])
                is_cond_p1.append(is_cond)

    nodes_p1 = np.array(nodes_p1)
    triangles_p1 = np.array(triangles_p1)
    is_cond_per_cell = np.array(is_cond_p1)

    # Build P2 nodes by adding mid-edge nodes (dedup'd via edge_to_mid).
    nodes_p2 = nodes_p1.tolist()
    edge_to_mid = {}
    triangles_p2 = []
    for tri in triangles_p1:
        v0, v1, v2 = tri

        def get_mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in edge_to_mid:
                mid = 0.5 * (nodes_p1[a] + nodes_p1[b])
                if abs(mid[0]) < 1e-12:
                    mid = np.array([0.0, mid[1]])
                edge_to_mid[key] = len(nodes_p2)
                nodes_p2.append(mid.tolist())
            return edge_to_mid[key]

        m01 = get_mid(v0, v1)
        m12 = get_mid(v1, v2)
        m20 = get_mid(v2, v0)
        triangles_p2.append([v0, v1, v2, m01, m12, m20])
    nodes_p2 = np.array(nodes_p2)
    triangles_p2 = np.array(triangles_p2)

    # Outer boundary by |position| = R_AIR.
    r = nodes_p2[:, 0]
    z = nodes_p2[:, 1]
    rho_node = np.sqrt(r * r + z * z)
    outer_idx = np.where(np.abs(rho_node - R_AIR) < 1e-9)[0]

    return nodes_p2, triangles_p2, is_cond_per_cell, outer_idx


def assemble(nodes, triangles, is_cond):
    N = nodes.shape[0]
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    skipped = 0
    for k, tri in enumerate(triangles):
        rn = nodes[tri, 0]
        zn = nodes[tri, 1]
        if np.all(rn < 1e-12):
            skipped += 1
            continue
        try:
            t = AxifemmP2Triangle(rn, zn)
        except RuntimeError:
            skipped += 1
            continue
        if t.area < 1e-20:
            skipped += 1
            continue
        Ke = t.stiffness(mu_r=MU0)
        if is_cond[k]:
            Me = t.sigma_mass(sigma=SIGMA)
        else:
            Me = np.zeros((6, 6))
        for i_l in range(6):
            gi = tri[i_l]
            for j_l in range(6):
                gj = tri[j_l]
                K[gi, gj] += Ke[i_l, j_l]
                M[gi, gj] += Me[i_l, j_l]
    return K, M, skipped


def schur_solve(K, M, fixed_idx, nodes):
    N = K.shape[0]
    axis_idx = np.where(nodes[:, 0] < 1e-12)[0]
    pinned = np.union1d(fixed_idx, axis_idx)
    free_idx = np.setdiff1d(np.arange(N), pinned)
    K_free = K[np.ix_(free_idx, free_idx)]
    M_free = M[np.ix_(free_idx, free_idx)]
    diag_M = np.diag(M_free)
    cond_local = np.where(diag_M > 1e-30 * (diag_M.max() + 1e-30))[0]
    air_local = np.where(diag_M <= 1e-30 * (diag_M.max() + 1e-30))[0]
    if len(air_local) > 0:
        K_cc = K_free[np.ix_(cond_local, cond_local)]
        K_ca = K_free[np.ix_(cond_local, air_local)]
        K_aa = K_free[np.ix_(air_local, air_local)]
        X = sla.solve(K_aa, K_ca.T, assume_a="sym")
        S = K_cc - K_ca @ X
    else:
        S = K_free[np.ix_(cond_local, cond_local)]
    M_cc = M_free[np.ix_(cond_local, cond_local)]
    S = 0.5 * (S + S.T)
    M_cc = 0.5 * (M_cc + M_cc.T)
    eigvals = sla.eigvalsh(S, M_cc, subset_by_index=[0, min(3, S.shape[0] - 1)])
    pos = eigvals[eigvals > 0]
    return pos


def main():
    print("=" * 64)
    print("Bisect 2: Python assembly on dedup'd mesh vs Python ref mesh")
    print("=" * 64)

    N_in, N_out, N_th = 4, 4, 12
    nodes, triangles, is_cond, outer_idx = build_dedup_mesh(N_in, N_out, N_th)
    print(f"  Dedup'd mesh: nodes={nodes.shape[0]} triangles={triangles.shape[0]}")
    print(f"  Outer boundary nodes: {len(outer_idx)}")
    K, M, skipped = assemble(nodes, triangles, is_cond)
    print(f"  Skipped triangles (degenerate): {skipped}")

    pos = schur_solve(K, M, outer_idx, nodes)
    tau = 1.0 / pos[0]
    MU0_local = 4 * np.pi * 1e-7
    stoll = MU0_local * SIGMA * R_SPHERE * R_SPHERE / (np.pi * np.pi)
    print(f"\n  tau_1 = {tau*1e6:.4f} us  (Stoll {stoll*1e6:.4f}, "
          f"gap = {(tau-stoll)/stoll*100:+.3f}%)")


if __name__ == "__main__":
    main()
