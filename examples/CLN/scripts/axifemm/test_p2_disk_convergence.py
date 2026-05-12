"""test_p2_disk_convergence.py — Phase B1b Item 2
P2 self-convergence on a Cu disk.

Setup: solid Cu disk [0, R] x [-t/2, t/2] in (r, z) plane,
  R = 10 mm, t = 2 mm, sigma = 5.8e7 S/m.

Dirichlet BC: A_phi = 0 on outer boundary (r=R, z=+/-t/2). On the axis r=0
the natural Henrotte basis already enforces regularity (no extra BC needed).

We solve the generalized eigenvalue problem K v = lambda M v restricted to
interior DOFs, and track the leading eigenvalue lambda_1 vs mesh size h.

For Henrotte axisym P2 (basis order p=2): expected h^4 convergence
(Strang-Fix theorem for eigenvalues of symmetric coercive forms with a
P2-conforming subspace; eigenvalue rate is 2p = 4).

Reference is the finest-mesh P2 eigenvalue (Richardson extrapolation in
spirit). P1 axifemm comparison is deferred to Phase B2 (C++ port) because
the existing axifemm-core P1 uses a different DOF convention (FEMM V_j
variable: V_j = A_phi(r_j, z_j) * 2*pi*r_j) and a direct apples-to-apples
comparison requires aligning the basis space.

NOTE on physical setup: with all outer DOFs pinned to 0 (Dirichlet) the
eigenmodes are interior diffusion modes; lambda_1 corresponds to the
slowest-decaying mode, which is the physical analog of the leading eddy-
current time constant. The mode is mesh-converged, not analytically known.
"""
from __future__ import annotations

import sys
sys.stdout.reconfigure(encoding="utf-8")

import time

import numpy as np
import scipy.linalg as sla

# Reuse the P1 element matrices from axifemm_core.
sys.path.insert(0, str(__file__.rsplit("\\", 1)[0] if "\\" in __file__
                       else __file__.rsplit("/", 1)[0]))

from axifemm_p2_triangle import AxifemmP2Triangle, MU0


R_DISK = 10.0e-3  # 10 mm
T_DISK = 2.0e-3   # 2 mm
SIGMA = 5.8e7
B0 = 1.0


def generate_structured_mesh(N_r, N_z):
    """Structured triangular mesh of disk cross-section [0, R] x [-t/2, t/2].

    Returns:
        nodes_p1 (N_p1, 2): vertex coords (r, z)
        triangles_p1 (M, 3): vertex indices per element
        nodes_p2 (N_p2, 2): vertex + mid-edge coords
        triangles_p2 (M, 6): full P2 nodal indices [v0, v1, v2, m01, m12, m20]
    """
    r_pts = np.linspace(0.0, R_DISK, N_r + 1)
    z_pts = np.linspace(-T_DISK / 2.0, T_DISK / 2.0, N_z + 1)
    nodes_p1 = []
    for z in z_pts:
        for r in r_pts:
            nodes_p1.append([r, z])
    nodes_p1 = np.array(nodes_p1)
    N_p1 = nodes_p1.shape[0]

    triangles_p1 = []
    for j in range(N_z):
        for i in range(N_r):
            n00 = j * (N_r + 1) + i
            n10 = j * (N_r + 1) + i + 1
            n01 = (j + 1) * (N_r + 1) + i
            n11 = (j + 1) * (N_r + 1) + i + 1
            # Split quad into 2 triangles (anti-diagonal).
            triangles_p1.append([n00, n10, n11])
            triangles_p1.append([n00, n11, n01])
    triangles_p1 = np.array(triangles_p1)

    # Build P2 nodes: copy P1, then add mid-edge nodes
    nodes_p2 = nodes_p1.tolist()
    edge_to_mid = {}
    triangles_p2 = []
    for tri in triangles_p1:
        v0, v1, v2 = tri

        def get_mid(a, b):
            key = tuple(sorted([a, b]))
            if key not in edge_to_mid:
                mid = 0.5 * (nodes_p1[a] + nodes_p1[b])
                edge_to_mid[key] = len(nodes_p2)
                nodes_p2.append(mid.tolist())
            return edge_to_mid[key]

        m01 = get_mid(v0, v1)
        m12 = get_mid(v1, v2)
        m20 = get_mid(v2, v0)
        triangles_p2.append([v0, v1, v2, m01, m12, m20])
    nodes_p2 = np.array(nodes_p2)
    triangles_p2 = np.array(triangles_p2)

    return nodes_p1, triangles_p1, nodes_p2, triangles_p2


def boundary_nodes(nodes, tol=1e-12):
    """Return indices of nodes on outer boundary (r=R, z=+/-t/2).
    Axis r=0 nodes are NOT boundary — left free for FEMM Henrotte basis.
    """
    r = nodes[:, 0]
    z = nodes[:, 1]
    on_outer_r = np.abs(r - R_DISK) < tol
    on_top_z = np.abs(z - T_DISK / 2.0) < tol
    on_bot_z = np.abs(z + T_DISK / 2.0) < tol
    return np.where(on_outer_r | on_top_z | on_bot_z)[0]


def assemble_p2(nodes, triangles, sigma):
    """Assemble K, M for P2 Henrotte basis (6 DOF / tri)."""
    N = nodes.shape[0]
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    for tri in triangles:
        rn = nodes[tri, 0]
        zn = nodes[tri, 1]
        t_elem = AxifemmP2Triangle(rn, zn)
        Ke = t_elem.stiffness(mu_r=MU0)
        Me = t_elem.sigma_mass(sigma=sigma)
        for i_l in range(6):
            for j_l in range(6):
                K[tri[i_l], tri[j_l]] += Ke[i_l, j_l]
                M[tri[i_l], tri[j_l]] += Me[i_l, j_l]
    return K, M


def solve_eigenvalue(K, M, fixed_idx):
    """Reduce system by removing fixed DOFs, then solve K v = lambda M v.

    Returns sorted (positive) eigenvalues [in seconds, since lambda has
    units that depend on the Henrotte normalization; we compare ratios.]
    """
    N = K.shape[0]
    free_idx = np.setdiff1d(np.arange(N), fixed_idx)
    K_ff = K[np.ix_(free_idx, free_idx)]
    M_ff = M[np.ix_(free_idx, free_idx)]
    # Symmetrize numerically.
    K_ff = 0.5 * (K_ff + K_ff.T)
    M_ff = 0.5 * (M_ff + M_ff.T)
    # Generalized eigenvalue (smallest first).
    eigvals = sla.eigvalsh(K_ff, M_ff, subset_by_index=[0, min(9, K_ff.shape[0] - 1)])
    pos = eigvals[eigvals > 0]
    return pos


def main():
    print("=" * 64)
    print("Phase B1b Item 2: P2 disk eigenvalue self-convergence")
    print("=" * 64)
    print(f"\nCu disk R={R_DISK*1e3} mm, t={T_DISK*1e3} mm, "
          f"sigma={SIGMA:.2e} S/m")

    # Mesh refinement levels (Nr, Nz). Keep aspect roughly constant.
    levels = [
        (5, 2),
        (10, 4),
        (20, 8),
        (40, 16),
        (60, 24),
    ]

    results = []

    for N_r, N_z in levels:
        h = R_DISK / N_r
        nodes_p1, tris_p1, nodes_p2, tris_p2 = generate_structured_mesh(N_r, N_z)
        bd_p2 = boundary_nodes(nodes_p2)

        t0 = time.time()
        K2, M2 = assemble_p2(nodes_p2, tris_p2, SIGMA)
        evals_p2 = solve_eigenvalue(K2, M2, bd_p2)
        t_p2 = time.time() - t0
        lam_p2 = evals_p2[0] if len(evals_p2) > 0 else np.nan
        # FEM Henrotte with sigma absorbed in M: K v = lambda M v gives
        # lambda = 1/tau (proved by axisymmetric diffusion weak form).
        tau_p2 = 1.0 / lam_p2 if lam_p2 > 0 else np.nan

        results.append((h, lam_p2, tau_p2, t_p2, nodes_p2.shape[0]))

        print(f"  Nr={N_r:>3d} Nz={N_z:>3d} (h={h*1e3:>5.2f} mm)  "
              f"P2: N={nodes_p2.shape[0]:>5d} tau={tau_p2*1e6:>10.5f} us "
              f"({t_p2:>5.2f} s)")

    # Convergence rate estimates: assume finest-mesh result is "exact"
    print("\nConvergence study (P2, finest mesh as reference):")
    print(f"  {'h [mm]':>8s}  {'tau [us]':>12s}  {'|err| [us]':>12s}  "
          f"{'rate':>6s}")
    ref_tau = results[-1][2]
    prev_err = None
    prev_h = None
    for (h, lam, tau, t_e, n) in results[:-1]:
        err = abs(tau - ref_tau)
        rate = (np.log(err / prev_err) / np.log(h / prev_h)
                if prev_err is not None and prev_h is not None
                   and prev_err > 0 and err > 0
                else None)
        rate_str = f"{rate:>5.2f}" if rate is not None else "  -  "
        print(f"  {h*1e3:>8.3f}  {tau*1e6:>12.6f}  "
              f"{err*1e6:>12.6e}  {rate_str}")
        prev_err = err
        prev_h = h

    print()
    print("Expected: P2 ~h^4 for SMOOTH domains; rectangular Dirichlet has")
    print("corner singularities (J_1 zeros + sin BC) that cap the rate to")
    print("roughly h^2 in practice.")

    # --- Analytical Foster tau closed form ---------------------------------
    # tau_{n, m} = mu_0 sigma / ((alpha_{1, n} / R)^2 + (m pi / t)^2)
    # alpha_{1, n} = n-th zero of Bessel J_1
    from scipy.special import jn_zeros
    j1_zeros = jn_zeros(1, 4)
    print()
    print("Analytical Foster tau for solid Cu disk with PEC outer boundary:")
    print(f"  tau_{{n,m}} = mu_0 sigma / ((alpha_{{1,n}}/R)^2 + (m pi / t)^2)")
    print()
    print(f"  {'n':>3s} {'m':>3s} {'alpha_1,n':>10s} {'k_r [/m]':>10s} "
          f"{'k_z [/m]':>10s} {'tau [us]':>10s}")
    analytical_taus = []
    for n_r in range(1, 4):
        for m_z in range(1, 4):
            kr = j1_zeros[n_r - 1] / R_DISK
            kz = m_z * np.pi / T_DISK
            tau_nm = MU0 * SIGMA / (kr ** 2 + kz ** 2)
            print(f"  {n_r:>3d} {m_z:>3d} {j1_zeros[n_r - 1]:>10.4f} "
                  f"{kr:>10.2f} {kz:>10.2f} {tau_nm*1e6:>10.4f}")
            analytical_taus.append((n_r, m_z, tau_nm))

    # Sort by tau descending (slowest decay first):
    analytical_taus.sort(key=lambda x: -x[2])
    tau_leading_analytical = analytical_taus[0][2]
    tau_leading_p2 = results[-1][2]
    print()
    print(f"  Leading mode: (n={analytical_taus[0][0]}, m={analytical_taus[0][1]})")
    print(f"  Analytical tau_{{1,1}} = {tau_leading_analytical*1e6:.4f} us")
    print(f"  P2 finest mesh tau     = {tau_leading_p2*1e6:.4f} us "
          f"({(tau_leading_p2/tau_leading_analytical - 1)*100:+.2f}%)")


if __name__ == "__main__":
    main()
