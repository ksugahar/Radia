"""Bisect: feed Python's K, M into C++ test's Schur reduction.

If the Schur logic from test_p2_cpp_sphere.run_one is correct, applying it
to Python reference K, M should give the same answer Python itself produces.
If it gives 50% off, the bug is in the Schur logic. If it agrees with Python,
the bug is in C++ assembly.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.linalg as sla

sys.path.insert(0, ".")
from test_p2_sphere_stoll import (
    generate_sphere_mesh, assemble_p2_split, boundary_outer_nodes,
    conductor_dofs, solve_schur_eigenvalue,
    SIGMA, R_AIR, R_SPHERE, EPS_AXIS, stoll_tau_n, MU0,
)


def my_schur_eigenvalue(K, M, fixed_idx, nodes):
    """Match the Schur logic from test_p2_cpp_sphere.run_one.

    free := all - fixed_idx (NGSolve dirichlet="outer" pins outer DOFs)
    cond := free DOFs with M_diag > 0
    air  := free DOFs with M_diag == 0 and K_diag > 0
    axis-zero-drop := M_diag == 0 and K_diag == 0
    S = K_cc - K_ca @ inv(K_aa) @ K_ac
    """
    N = K.shape[0]
    # In my C++ test, the only "pinned" DOFs are outer-boundary (fixed_idx).
    # Axis DOFs are auto-decoupled because their K_diag = 0 via the (r_i)
    # factor in N_A_i = r_i psi_i / r, so they end up in "axis-zero-drop"
    # below, not actively pinned.
    free_idx = np.setdiff1d(np.arange(N), fixed_idx)

    K_free = K[np.ix_(free_idx, free_idx)]
    M_free = M[np.ix_(free_idx, free_idx)]
    K_free = 0.5 * (K_free + K_free.T)
    M_free = 0.5 * (M_free + M_free.T)

    diag_M = np.diag(M_free)
    diag_K = np.diag(K_free)
    cond_mask = diag_M > 1e-30 * (np.max(diag_M) + 1e-30)
    nonzero_K = diag_K > 1e-30 * (np.max(diag_K) + 1e-30)
    air_mask = (~cond_mask) & nonzero_K
    cond_idx = np.where(cond_mask)[0]
    air_idx = np.where(air_mask)[0]
    print(f"  Schur split: cond={len(cond_idx)}  air={len(air_idx)}  "
          f"axis-zero={(~nonzero_K).sum()}")

    if len(air_idx) > 0:
        K_cc = K_free[np.ix_(cond_idx, cond_idx)]
        K_ca = K_free[np.ix_(cond_idx, air_idx)]
        K_aa = K_free[np.ix_(air_idx, air_idx)]
        X = sla.solve(K_aa, K_ca.T, assume_a="sym")
        S = K_cc - K_ca @ X
    else:
        S = K_free[np.ix_(cond_idx, cond_idx)]
    M_cc = M_free[np.ix_(cond_idx, cond_idx)]
    S = 0.5 * (S + S.T)
    M_cc = 0.5 * (M_cc + M_cc.T)
    eigvals = sla.eigvalsh(S, M_cc, subset_by_index=[0, min(3, S.shape[0]-1)])
    pos = eigvals[eigvals > 0]
    return pos


def main():
    print("=" * 64)
    print("Bisect: feed Python K, M into C++ Schur logic")
    print("=" * 64)
    tau_stoll = stoll_tau_n(1)
    print(f"Stoll tau_1 = {tau_stoll*1e6:.4f} us")

    N_in, N_out, N_th = 4, 4, 12
    nodes, triangles, is_cond_cell = generate_sphere_mesh(N_in, N_out, N_th)
    bd_idx = boundary_outer_nodes(nodes)
    cond_idx_py = conductor_dofs(nodes, triangles, is_cond_cell)
    K, M = assemble_p2_split(nodes, triangles, is_cond_cell, SIGMA)
    print(f"  Python (N_in, N_out, N_th)=({N_in}, {N_out}, {N_th}): "
          f"nodes={nodes.shape[0]}, triangles={triangles.shape[0]}")

    # Python's own Schur (the reference)
    print(f"\nPython reference Schur:")
    pos_ref, n_dropped = solve_schur_eigenvalue(K, M, bd_idx, cond_idx_py,
                                                  nodes, n_eigs=3)
    tau_ref = 1.0 / pos_ref[0]
    print(f"  tau_1 = {tau_ref*1e6:.4f} us  gap = {(tau_ref - tau_stoll)/tau_stoll*100:+.3f}%")

    # My C++ test's Schur logic, applied to Python K, M.
    print(f"\nC++ test Schur logic on Python K, M:")
    pos_my = my_schur_eigenvalue(K, M, bd_idx, nodes)
    if len(pos_my) > 0:
        tau_my = 1.0 / pos_my[0]
        print(f"  tau_1 = {tau_my*1e6:.4f} us  gap = {(tau_my - tau_stoll)/tau_stoll*100:+.3f}%")
    else:
        print(f"  No positive eigenvalues")

    print(f"\n=== Conclusion ===")
    if len(pos_my) > 0:
        ratio = tau_my / tau_ref
        if abs(ratio - 1) < 0.05:
            print(f"  My Schur AGREES with Python ref. So C++ assembly is wrong.")
        else:
            print(f"  My Schur DISAGREES with Python ref (ratio={ratio:.3f}). "
                  f"Schur logic is wrong.")


if __name__ == "__main__":
    main()
