"""extract_tau_A1_c4v.py — A1 cuboid Foster + Cauer using C4v block reduction.

Wraps extract_tau_A1_square_mp's K, M, b assembly with the c4v_symmetry
projection onto an irrep block (default A2g, the one excited by B∥z).  All
heavy mpmath operations (K assembly, eigenvalue) are then performed on the
much smaller A2g sub-block:

  p=3  : 81 dof  →  7 dof    (1547× faster eigen)
  p=5  : 325 dof →  24 dof   (2500× faster eigen)
  p=7  : 833 dof →  58 dof   (2970× faster eigen)

K assembly itself is still done in full at the moment (a future direct-
projection assembler would shrink that too).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
from mpmath import mp, mpf

# Reuse VIM building blocks.
from extract_tau_A1_square_mp import (
    HDivDivFreeHexBasisMP,
    assemble_mass_matrix_mp,
    assemble_b_vectors_mp,
    combine_b_vector,
    assemble_K_bare_mp,
    assemble_K_bare_mp_hoisted,
    extract_foster_mp,
    cauer_from_spectrum_mp,
    print_cauer_rungs,
    diagnose_cauer_rungs,
    print_asymptote_fit,
    save_foster_spectrum_json,
    A_M, B_M, C_M, SIGMA, MU0,
)

import c4v_symmetry as c4v


def numpy_to_mp(U: np.ndarray):
    """Convert a numpy float matrix to an mpmath matrix (column-major)."""
    rows, cols = U.shape
    M = mp.matrix(rows, cols)
    for i in range(rows):
        for j in range(cols):
            M[i, j] = mpf(repr(float(U[i, j])))
    return M


def project_KMb(K, M, b_vec, U_mp):
    """Compute K_block = U^T K U, M_block = U^T M U, b_block = U^T b in mpmath."""
    Ut = U_mp.T
    K_block = Ut * K * U_mp
    M_block = Ut * M * U_mp
    b_block = Ut * b_vec
    return K_block, M_block, b_block


def run_c4v(*, order: int, irrep: str = "A2g",
            n_t: int = 8, n_p: int = 16, n_r: int = 8, n_c: int = 4,
            n_gauss_mass: int = 8,
            B_dir=(0.0, 0.0, 1.0),
            n_processes: int = 1,
            n_cauer_stages: int = 6,
            use_hoisted: bool = True,
            tau_min_us: float = 0.0,
            save_foster_path: str = None):
    print("=" * 72)
    print(f" A1 cuboid C4v-reduced Foster + Cauer extraction")
    print(f"   order={order}, irrep={irrep}, mpmath dps={mp.dps}")
    print(f"   B={B_dir}, hoisted={use_hoisted}, n_processes={n_processes}")
    print("=" * 72)

    # Selection rule check
    excited = c4v.b_irrep_assignment(B_dir)
    if irrep not in excited:
        print(f"  WARNING: B={B_dir} does not couple to irrep {irrep};"
              f" expected: {excited}")

    basis = HDivDivFreeHexBasisMP(order=order)
    n = basis.n_dofs
    print(f"  full basis: n_dofs = {n}")

    # Phase 1-2: build C4v irrep basis (numpy, FP64)
    print(f"  building C4v irrep basis ({irrep})...")
    t0 = time.time()
    U_np = c4v.get_irrep_basis(order, irrep)
    n_block = U_np.shape[1]
    print(f"    {irrep} block: {n_block} dof  "
          f"(eigen speedup ≈ {(n/n_block)**3:.0f}×)")
    print(f"    C4v projection setup: {time.time()-t0:.2f} s")

    # Convert to mpmath
    print(f"  converting U_{irrep} to mpmath...")
    t0 = time.time()
    U_mp = numpy_to_mp(U_np)
    print(f"    U_mp shape: ({n}, {n_block}), {time.time()-t0:.2f} s")

    # Assemble M (full, then project)
    print("Assembling mass matrix M (full)...")
    t0 = time.time()
    M = assemble_mass_matrix_mp(basis, A_M, B_M, C_M, n_gauss=n_gauss_mass)
    print(f"  M done in {time.time()-t0:.1f}s")

    # Assemble b-vectors (full, then combine and project)
    print("Assembling b-vectors (3 directions)...")
    t0 = time.time()
    b_x, b_y, b_z = assemble_b_vectors_mp(basis, A_M, B_M, C_M,
                                          n_gauss=n_gauss_mass)
    print(f"  b_x,b_y,b_z done in {time.time()-t0:.1f}s")
    b_vec = combine_b_vector(b_x, b_y, b_z, B_dir)

    # Assemble K (full, costly step)
    print(f"Assembling K bare ({'hoisted' if use_hoisted else 'per-entry'}, "
          f"spherical Duffy)...")
    t0 = time.time()
    if use_hoisted:
        K = assemble_K_bare_mp_hoisted(
            basis, A_M, B_M, C_M,
            n_omega_theta=n_t, n_omega_phi=n_p, n_rho=n_r, n_c=n_c,
            n_processes=n_processes,
        )
    else:
        K = assemble_K_bare_mp(
            basis, A_M, B_M, C_M,
            n_omega_theta=n_t, n_omega_phi=n_p, n_rho=n_r, n_c=n_c,
            n_processes=n_processes,
        )
    print(f"  K done in {time.time()-t0:.1f}s")

    prefactor = MU0 / (4 * mp.pi)
    for i in range(n):
        for j in range(n):
            K[i, j] = prefactor * K[i, j]

    # Phase 3: project K, M, b into A2g block
    print(f"Projecting K, M, b into irrep {irrep} block...")
    t0 = time.time()
    K_block, M_block, b_block = project_KMb(K, M, b_vec, U_mp)
    print(f"  K_{irrep}: ({n_block}, {n_block}), {time.time()-t0:.1f}s")

    # Diagnostic: |b_block| should be a good fraction of |b_vec|
    b_full_norm = mp.sqrt(sum(b_vec[i, 0]**2 for i in range(n)))
    b_block_norm = mp.sqrt(sum(b_block[i, 0]**2 for i in range(n_block)))
    ratio = float(b_block_norm / b_full_norm) if b_full_norm > 0 else 0
    print(f"  |b_block| / |b_full| = {ratio:.6f} "
          f"(expect ~1 for fully-selected irrep, ~0 otherwise)")

    # Foster decomposition on the block
    print(f"\nFoster decomposition on {irrep} block (mp.eig generalized)...")
    t0 = time.time()
    result = extract_foster_mp(K_block, M_block, b_block, SIGMA,
                               n_show=min(10, n_block))
    print(f"  Foster done in {time.time()-t0:.1f}s")

    if save_foster_path:
        save_foster_spectrum_json(
            save_foster_path,
            order=order, B_dir=B_dir,
            tau_us_mp=result["tau_us"], g2_mp=result["g2"],
            leading_tau_us=result["leading_tau_us"],
            dps=mp.dps, n_t=n_t, n_p=n_p, n_r=n_r, n_c=n_c,
        )

    # Cauer extraction
    print("\nCauer-I extraction (Hankel-Padé QD-Wynn, mpmath)...")
    if tau_min_us > 0:
        print(f"  τ-filter: tau_min_us = {tau_min_us}")
    t0 = time.time()
    rungs, c0_rel, p_array = cauer_from_spectrum_mp(
        result["tau_us"], result["g2"],
        n_use=min(80, n_block),
        n_moments=max(20, 4 * n_cauer_stages),
        max_stages=2 * n_cauer_stages,
        tau_min_us=tau_min_us,
    )
    print(f"  Cauer done in {time.time()-t0:.1f}s")
    print_cauer_rungs(rungs)
    diagnose_cauer_rungs(rungs, c0_rel)
    print_asymptote_fit(rungs)

    print()
    print("=" * 72)
    print(f" C4v-reduced (irrep={irrep}, order={order}) summary:")
    print(f"   leading τ = {mp.nstr(result['leading_tau_us'], 12)} μs")
    if rungs:
        print(f"   Cauer rung-0 τ_pair = "
              f"{mp.nstr(rungs[0]['tau_pair_us'], 12)} μs")
    print("=" * 72)

    return result, rungs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--dps", type=int, default=30)
    ap.add_argument("--irrep", type=str, default="A2g",
                    help="D_4h irrep to project onto (A2g for B∥z)")
    ap.add_argument("--B-dir", type=str, default="0,0,1",
                    help="Comma-separated direction of applied B")
    ap.add_argument("--n-processes", type=int, default=8)
    ap.add_argument("--n-cauer-stages", type=int, default=6)
    ap.add_argument("--no-hoist", action="store_true")
    ap.add_argument("--tau-min-us", type=float, default=0.0)
    ap.add_argument("--save-foster", type=str, default=None)
    args = ap.parse_args()

    mp.dps = args.dps

    B_dir = tuple(float(x) for x in args.B_dir.split(","))
    if len(B_dir) != 3:
        raise SystemExit("--B-dir must be 'Bx,By,Bz'")

    run_c4v(
        order=args.order,
        irrep=args.irrep,
        B_dir=B_dir,
        n_processes=args.n_processes,
        n_cauer_stages=args.n_cauer_stages,
        use_hoisted=not args.no_hoist,
        tau_min_us=args.tau_min_us,
        save_foster_path=args.save_foster,
    )


if __name__ == "__main__":
    main()
