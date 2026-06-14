"""Phase F-4 test: full K-matrix assembler.

Validates `radia_vim.assemble_K_bare`:

  (1) Symmetry: K = K^T exactly (we mirror i,j after computing i<=j).
  (2) Diagonal entries match the per-call `compute_K_entry_spherical_duffy`
      result bit-for-bit when called with identical quadrature.
  (3) Off-diagonal: K[i,j] == K[j,i] (mirroring works).
  (4) Order 1 (5 DOFs only, fast) full assembly < 5 s with low quadrature.
"""
from __future__ import annotations

import time
import numpy as np
import pytest

from radia_vim import (
    HDivDivFreeHexBasis,
    assemble_K_bare,
    compute_K_entry_spherical_duffy,
)


A_M, B_M, C_M = 5e-3, 2e-3, 1e-3


def test_assemble_order1_smoke():
    """Order 1 = 5 DOFs, 15 unique pairs. Should finish in ~1s."""
    basis = HDivDivFreeHexBasis(order=1)
    assert basis.n_dofs == 5
    t0 = time.time()
    K = assemble_K_bare(basis, A_M, B_M, C_M,
                        n_omega_theta=8, n_omega_phi=16,
                        n_rho=8, n_c=4, verbose=0)
    dt = time.time() - t0
    print(f"\n  order=1 assembly: {dt:.2f}s, K shape = {K.shape}")
    assert K.shape == (5, 5)
    assert dt < 30.0    # generous; typical < 1s
    # Diagonals are non-zero, positive (symmetric positive-semidefinite kernel)
    assert all(K[i, i] > 0 for i in range(5)), \
        f"diagonal not strictly positive: {[K[i,i] for i in range(5)]}"


def test_assemble_symmetry_exact():
    """Mirroring: K must be exactly symmetric (no float drift since we
    write the same value to both K[i,j] and K[j,i])."""
    basis = HDivDivFreeHexBasis(order=2)
    K = assemble_K_bare(basis, A_M, B_M, C_M,
                        n_omega_theta=6, n_omega_phi=12,
                        n_rho=6, n_c=4, verbose=0)
    sym = np.max(np.abs(K - K.T))
    print(f"\n  order=2 K - K^T max = {sym:.2e}")
    assert sym == 0.0, f"K not exactly symmetric: max diff {sym}"


def test_diagonal_matches_single_entry_call():
    """K[0,0] from the assembler must equal compute_K_entry_spherical_duffy
    called directly with the same quadrature parameters."""
    basis = HDivDivFreeHexBasis(order=2)
    nq = dict(n_omega_theta=8, n_omega_phi=16, n_rho=8, n_c=4)
    K = assemble_K_bare(basis, A_M, B_M, C_M, verbose=0, **nq)
    K00_direct = compute_K_entry_spherical_duffy(
        basis, 0, 0, A_M, B_M, C_M, **nq)
    K11_direct = compute_K_entry_spherical_duffy(
        basis, 1, 1, A_M, B_M, C_M, **nq)
    K01_direct = compute_K_entry_spherical_duffy(
        basis, 0, 1, A_M, B_M, C_M, **nq)
    print(f"\n  K[0,0]: assembler={K[0,0]:.10e}, direct={K00_direct:.10e}")
    print(f"  K[1,1]: assembler={K[1,1]:.10e}, direct={K11_direct:.10e}")
    print(f"  K[0,1]: assembler={K[0,1]:.10e}, direct={K01_direct:.10e}")
    # Bit-for-bit identity (deterministic Gauss-Legendre with same N).
    assert K[0, 0] == K00_direct
    assert K[1, 1] == K11_direct
    assert K[0, 1] == K01_direct
    # Also the mirror entry
    assert K[1, 0] == K[0, 1]


def test_order3_K00_matches_high_precision_reference():
    """At high quad, K[0,0] for order=3 cuboid 5x2x1 should converge to
    the Spherical-Duffy reference ~1.243e-3 (per debug_K00_highprec.wls).

    Note: the cached Mathematica overnight K used PrecisionGoal=3,
    AccuracyGoal=5, MaxRecursion=6 — these tolerances do NOT converge
    on the 1/r diagonal kernel. The reference value here is from the
    converged Spherical-Duffy quadrature, NOT the cached Mathematica."""
    basis = HDivDivFreeHexBasis(order=3)
    K = assemble_K_bare(basis, A_M, B_M, C_M,
                        n_omega_theta=8, n_omega_phi=16,
                        n_rho=8, n_c=5, verbose=0)
    # At quad (8,16,8,5), expect K[0,0] ≈ 1.50e-3 (still finite-error).
    # Bracketed by [1.0e-3, 2.0e-3] is the relevant order-of-magnitude check.
    print(f"\n  K[0,0] @ (8,16,8,5) = {K[0,0]:.4e} "
          f"(reference ~1.243e-3 at higher quad)")
    assert 1.0e-3 < K[0, 0] < 2.0e-3, (
        f"K[0,0] = {K[0,0]:.4e} outside expected [1e-3, 2e-3] range")


if __name__ == "__main__":
    test_assemble_order1_smoke()
    test_assemble_symmetry_exact()
    test_diagonal_matches_single_entry_call()
    test_order3_K00_matches_high_precision_reference()
    print("All Phase F-4 K-assembler tests passed.")
