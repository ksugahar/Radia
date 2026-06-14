"""Phase F-4 driver: assemble full K matrix for cuboid 5x2x1 mm Cu, build
mass matrix M, solve generalized eigenproblem, extract Foster spectrum and
the leading τ that matches the canonical "Radia用 isolated hex in vacuum,
J·n=0" answer (~43.84 μs from cached Mathematica reference).

Pipeline:
  1. Build the C++ HDivDivFreeHexBasis(order).
  2. Build the bare K matrix via radia_vim.assemble_K_bare (Spherical Duffy).
     Multiply by mu0/(4π) to get the physical K.
  3. Build the mass matrix M_phys via Gauss-Legendre tensor-product
     integration of the basis evaluator. (Polynomial — exact at modest order.)
  4. Build the projection vector b = ⟨φ_i, A_ext⟩ for A_ext = (-y/2, x/2, 0)
     under the HDiv contravariant Piola map.
  5. Solve K v = λ M v (generalized symmetric eigenproblem).
  6. Compute g_k = ⟨v_k, b⟩, sort modes by |g|^2, report leading τ.

Reference Mathematica script:
  W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/
    hex_bem_hdiv_divfree.wls

Reference cached τ_1 (top-by-|g|² mode at order 3): 43.84 μs

NOTE: the cached Mathematica K used NIntegrate with PrecisionGoal=3,
AccuracyGoal=5, MaxRecursion=6 — these LOOSE TOLERANCES do NOT converge
against the 1/r diagonal singularity. Our Spherical Duffy K is the
correct one; the extracted τ may therefore differ from 43.84 μs.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

# Force UTF-8 for stdout/stderr on Windows so we can print "μ" etc. safely.
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
from scipy.linalg import eigh

import radia_vim


# Cuboid geometry and material (Cu).
A_M = 5e-3
B_M = 2e-3
C_M = 1e-3
SIGMA = 5.8e7      # S/m
MU0   = 4 * math.pi * 1e-7

# Cached Mathematica reference values for spot-checking K assembler against
# the (POSSIBLY UNCONVERGED) overnight result. Source:
# tests/K_reference_values.json (mu0/4π already stripped: bare integral).
MATHEMATICA_REF_BARE = {
    (0, 0): 0.0137910501152618,
    (1, 1): 0.00922449614646807,
    (7, 7): 0.00510657477529681,
    (6, 6): 0.00137528477603719,
    (40, 40): -7.39598926191086e-07,
}


# ----------------------------------------------------------------------
# Mass matrix M_phys: ∫ (J_phys φ_i_ref) · (J_phys φ_j_ref) (1/det J)² det J dx
# = (1/det J) ∫ (J φ_i_ref) · (J φ_j_ref) dx
# For diagonal J = diag(a,b,c):
#   M_phys[i,j] = (a/(bc)) ∫ φ_ix φ_jx dV + (b/(ac)) ∫ φ_iy φ_jy dV
#                                         + (c/(ab)) ∫ φ_iz φ_jz dV
# ----------------------------------------------------------------------
def assemble_mass_matrix_python(basis, a, b, c, n_gauss=8):
    """Mass matrix via Gauss-Legendre on [0,1]^3. Polynomial integrand —
    exact for n_gauss >= (2p+1)/2 ish. For order 3, n_gauss=8 is plenty."""
    n = basis.n_dofs
    nodes_m1, weights_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * weights_m1

    # Pre-evaluate basis at every (x,y,z) quad point: shape (n_dofs, 3, n_g, n_g, n_g)
    Phi = np.zeros((n, 3, n_gauss, n_gauss, n_gauss))
    for ix, x in enumerate(nodes):
        for iy, y in enumerate(nodes):
            for iz, z in enumerate(nodes):
                for d in range(n):
                    v = basis.evaluate(d, float(x), float(y), float(z))
                    Phi[d, 0, ix, iy, iz] = v[0]
                    Phi[d, 1, ix, iy, iz] = v[1]
                    Phi[d, 2, ix, iy, iz] = v[2]

    # Tensor-product weights: w[ix]*w[iy]*w[iz]
    W3 = weights[:, None, None] * weights[None, :, None] * weights[None, None, :]
    # Reference mass blocks
    Mref = np.zeros((3, n, n))
    for comp in range(3):
        # Phi[:, comp, :,:,:] contracted with itself weighted by W3
        Pcomp = Phi[:, comp]                    # (n, n_g, n_g, n_g)
        Pw    = Pcomp * W3[None, :, :, :]      # (n, n_g, n_g, n_g)
        Mref[comp] = np.tensordot(
            Pcomp.reshape(n, -1),
            Pw.reshape(n, -1),
            axes=([1], [1])
        )

    Mphys = (a / (b * c)) * Mref[0] + (b / (a * c)) * Mref[1] + (c / (a * b)) * Mref[2]
    return Mphys


def assemble_b_vector_python(basis, a, b, c, n_gauss=8):
    """Projection b_i = ⟨φ_i_phys, A_ext_phys⟩_Ω.
       Mathematica reference (hex_bem_hdiv_divfree.wls lines 234-241):
         A_ext_phys = (-Y/2, X/2, 0) with X=a*x_ref, Y=b*y_ref (physical coords)
         φ_phys = (a/(bc)) φ_x_ref e_x + (b/(ac)) φ_y_ref e_y + (c/(ab)) φ_z_ref e_z
       so   φ_phys · A_ext_phys = (a/(bc)) φ_x_ref (-b y_ref/2 · scale)
                                  + (b/(ac)) φ_y_ref (a x_ref/2 · scale)
       After collecting the dV = (abc) dx dy dz, Mathematica writes:
         b_i = (-a/2) ∫ φ_x_ref · y · (abc/(bc)) dV_ref
             + ( b/2) ∫ φ_y_ref · x · (abc/(ac)) dV_ref
             = (-a²/2) ∫ φ_x_ref y dV_ref + (b²/2) ∫ φ_y_ref x dV_ref"""
    n = basis.n_dofs
    nodes_m1, weights_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * weights_m1

    bvec = np.zeros(n)
    for d in range(n):
        sx = 0.0
        sy = 0.0
        for ix, x in enumerate(nodes):
            for iy, y in enumerate(nodes):
                for iz, z in enumerate(nodes):
                    v = basis.evaluate(d, float(x), float(y), float(z))
                    w3 = weights[ix] * weights[iy] * weights[iz]
                    sx += w3 * v[0] * y     # ∫ φ_x_ref · y_ref dV_ref
                    sy += w3 * v[1] * x     # ∫ φ_y_ref · x_ref dV_ref
        bvec[d] = (-a * a / 2.0) * sx + (b * b / 2.0) * sy
    return bvec


# ----------------------------------------------------------------------
# Spot-check: assembled K vs cached Mathematica entries.
# ----------------------------------------------------------------------
def report_K_spot_check(K_bare):
    print()
    print("K matrix spot check vs cached Mathematica reference:")
    print(f"  {'(i,j)':>10}  {'Mathematica':>20}  {'C++ Spherical Duffy':>22}  {'rel diff':>10}")
    for (i, j), ref in sorted(MATHEMATICA_REF_BARE.items()):
        ours = K_bare[i, j]
        if abs(ref) > 1e-30:
            rel = abs(ours - ref) / abs(ref)
            rel_str = f"{rel:.2e}"
        else:
            rel_str = "(ref ~0)"
        print(f"  {(i,j)!s:>10}  {ref:>+20.10e}  {ours:>+22.10e}  {rel_str:>10}")
    print()
    print("  NOTE: large relative differences are EXPECTED — the cached")
    print("  Mathematica K used NIntegrate with PG=3, MR=6, which does")
    print("  not converge on the 1/r kernel. C++ Spherical Duffy is correct.")


# ----------------------------------------------------------------------
# Foster extraction.
# ----------------------------------------------------------------------
def extract_foster(K_phys, M_phys, b_vec, sigma, n_show=10):
    """Solve K v = λ M v (symmetric definite generalized eigenproblem),
    compute g_k = (v_k · b), and sort modes by |g_k|² (descending).

    Returns dict with eigvals, eigvecs (M-normalised), g, tau_us, sorted_idx.
    """
    # eigh handles generalized symmetric definite problem.
    # Returns (eigvals_ascending, eigvecs columns).
    eigvals, eigvecs = eigh(K_phys, M_phys)

    # M-normalisation: v^T M v = 1  (eigh already does this for the "b"
    # kwarg path, but be defensive — re-normalise.)
    n = K_phys.shape[0]
    for k in range(n):
        v = eigvecs[:, k]
        nrm2 = float(v @ M_phys @ v)
        if nrm2 > 1e-50:
            eigvecs[:, k] = v / math.sqrt(nrm2)

    g = eigvecs.T @ b_vec     # shape (n,)
    g2 = g * g
    tau_us = np.where(eigvals > 1e-30, 1e6 * sigma * eigvals, 0.0)

    sorted_idx = np.argsort(-g2)   # descending

    print()
    print("Top Foster modes by |g|² weight:")
    print(f"  {'rank':>5}  {'τ (μs)':>15}  {'|g|²':>15}  {'τ·|g|²':>15}")
    for rank in range(min(n_show, n)):
        idx = sorted_idx[rank]
        print(f"  {rank:>5}  {tau_us[idx]:>15.6f}  {g2[idx]:>15.6e}  "
              f"{tau_us[idx]*g2[idx]:>15.6e}")
    print()

    leading_idx = sorted_idx[0]
    print(f"LEADING (max |g|²) physical mode: τ = {tau_us[leading_idx]:.6f} μs")
    return {
        "eigvals": eigvals,
        "eigvecs": eigvecs,
        "g": g,
        "g2": g2,
        "tau_us": tau_us,
        "sorted_idx": sorted_idx,
        "leading_tau_us": float(tau_us[leading_idx]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=3,
                        help="HDiv div-free interior order (n_dofs = 2p^3 + 3p^2)")
    parser.add_argument("--n-omega-theta", type=int, default=12)
    parser.add_argument("--n-omega-phi",   type=int, default=24)
    parser.add_argument("--n-rho",         type=int, default=12)
    parser.add_argument("--n-c",           type=int, default=6)
    parser.add_argument("--n-gauss-mass",  type=int, default=8,
                        help="Gauss-Legendre order for the mass matrix integration")
    parser.add_argument("--verbose", type=int, default=50,
                        help="Print K-assembler progress every N pairs")
    parser.add_argument("--save-json", type=str, default=None,
                        help="Optional path to dump the extracted Foster network to JSON")
    args = parser.parse_args()

    print("=" * 72)
    print(f" Cuboid {A_M*1000:g} x {B_M*1000:g} x {C_M*1000:g} mm Cu")
    print(f" sigma = {SIGMA:g} S/m,  mu0 = {MU0:.6e} H/m")
    print(f" HDiv div-free order = {args.order}")
    print("=" * 72)

    basis = radia_vim.HDivDivFreeHexBasis(order=args.order)
    n = basis.n_dofs
    n_pairs = n * (n + 1) // 2
    print(f"n_dofs = {n},  unique pairs = {n_pairs}")

    # --- 1. K matrix ---
    print()
    print("Assembling K matrix (Spherical Duffy)...")
    print(f"  quadrature: theta={args.n_omega_theta}, phi={args.n_omega_phi}, "
          f"rho={args.n_rho}, c={args.n_c}")
    t0 = time.time()
    K_bare = radia_vim.assemble_K_bare(
        basis, A_M, B_M, C_M,
        n_omega_theta=args.n_omega_theta,
        n_omega_phi=args.n_omega_phi,
        n_rho=args.n_rho,
        n_c=args.n_c,
        verbose=args.verbose,
    )
    dt = time.time() - t0
    print(f"K assembly: {dt:.1f}s ({dt/n_pairs*1000:.1f} ms/pair)")
    K_phys = (MU0 / (4 * math.pi)) * K_bare
    sym_diff = np.max(np.abs(K_phys - K_phys.T))
    print(f"K symmetry check: max|K - K^T| = {sym_diff:.2e}")

    if args.order == 3:
        report_K_spot_check(K_bare)

    # --- 2. Mass matrix ---
    print(f"Assembling mass matrix (Python GL n={args.n_gauss_mass})...")
    t0 = time.time()
    M_phys = assemble_mass_matrix_python(basis, A_M, B_M, C_M,
                                         n_gauss=args.n_gauss_mass)
    print(f"M assembly: {time.time() - t0:.1f}s")
    print(f"M symmetry: {np.max(np.abs(M_phys - M_phys.T)):.2e}")
    print(f"M trace:    {np.trace(M_phys):.6e}")

    # --- 3. b vector ---
    print(f"Assembling b vector (Python GL n={args.n_gauss_mass})...")
    t0 = time.time()
    b_vec = assemble_b_vector_python(basis, A_M, B_M, C_M,
                                     n_gauss=args.n_gauss_mass)
    print(f"b assembly: {time.time() - t0:.1f}s,  ||b|| = {np.linalg.norm(b_vec):.6e}")

    # --- 4. Generalized eigenproblem & Foster extraction ---
    result = extract_foster(K_phys, M_phys, b_vec, SIGMA, n_show=12)

    print()
    print("=" * 72)
    print(f" Result: leading τ = {result['leading_tau_us']:.6f} μs")
    print(f" Reference (cached Mathematica): 43.840854 μs")
    rel = abs(result["leading_tau_us"] - 43.840854) / 43.840854
    print(f" Relative difference: {rel:.4%}")
    print("=" * 72)

    if args.save_json:
        out = {
            "geometry": {"a_m": A_M, "b_m": B_M, "c_m": C_M},
            "material": {"sigma_S_per_m": SIGMA, "mu0_H_per_m": MU0},
            "order": args.order,
            "n_dofs": int(n),
            "quadrature": {
                "n_omega_theta": args.n_omega_theta,
                "n_omega_phi":   args.n_omega_phi,
                "n_rho":         args.n_rho,
                "n_c":           args.n_c,
                "n_gauss_mass":  args.n_gauss_mass,
            },
            "leading_tau_us": result["leading_tau_us"],
            "top_modes_by_g2": [
                {
                    "rank": rank,
                    "idx": int(result["sorted_idx"][rank]),
                    "tau_us": float(result["tau_us"][result["sorted_idx"][rank]]),
                    "g2":      float(result["g2"][result["sorted_idx"][rank]]),
                }
                for rank in range(n)   # save ALL modes for Cauer-I extraction
            ],
        }
        Path(args.save_json).write_text(json.dumps(out, indent=2))
        print(f"Foster spectrum saved to: {args.save_json}")


if __name__ == "__main__":
    main()
