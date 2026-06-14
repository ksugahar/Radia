"""Phase F-4 driver: assemble K matrix for A1 = 17.72 x 17.72 x 2 mm Cu
square prism (V = 628 mm^3, matched to the working cylinder R=10mm t=2mm
case for the 2026-05-08 cuboid-CLN bug isolation experiment).

Identical to extract_tau_cuboid_521.py except A_M, B_M, C_M.
Reference τ_pair[0]:
  - Cylinder (V=628, BEM-Foster): 219.32 μs
  - Cylinder (3D HCurl FEM v5): 218.67 μs
  - A1 (3D HCurl FEM v5): 65.29 μs (gauge collapse, broken)
  - A1 (this script, BEM-Foster): TBD

If A1 BEM-Foster gives ~200-220 μs (close to cylinder, modest corner
effect), then BEM is consistent across geometries → FEM v5 has a
rectangular-corner bug.

If A1 BEM-Foster gives ~60 μs (close to FEM), then BEM and FEM both
collapse on rectangular cross-sections → the user's hypothesis "BEM
might be wrong" is supported.

If A1 BEM-Foster gives an intermediate value, more investigation needed.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
from scipy.linalg import eigh

import radia_vim


# A1 geometry (V=628 mm^3, square cross-section, t=2mm, matched to cylinder)
A_M = 17.72e-3
B_M = 17.72e-3
C_M = 2e-3
SIGMA = 5.8e7
MU0   = 4 * math.pi * 1e-7


def assemble_mass_matrix_python(basis, a, b, c, n_gauss=8):
    n = basis.n_dofs
    nodes_m1, weights_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * weights_m1

    Phi = np.zeros((n, 3, n_gauss, n_gauss, n_gauss))
    for ix, x in enumerate(nodes):
        for iy, y in enumerate(nodes):
            for iz, z in enumerate(nodes):
                for d in range(n):
                    v = basis.evaluate(d, float(x), float(y), float(z))
                    Phi[d, 0, ix, iy, iz] = v[0]
                    Phi[d, 1, ix, iy, iz] = v[1]
                    Phi[d, 2, ix, iy, iz] = v[2]

    W3 = weights[:, None, None] * weights[None, :, None] * weights[None, None, :]
    Mref = np.zeros((3, n, n))
    for comp in range(3):
        Pcomp = Phi[:, comp]
        Pw    = Pcomp * W3[None, :, :, :]
        Mref[comp] = np.tensordot(
            Pcomp.reshape(n, -1),
            Pw.reshape(n, -1),
            axes=([1], [1])
        )

    Mphys = (a / (b * c)) * Mref[0] + (b / (a * c)) * Mref[1] + (c / (a * b)) * Mref[2]
    return Mphys


def assemble_b_vector_python(basis, a, b, c, n_gauss=8):
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
                    sx += w3 * v[0] * y
                    sy += w3 * v[1] * x
        bvec[d] = (-a * a / 2.0) * sx + (b * b / 2.0) * sy
    return bvec


def extract_foster(K_phys, M_phys, b_vec, sigma, n_show=10):
    eigvals, eigvecs = eigh(K_phys, M_phys)
    n = K_phys.shape[0]
    for k in range(n):
        v = eigvecs[:, k]
        nrm2 = float(v @ M_phys @ v)
        if nrm2 > 1e-50:
            eigvecs[:, k] = v / math.sqrt(nrm2)

    g = eigvecs.T @ b_vec
    g2 = g * g
    tau_us = np.where(eigvals > 1e-30, 1e6 * sigma * eigvals, 0.0)
    sorted_idx = np.argsort(-g2)

    print()
    print("Top Foster modes by |g|^2 weight:")
    print(f"  {'rank':>5}  {'tau (us)':>15}  {'|g|^2':>15}  {'tau*|g|^2':>15}")
    for rank in range(min(n_show, n)):
        idx = sorted_idx[rank]
        print(f"  {rank:>5}  {tau_us[idx]:>15.6f}  {g2[idx]:>15.6e}  "
              f"{tau_us[idx]*g2[idx]:>15.6e}")
    print()

    leading_idx = sorted_idx[0]
    print(f"LEADING (max |g|^2) physical mode: tau = {tau_us[leading_idx]:.6f} us")
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
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument("--n-omega-theta", type=int, default=12)
    parser.add_argument("--n-omega-phi",   type=int, default=24)
    parser.add_argument("--n-rho",         type=int, default=12)
    parser.add_argument("--n-c",           type=int, default=6)
    parser.add_argument("--n-gauss-mass",  type=int, default=8)
    parser.add_argument("--verbose", type=int, default=50)
    parser.add_argument("--save-json", type=str, default=None)
    args = parser.parse_args()

    print("=" * 72)
    print(f" A1 square prism {A_M*1000:g} x {B_M*1000:g} x {C_M*1000:g} mm Cu")
    print(f" V = {A_M*B_M*C_M*1e9:.3f} mm^3 (matched to cylinder R=10mm t=2mm: 628.319 mm^3)")
    print(f" sigma = {SIGMA:g} S/m,  mu0 = {MU0:.6e} H/m")
    print(f" HDiv div-free order = {args.order}")
    print("=" * 72)

    basis = radia_vim.HDivDivFreeHexBasis(order=args.order)
    n = basis.n_dofs
    n_pairs = n * (n + 1) // 2
    print(f"n_dofs = {n},  unique pairs = {n_pairs}")

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

    print(f"Assembling mass matrix (Python GL n={args.n_gauss_mass})...")
    t0 = time.time()
    M_phys = assemble_mass_matrix_python(basis, A_M, B_M, C_M,
                                         n_gauss=args.n_gauss_mass)
    print(f"M assembly: {time.time() - t0:.1f}s")

    print(f"Assembling b vector (Python GL n={args.n_gauss_mass})...")
    t0 = time.time()
    b_vec = assemble_b_vector_python(basis, A_M, B_M, C_M,
                                     n_gauss=args.n_gauss_mass)
    print(f"b assembly: {time.time() - t0:.1f}s,  ||b|| = {np.linalg.norm(b_vec):.6e}")

    result = extract_foster(K_phys, M_phys, b_vec, SIGMA, n_show=12)

    print()
    print("=" * 72)
    print(f" A1 leading tau = {result['leading_tau_us']:.6f} us")
    print(f" Comparison:")
    print(f"   Cylinder (V=628, BEM-Foster): 219.32 us")
    print(f"   Cylinder (FEM v5):            218.67 us")
    print(f"   A1 (FEM v5, broken):           65.29 us")
    print(f"   A1 (this BEM):                {result['leading_tau_us']:>9.4f} us")
    print("=" * 72)

    if args.save_json:
        out = {
            "geometry": {"a_m": A_M, "b_m": B_M, "c_m": C_M,
                          "label": "A1 square prism, V=628 mm^3"},
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
                for rank in range(min(20, n))
            ],
        }
        Path(args.save_json).write_text(json.dumps(out, indent=2))
        print(f"Foster spectrum saved to: {args.save_json}")


if __name__ == "__main__":
    main()
