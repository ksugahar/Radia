"""extract_tau_cuboid_521_gpu.py — GPU port of extract_tau_cuboid_521.py.

Replaces radia_vim.assemble_K_bare (C++) with assemble_K_cupy (GPU)
for K matrix. M and b use the existing Python (CPU) functions from
the original script (these are small for single hex).

Provides quadrature density sweep: baseline (12,24,12,6) vs doubled
(24,48,24,12) to test if Spherical Duffy converges with refinement,
analogous to sphere self-cell subdivision.
"""
from __future__ import annotations

import argparse
import json
import math
import time
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

# Import the GPU K assembly + Python M, b from the existing scripts.
script_dir = Path(__file__).parent
legacy_dir = script_dir.parent.parent.parent / "examples" / "CLN" / "scripts" / "ngsolve_validation" / "legacy_fp64"
sys.path.insert(0, str(legacy_dir))
from hex_vim_cupy_kassembly import assemble_K_cupy

# extract_tau_cuboid_521.py's M, b assembly functions
from extract_tau_cuboid_521 import (
    assemble_mass_matrix_python,
    assemble_b_vector_python,
    extract_foster,
)
import radia_vim


A_M = 5e-3
B_M = 2e-3
C_M = 1e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7


def run_one(order, n_th, n_ph, n_rh, n_c, save_path=None, verbose=True):
    print("=" * 72, flush=True)
    print(f" Cuboid 5x2x1 mm Cu order={order} quad=({n_th},{n_ph},{n_rh},{n_c})",
          flush=True)
    print("=" * 72, flush=True)

    basis = radia_vim.HDivDivFreeHexBasis(order=order)
    n = basis.n_dofs

    # K via GPU
    print(f"Assembling K matrix (GPU Spherical Duffy)...", flush=True)
    t0 = time.time()
    K_bare = assemble_K_cupy(order, A_M, B_M, C_M,
                              n_th=n_th, n_ph=n_ph, n_rh=n_rh, n_c=n_c,
                              use_gpu=True, dtype=np.float64,
                              verbose=verbose)
    t_K = time.time() - t0
    print(f"K assembly (GPU): {t_K:.2f}s", flush=True)
    K_phys = (MU0 / (4 * math.pi)) * K_bare
    sym_diff = np.max(np.abs(K_phys - K_phys.T))
    print(f"K symmetry: max|K - K^T| = {sym_diff:.2e}", flush=True)

    # Mass matrix (CPU)
    print(f"Assembling M (CPU, n_gauss=8)...", flush=True)
    t0 = time.time()
    M_phys = assemble_mass_matrix_python(basis, A_M, B_M, C_M, n_gauss=8)
    t_M = time.time() - t0

    # b vector (CPU)
    print(f"Assembling b (CPU)...", flush=True)
    t0 = time.time()
    b_vec = assemble_b_vector_python(basis, A_M, B_M, C_M, n_gauss=8)
    t_b = time.time() - t0

    # Foster + leading tau
    result = extract_foster(K_phys, M_phys, b_vec, SIGMA, n_show=12)

    print(flush=True)
    print(f" Result: leading τ = {result['leading_tau_us']:.6f} μs",
          flush=True)

    out = {
        "order": order, "n_dofs": int(n),
        "quadrature": {"n_omega_theta": n_th, "n_omega_phi": n_ph,
                       "n_rho": n_rh, "n_c": n_c},
        "t_K_assembly_s": float(t_K),
        "t_M_assembly_s": float(t_M),
        "t_b_assembly_s": float(t_b),
        "leading_tau_us": float(result["leading_tau_us"]),
        "all_tau_us": [float(t) for t in result["tau_us"]],
        "all_g2": [float(g2) for g2 in result["g2"]],
    }
    if save_path:
        Path(save_path).write_text(json.dumps(out, indent=2))
        print(f"Saved: {save_path}", flush=True)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument("--sweep", action="store_true",
                        help="Run baseline vs doubled quadrature comparison")
    args = parser.parse_args()

    if args.sweep:
        # Run baseline + intermediate for comparison
        # Doubled (24,48,24,12) needs Phi tensor 93 GB - OOM on local + A100.
        cases = [
            ("baseline (12,24,12,6)", 12, 24, 12, 6),
            ("intermed (16,32,16,8)", 16, 32, 16, 8),     # Phi ~ 8 GB
            ("refined  (20,40,20,10)", 20, 40, 20, 10),   # Phi ~ 16 GB local lim
        ]
    else:
        cases = [
            (f"default order={args.order}", 12, 24, 12, 6),
        ]

    results = []
    for label, n_th, n_ph, n_rh, n_c in cases:
        print(f"\n### {label} ###", flush=True)
        r = run_one(args.order, n_th, n_ph, n_rh, n_c,
                    save_path=f"C:/Users/ADMINI~1/AppData/Local/Temp/1/cuboid_p{args.order}_"
                              f"q{n_th}_{n_ph}_{n_rh}_{n_c}.json",
                    verbose=False)
        results.append((label, r))

    # Summary
    print()
    print("=" * 72, flush=True)
    print(f"Summary order={args.order}:", flush=True)
    print("=" * 72, flush=True)
    print(f"  {'case':>30s} | {'t_K [s]':>8s} | {'leading τ [μs]':>16s}",
          flush=True)
    for label, r in results:
        print(f"  {label:>30s} | {r['t_K_assembly_s']:>8.2f} | "
              f"{r['leading_tau_us']:>16.6f}",
              flush=True)
    if len(results) == 2:
        base = results[0][1]["leading_tau_us"]
        doub = results[1][1]["leading_tau_us"]
        delta = doub - base
        print(f"\n  Δ (doubled - baseline) = {delta:.6f} μs ({delta/base*100:.4f}%)",
              flush=True)


if __name__ == "__main__":
    main()
