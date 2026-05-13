"""extract_tau_A1_square_gpu.py — GPU port of extract_tau_A1_square.py.

Same structure as extract_tau_cuboid_521_gpu.py but A1 cuboid
dimensions (17.72 x 17.72 x 2 mm).
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

script_dir = Path(__file__).parent
legacy_dir = script_dir.parent.parent.parent / "examples" / "CLN" / "scripts" / "ngsolve_validation" / "legacy_fp64"
sys.path.insert(0, str(legacy_dir))
from hex_vim_cupy_kassembly import assemble_K_cupy

from extract_tau_cuboid_521 import (
    assemble_mass_matrix_python,
    assemble_b_vector_python,
    extract_foster,
)
import radia_vim


A_M = 17.72e-3
B_M = 17.72e-3
C_M = 2e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7


def run_one(order, n_th, n_ph, n_rh, n_c, save_path=None):
    print("=" * 72, flush=True)
    print(f" A1 cuboid 17.72x17.72x2 mm Cu order={order} "
          f"quad=({n_th},{n_ph},{n_rh},{n_c})", flush=True)
    print("=" * 72, flush=True)

    basis = radia_vim.HDivDivFreeHexBasis(order=order)
    n = basis.n_dofs

    print(f"Assembling K matrix (GPU Spherical Duffy)...", flush=True)
    t0 = time.time()
    K_bare = assemble_K_cupy(order, A_M, B_M, C_M,
                              n_th=n_th, n_ph=n_ph, n_rh=n_rh, n_c=n_c,
                              use_gpu=True, dtype=np.float64,
                              verbose=False)
    t_K = time.time() - t0
    print(f"K assembly (GPU): {t_K:.2f}s", flush=True)
    K_phys = (MU0 / (4 * math.pi)) * K_bare

    M_phys = assemble_mass_matrix_python(basis, A_M, B_M, C_M, n_gauss=8)
    b_vec = assemble_b_vector_python(basis, A_M, B_M, C_M, n_gauss=8)
    result = extract_foster(K_phys, M_phys, b_vec, SIGMA, n_show=12)
    print(f" Result: leading τ = {result['leading_tau_us']:.6f} μs",
          flush=True)

    out = {
        "shape": "A1",
        "order": order, "n_dofs": int(n),
        "quadrature": {"n_omega_theta": n_th, "n_omega_phi": n_ph,
                       "n_rho": n_rh, "n_c": n_c},
        "t_K_assembly_s": float(t_K),
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
    args = parser.parse_args()

    cases = [
        ("baseline (12,24,12,6)", 12, 24, 12, 6),
        ("intermed (16,32,16,8)", 16, 32, 16, 8),
    ]

    results = []
    for label, n_th, n_ph, n_rh, n_c in cases:
        print(f"\n### {label} ###", flush=True)
        r = run_one(args.order, n_th, n_ph, n_rh, n_c,
                    save_path=f"C:/Users/ADMINI~1/AppData/Local/Temp/1/A1_p{args.order}_"
                              f"q{n_th}_{n_ph}_{n_rh}_{n_c}.json")
        results.append((label, r))

    print()
    print("=" * 72, flush=True)
    print(f"A1 Summary order={args.order}:", flush=True)
    print("=" * 72, flush=True)
    print(f"  {'case':>30s} | {'t_K [s]':>8s} | {'leading τ [μs]':>16s}",
          flush=True)
    for label, r in results:
        print(f"  {label:>30s} | {r['t_K_assembly_s']:>8.2f} | "
              f"{r['leading_tau_us']:>16.6f}",
              flush=True)
    if len(results) == 2:
        base = results[0][1]["leading_tau_us"]
        inter = results[1][1]["leading_tau_us"]
        delta = inter - base
        print(f"\n  Δ (intermed - baseline) = {delta:+.6f} μs ({delta/base*100:+.4f}%)",
              flush=True)


if __name__ == "__main__":
    main()
