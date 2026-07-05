"""Benchmark: (ACA+)+TSVD vs naive dense TSVD on the SAME matrix kernel.

Both methods produce the same truncated SVD  A ~= U diag(S) V^T, but:

  naive TSVD : build the full dense A -- M*N kernel evaluations -- then take a
               full numpy.linalg.svd (O(N M^2)).
  (ACA+)+TSVD: ACA+ evaluates only ~k_aca*(M+N) kernel entries, then TSVDs the
               small factors (negligible).

So there are two wins, both shown here:
  1. kernel-evaluation count:   M*N        ->  ~k_aca*(M+N)
  2. dense SVD cost:            O(N M^2)   ->  small-factor SVDs

The matrix entry is computed per call by the SAME smooth (numerically low-rank)
kernel for both methods, so the comparison is apples-to-apples and the
eval-count reduction shows up in wall-clock time when the kernel is non-trivial
(as a real Biot-Savart / magnetic-material field kernel is).

Output: results_aca_vs_dense.json (per the repository benchmark policy).
Run:    python bench_aca_vs_dense.py [--max-n N]
"""
import argparse
import json
import os
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import psutil

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from radia.stream_function import aca_tsvd

_N_EVAL = 0  # global kernel-evaluation counter


def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return (mem.peak_wset / (1024 * 1024) if hasattr(mem, "peak_wset")
            else mem.rss / (1024 * 1024))


def make_kernel(obs, src, alpha):
    """A(i,j) = 1/(1 + alpha |obs_i - src_j|^2): smooth, numerically low rank.
    Computed per call (counts as one kernel evaluation) so the eval-count
    difference between the two methods is reflected in the wall-clock time."""
    def kernel(i, j):
        global _N_EVAL
        _N_EVAL += 1
        d = obs[i] - src[j]
        return 1.0 / (1.0 + alpha * float(d @ d))
    return kernel


def run_case(N, modes, alpha=400.0, aca_eps=1.0e-8, seed=0):
    global _N_EVAL
    M = max(4, N // 4)
    rng = np.random.default_rng(seed)
    obs = np.column_stack([rng.uniform(-0.1, 0.1, M), rng.uniform(-0.1, 0.1, M),
                           np.full(M, 0.30)])
    src = np.column_stack([rng.uniform(-0.1, 0.1, N), rng.uniform(-0.1, 0.1, N),
                           np.zeros(N)])
    kernel = make_kernel(obs, src, alpha)

    # --- naive dense TSVD: build full A via the kernel, then full SVD ---
    _N_EVAL = 0
    t0 = time.perf_counter()
    A = np.empty((M, N))
    for i in range(M):
        for j in range(N):
            A[i, j] = kernel(i, j)
    t_build = time.perf_counter() - t0
    n_eval_naive = _N_EVAL                      # == M * N
    t0 = time.perf_counter()
    Ud, Sd, Vtd = np.linalg.svd(A, full_matrices=False)
    t_svd = time.perf_counter() - t0
    t_naive = t_build + t_svd
    mem_naive = get_peak_memory_mb()

    # --- (ACA+)+TSVD: ACA evaluates only k_aca*(M+N)-ish entries ---
    _N_EVAL = 0
    t0 = time.perf_counter()
    res = aca_tsvd(M, N, kernel, modes=modes, kmax=min(M, N),
                   aca_eps=aca_eps, method=3)
    t_aca = time.perf_counter() - t0
    n_eval_aca = _N_EVAL
    mem_aca = get_peak_memory_mb()

    ns = min(res.k_aca, modes, len(Sd))
    s_rel = float(np.linalg.norm(Sd[:ns] - res.S[:ns]) / np.linalg.norm(Sd[:ns]))
    eval_ratio = n_eval_naive / max(1, n_eval_aca)
    speedup = t_naive / t_aca if t_aca > 0 else float("inf")

    print(f"N={N:5d} M={M:4d} k_aca={res.k_aca:3d} | "
          f"evals naive={n_eval_naive:>8d} aca={n_eval_aca:>7d} "
          f"({eval_ratio:5.1f}x) | "
          f"t naive={t_naive*1e3:8.1f}ms aca={t_aca*1e3:7.1f}ms "
          f"({speedup:5.1f}x) | dS/S={s_rel:.1e}")

    return [
        {"method": "naive_dense_tsvd", "N": N, "M": M, "modes": modes,
         "t_setup": t_build, "t_solve": t_svd, "t_total": t_naive,
         "iterations": int(min(M, N)), "converged": True,
         "peak_memory_mb": mem_naive, "n_kernel_evals": n_eval_naive,
         "k_aca": int(min(M, N))},
        {"method": "aca_tsvd", "N": N, "M": M, "modes": modes,
         "t_setup": 0.0, "t_solve": t_aca, "t_total": t_aca,
         "iterations": int(res.k_aca), "converged": True,
         "peak_memory_mb": mem_aca, "n_kernel_evals": n_eval_aca,
         "k_aca": int(res.k_aca), "eval_reduction_vs_naive": eval_ratio,
         "speedup_vs_naive": speedup,
         "singular_value_rel_err_vs_naive": s_rel}]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=2048)
    args = ap.parse_args()

    print("=== (ACA+)+TSVD vs naive dense TSVD (same per-call kernel) ===")
    # Warm up (first aca_tsvd call pays one-off pybind/HACApK setup).
    _ow = np.zeros((4, 3))
    aca_tsvd(4, 8, lambda i, j: 1.0 / (1.0 + i + j), modes=2, kmax=4, aca_eps=1e-8)

    results = []
    N = 256
    while N <= args.max_n:
        results.extend(run_case(N, modes=20))
        N *= 2

    data = {
        "timestamp": datetime.now().isoformat(),
        "hostname": platform.node(),
        "benchmark": "stream_function_aca_vs_naive_tsvd",
        "problem": {"kernel": "1/(1+alpha r^2)", "alpha": 400.0,
                    "M_over_N": 0.25, "modes": 20, "aca_eps": 1.0e-8,
                    "method": 3},
        "results": results,
    }
    out = HERE / "results_aca_vs_dense.json"
    with open(out, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\nResults saved to {out.name}")


if __name__ == "__main__":
    main()
