"""bench_peec_mna_crossover.py -- Step 4 stage 5

Measures the crossover N at which the HACApK-MNA path
(PEECCircuitSolver(..., use_hacapk=True)) beats the dense C++ MNASolver
on the full nodal impedance solve.

FINDINGS (LAB, 2026-04-16)
--------------------------
Stage 5 showed that naive nested BiCGSTAB (scipy bicgstab outer +
HACApK BiCGSTAB inner) stalls above N ~ 200 on a strongly coupled
helical coil.  Stage 6 (this revision) replaces the outer BiCGSTAB
with scipy.sparse.linalg.lgmres, which tolerates the varying-operator
preconditioner introduced by the inner Krylov -- so convergence is
now robust across DC..10 MHz and N up to ~2000.

However, *robust* does not yet mean *fast*:

  * Dense LU C++ @ N=1800:      1.07 s    (reference)
  * HACApK-MNA lgmres @ N=1800:   see results_bench_peec_mna_crossover.json

Even with a flexible outer, each outer matvec triggers a full inner
HACApK BiCGSTAB (20-50 real matvecs at O(N log N)), so the overall
wall time per port impedance solve remains dominated by the inner
cost.  The stage-6 lgmres path still does not beat dense LU at
benchmarked sizes; the expected crossover is deferred to stage 7
(saddle-point direct formulation or H-matrix assembly of Y).

Invariants preserved vs stage 5:
  - Pure L matvec crossover at N ~ 500 (stage 3, unchanged)
  - Validate at N=72: rel err ~3.6e-08 across DC..10 MHz (still PASS)
  - `bicgstab` outer remains available via outer_method=... for
    comparison, even though it stalls.

This script is the stage-6 baseline.  Stage 7 changes will re-run
it and diff against results_bench_peec_mna_crossover.json.

Geometry mirrors bench_peec_hacapk.py exactly (20 turns x 25 seg/turn
helical solenoid, nwinc x nhinc sub-filaments) so rows match on N_fil.

Per case (nwinc x mode) the subprocess reports:
  t_topology    PEECBuilder + build_topology (shared, same for both modes)
  t_solver_init PEECCircuitSolver(...) constructor
  t_solve       compute_port_impedance(f) at 10 kHz
  Z11           sanity check (same value across modes within rtol)
  peak_memory_mb psutil peak_wset (captures C++ heap + H-matrix)
  mode-specific stats:
    dense:      n/a (single C++ LU solve, no iters)
    hacapk:     outer_iters, inner_iters_total, inner_matvec_real

One subprocess per (nwinc, mode) pair for clean memory measurement.

Usage:
    python bench_peec_mna_crossover.py                              # run all
    python bench_peec_mna_crossover.py --case 0 --mode dense        # subprocess
    python bench_peec_mna_crossover.py --case 0 --mode hacapk
    python bench_peec_mna_crossover.py --ncases

Output: validation_test/solver_benchmarks/results_bench_peec_mna_crossover.json.
"""

from __future__ import annotations

import argparse
import importlib.metadata as im
import json
import math
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np
import psutil


# ----------------------------------------------------------------------
# Benchmark configuration
# ----------------------------------------------------------------------

BASE_TURNS = 10
SEG_PER_TURN = 20                    # -> 200 connected segments
COIL_RADIUS_M = 0.050
COIL_PITCH_M = 0.006
WIRE_W_M = 0.003
WIRE_H_M = 0.003
SIGMA_CU = 5.8e7
FREQ_HZ = 10.0e3

# nwinc = nhinc each case; N_fil = 200 * nwinc^2 (with base above)
# Dense covers [1..5] (up to 5000 filaments).  Stage-6b saddle point
# is much faster than the stage-6a nested lgmres so extending HACApK
# coverage to [1..5] to match dense and expose the crossover.
SUBDIVISIONS_DENSE = [1, 2, 3, 4, 5]
SUBDIVISIONS_HACAPK = [1, 2, 3, 4, 5]

MODES = ("dense", "hacapk")

# Outer Krylov method for the HACApK MNA path.
# "saddle"  : stage-6b block saddle-point (default). Single lgmres on
#             the combined [I;V] system with block-diag preconditioner
#             (diag(Z)^-1 on upper, sparse LU of Schur on lower).
#             No nested Krylov.  ~50 matvecs per port solve at N=72.
# "lgmres"  : stage-6a nested (lgmres outer, BiCGSTAB inner).  Robust
#             but slow; retained for comparison.
# "gcrotmk" : 10x faster than lgmres when it works, but NaN-divides
#             at N=72 + 10 kHz on helix.
# "bicgstab": stalls above N ~ 200 (stage-5 finding).
OUTER_METHOD = "saddle"

# Which outer-method schedule to run.  Set to a list to sweep methods
# in a single bench pass (useful for paper plots).
OUTER_METHOD_SWEEP = None   # e.g. ["lgmres", "saddle"] to compare

# Tolerances:
#   Dense C++ MNA uses LU -- tolerances ignored.
#   HACApK Python MNA: loose tolerances because the inner BiCGSTAB
#   breakdown floor (residual around 1e-7 on a strongly coupled
#   helical coil) bounds the usable outer tolerance from below.
BICG_OUTER_RTOL = 1.0e-4
BICG_INNER_RTOL = 1.0e-8
BICG_MAXITER = 5000

BENCH_NAME = "bench_peec_mna_crossover"
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESULT_DIR = os.path.join(REPO_ROOT, "validation_test", "solver_benchmarks")
OUTFILE = os.path.join(
    RESULT_DIR,
    f"results_{BENCH_NAME}.json",
)


# ----------------------------------------------------------------------
# Geometry: helical solenoid via PEECBuilder (matches bench_peec_hacapk.py
# sub-filament layout because PEECMatrixBuilder performs the same
# nwinc x nhinc cross-section expansion internally).
# ----------------------------------------------------------------------

def _ensure_radia_on_path() -> None:
    """Insert repo's src/ and src/radia on sys.path and trigger MKL load."""
    here = os.path.dirname(os.path.abspath(__file__))
    root_src = os.path.join(here, "..", "..", "src")
    inner = os.path.join(root_src, "radia")
    for p in (root_src, inner):
        if p not in sys.path:
            sys.path.insert(0, p)
    import radia  # noqa: F401 -- MKL DLL registration


def _build_topology(nwinc: int):
    """Build the per-filament topology dict via PEECBuilder."""
    _ensure_radia_on_path()
    from peec_matrices import PEECBuilder

    builder = PEECBuilder()
    nodes = []
    for i in range(BASE_TURNS * SEG_PER_TURN + 1):
        t = 2.0 * math.pi * i / SEG_PER_TURN
        z = COIL_PITCH_M * (i / SEG_PER_TURN)
        nodes.append(builder.add_node_at(
            COIL_RADIUS_M * math.cos(t),
            COIL_RADIUS_M * math.sin(t),
            z,
        ))
    for i in range(BASE_TURNS * SEG_PER_TURN):
        builder.add_connected_segment(
            nodes[i], nodes[i + 1],
            WIRE_W_M, WIRE_H_M, sigma=SIGMA_CU,
            nwinc=nwinc, nhinc=nwinc,
        )
    builder.add_port(nodes[0], nodes[-1])
    return builder.build_topology()


# ----------------------------------------------------------------------
# Per-mode case runner
# ----------------------------------------------------------------------

def run_case(nwinc: int, mode: str) -> dict:
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")

    _ensure_radia_on_path()
    from peec_topology import PEECCircuitSolver

    proc = psutil.Process(os.getpid())

    # ---- Topology (shared cost) ----------------------------------
    t0 = time.perf_counter()
    topo = _build_topology(nwinc)
    t_topology = time.perf_counter() - t0

    n_fil = int(topo["n_loop"])
    n_nodes = int(topo["n_nodes"])

    # ---- Solver init ---------------------------------------------
    t0 = time.perf_counter()
    if mode == "dense":
        solver = PEECCircuitSolver(topo)
    else:
        solver = PEECCircuitSolver(
            topo,
            use_hacapk=True,
            hacapk_params={
                "aca_eps": 1.0e-8,
                "leaf_size": 128,
                "eta": 3.0,
                "max_rank": 400,
                "print_level": 0,
            },
            bicgstab_outer_tol=BICG_OUTER_RTOL,
            bicgstab_inner_tol=BICG_INNER_RTOL,
            bicgstab_maxiter=BICG_MAXITER,
            outer_method=OUTER_METHOD,
        )
    t_solver_init = time.perf_counter() - t0

    # ---- Port impedance solve ------------------------------------
    t0 = time.perf_counter()
    Z11 = solver.compute_port_impedance(FREQ_HZ)
    t_solve = time.perf_counter() - t0

    mem = proc.memory_info()
    peak_memory_mb = (
        mem.peak_wset / (1024 * 1024)
        if hasattr(mem, "peak_wset")
        else mem.rss / (1024 * 1024)
    )

    result = {
        "nwinc": nwinc,
        "nhinc": nwinc,
        "mode": mode,
        "n_fil": n_fil,
        "n_nodes": n_nodes,
        "freq_hz": FREQ_HZ,
        "t_topology": t_topology,
        "t_solver_init": t_solver_init,
        "t_solve": t_solve,
        "t_total_after_topology": t_solver_init + t_solve,
        "Z11_real": float(Z11.real),
        "Z11_imag": float(Z11.imag),
        "abs_Z11": float(abs(Z11)),
        "peak_memory_mb": peak_memory_mb,
        "converged": True,  # non-convergence raises
    }

    if mode == "hacapk":
        nodal_stats = getattr(solver, "_last_nodal_stats", {})
        hac_solver = solver._hacapk_solver
        hac_stats = hac_solver.hacapk_stats if hac_solver is not None else {}
        result.update({
            "outer_method": OUTER_METHOD,
            "hacapk_build_time": float(hac_stats.get("build_time", -1.0)),
            "hmatrix_mb": float(hac_stats.get("memory_mb", -1.0)),
            "compression": float(hac_stats.get("compression", -1.0)),
            "outer_iters": int(nodal_stats.get("outer", -1)),
            "inner_iters_total": int(nodal_stats.get("inner_iters", -1)),
            "inner_matvec_real": int(nodal_stats.get("inner_matvec_real", -1)),
        })

    return result


# ----------------------------------------------------------------------
# Driver: one subprocess per (nwinc, mode)
# ----------------------------------------------------------------------

def run_all() -> None:
    # Build the (nwinc, mode) schedule.  HACApK cases are restricted to
    # SUBDIVISIONS_HACAPK because nested Krylov stalls past N ~ 500.
    schedule = []
    for nwinc in SUBDIVISIONS_DENSE:
        schedule.append((nwinc, "dense"))
    for nwinc in SUBDIVISIONS_HACAPK:
        schedule.append((nwinc, "hacapk"))
    n_cases = len(schedule)

    print(f"[{BENCH_NAME}] {n_cases} subprocess runs")
    print(f"[{BENCH_NAME}] base: {BASE_TURNS} turns x {SEG_PER_TURN} seg/turn "
          f"= {BASE_TURNS * SEG_PER_TURN} connected segments")
    print(f"[{BENCH_NAME}] freq: {FREQ_HZ * 1e-3:.1f} kHz, "
          f"outer/inner rtol={BICG_OUTER_RTOL:.0e}/{BICG_INNER_RTOL:.0e}")
    print()

    # For --case indexing we union the two subdivision lists and pass
    # the union index; the subprocess picks up nwinc from that index.
    all_nwincs = sorted(set(SUBDIVISIONS_DENSE) | set(SUBDIVISIONS_HACAPK))
    script = os.path.abspath(__file__)
    results = []
    for done, (nwinc, mode) in enumerate(schedule, start=1):
        n_fil_expected = BASE_TURNS * SEG_PER_TURN * nwinc * nwinc
        dense_gb = n_fil_expected ** 2 * 8 / (1024 ** 3)
        print(f"  [{done}/{n_cases}] nwinc={nwinc} mode={mode}  "
              f"N_fil~{n_fil_expected}  (dense L ~ {dense_gb:.2f} GB)",
              flush=True)
        t0 = time.perf_counter()
        pr = subprocess.run(
            [sys.executable, script,
             "--case", str(all_nwincs.index(nwinc)),
             "--mode", mode,
             "--json-stdout"],
            capture_output=True, text=True,
        )
        wall = time.perf_counter() - t0
        if pr.returncode != 0:
            print(f"     FAILED (rc={pr.returncode}, wall={wall:.1f}s)")
            print("     stderr:", pr.stderr[-500:])
            results.append({
                "nwinc": nwinc, "mode": mode,
                "error": pr.stderr[-2000:],
                "wall_time": wall, "converged": False,
            })
            continue
        last_json = None
        for line in pr.stdout.splitlines():
            line = line.strip()
            if line.startswith("{"):
                last_json = line
        if last_json is None:
            print(f"     FAILED: no JSON in stdout")
            results.append({
                "nwinc": nwinc, "mode": mode,
                "error": "no JSON",
                "converged": False,
            })
            continue
        res = json.loads(last_json)
        res["wall_time"] = wall
        results.append(res)
        if mode == "hacapk":
            print(f"     N={res['n_fil']}  "
                  f"t_init={res['t_solver_init']:.2f}s  "
                  f"t_solve={res['t_solve']:.2f}s  "
                  f"outer={res.get('outer_iters', '?')}  "
                  f"inner_real={res.get('inner_matvec_real', '?')}  "
                  f"peak={res['peak_memory_mb']:.0f}MB  "
                  f"|Z11|={res['abs_Z11']:.3e}")
        else:
            print(f"     N={res['n_fil']}  "
                  f"t_init={res['t_solver_init']:.2f}s  "
                  f"t_solve={res['t_solve']:.2f}s  "
                  f"peak={res['peak_memory_mb']:.0f}MB  "
                  f"|Z11|={res['abs_Z11']:.3e}")

    generated_at_utc = datetime.now(timezone.utc).isoformat()
    data = {
        "schema": "radia.validation.peec-mna-crossover-benchmark.v1",
        "generated_at_utc": generated_at_utc,
        "timestamp": generated_at_utc,
        "hostname": platform.node(),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "radia": im.version("radia"),
        },
        "benchmark": BENCH_NAME,
        "problem": {
            "base_turns": BASE_TURNS,
            "seg_per_turn": SEG_PER_TURN,
            "n_connected_segments": BASE_TURNS * SEG_PER_TURN,
            "coil_radius_m": COIL_RADIUS_M,
            "coil_pitch_m": COIL_PITCH_M,
            "wire_w_m": WIRE_W_M,
            "wire_h_m": WIRE_H_M,
            "sigma_S_per_m": SIGMA_CU,
            "freq_hz": FREQ_HZ,
            "subdivisions_dense": SUBDIVISIONS_DENSE,
            "subdivisions_hacapk": SUBDIVISIONS_HACAPK,
            "modes": list(MODES),
            "outer_method": OUTER_METHOD,
            "bicg_outer_rtol": BICG_OUTER_RTOL,
            "bicg_inner_rtol": BICG_INNER_RTOL,
            "bicg_maxiter": BICG_MAXITER,
        },
        "results": results,
    }
    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
    with open(OUTFILE, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print()
    print(f"[{BENCH_NAME}] results saved to {OUTFILE}")

    _print_crossover_table(results)


def _print_crossover_table(results: list) -> None:
    """Plain-text crossover table (also persisted inside the JSON)."""
    by_n = {}
    for r in results:
        if r.get("converged"):
            by_n.setdefault(int(r["n_fil"]), {})[r["mode"]] = r

    if not by_n:
        return
    print()
    print("  N_fil    dense t_solve    hacapk t_solve    speedup    "
          "dense peak MB    hacapk peak MB    |dZ11/Z11|")
    print("  -----    -------------    --------------    -------    "
          "-------------    --------------    ----------")
    for n in sorted(by_n):
        dense = by_n[n].get("dense")
        hac = by_n[n].get("hacapk")
        if dense is None or hac is None:
            continue
        sp = dense["t_solve"] / hac["t_solve"] if hac["t_solve"] > 0 else 0.0
        dZ = abs(complex(hac["Z11_real"], hac["Z11_imag"])
                 - complex(dense["Z11_real"], dense["Z11_imag"]))
        rel = dZ / abs(complex(dense["Z11_real"], dense["Z11_imag"]))
        print(f"  {n:5d}    {dense['t_solve']:9.2f} s      "
              f"{hac['t_solve']:9.2f} s     {sp:6.1f}x    "
              f"{dense['peak_memory_mb']:9.0f}        "
              f"{hac['peak_memory_mb']:9.0f}        "
              f"{rel:.2e}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, default=None,
                        help="Run one case index (subprocess mode)")
    parser.add_argument("--mode", choices=MODES, default=None,
                        help="dense or hacapk (subprocess mode)")
    parser.add_argument("--json-stdout", action="store_true")
    parser.add_argument("--ncases", action="store_true")
    args = parser.parse_args()

    all_nwincs = sorted(set(SUBDIVISIONS_DENSE) | set(SUBDIVISIONS_HACAPK))
    if args.ncases:
        print(len(SUBDIVISIONS_DENSE) + len(SUBDIVISIONS_HACAPK))
        return

    if args.case is not None:
        if args.mode is None:
            parser.error("--mode is required with --case")
        nwinc = all_nwincs[args.case]
        res = run_case(nwinc, args.mode)
        if args.json_stdout:
            print(json.dumps(res))
        else:
            print(json.dumps(res, indent=2))
        return

    run_all()


if __name__ == "__main__":
    main()
