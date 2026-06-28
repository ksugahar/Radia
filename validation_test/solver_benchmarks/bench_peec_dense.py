"""
bench_peec_dense.py -- Baseline for HACApK-PEEC (Phase 1)

Dense PEEC L assembly + LAPACK LU port solve on an induction-heating-style
helical solenoid. Scales filament count N via nwinc x nhinc subdivision of a
fixed base geometry (500 connected segments, 20 turns x 25 seg/turn).

Measures, per case:
  - t_setup   : geometry construction (add_node_at + add_connected_segment)
  - t_build   : dense L assembly via PEECBuilder.build() [core HACApK target]
  - t_solve   : single-frequency LAPACK LU port impedance solve
  - iterations: always 1 (direct LU, no Krylov)
  - converged : True if Z is finite
  - peak_memory_mb : psutil peak_wset (captures C++ heap)

Per repository benchmark policy:
  - One process per case (fresh Python interpreter)
  - JSON output: docs/solver_benchmarks/results_bench_peec_dense.json

Baseline for Phase 1 HACApK-PEEC paper (see docs/research/HACAPK_PEEC_PRIMA_PAPER.md).
Run this once to establish the dense-PEEC operating point, then compare
after the HACApK-PEEC adapter is wired in.

Usage:
    python bench_peec_dense.py                 # run all cases, save JSON
    python bench_peec_dense.py --case 1        # run case index 1 only (subprocess mode)
    python bench_peec_dense.py --ncases        # print number of cases and exit
"""

import argparse
import json
import math
import os
import platform
import subprocess
import sys
import time
from datetime import datetime

import numpy as np
import psutil


# ----------------------------------------------------------------------
# Benchmark configuration
# ----------------------------------------------------------------------

BASE_TURNS = 20
SEG_PER_TURN = 25                      # -> 500 connected segments
COIL_RADIUS_M = 0.050                  # 100 mm diameter
COIL_PITCH_M = 0.006                   # 6 mm turn pitch
WIRE_W_M = 0.003                       # 3 mm square cross-section
WIRE_H_M = 0.003
SIGMA_CU = 5.8e7
FREQ_HZ = 10.0e3                       # typical IH frequency

# nwinc = nhinc for each case; N_fil = 500 * nwinc^2
SUBDIVISIONS = [1, 2, 3, 4, 5]

BENCH_NAME = "bench_peec_dense"
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESULT_DIR = os.path.join(REPO_ROOT, "docs", "solver_benchmarks")
OUTFILE = os.path.join(
    RESULT_DIR,
    f"results_{BENCH_NAME}.json",
)


# ----------------------------------------------------------------------
# Single-case runner
# ----------------------------------------------------------------------

def run_case(nwinc: int) -> dict:
    """Run one benchmark case in the current (fresh) process."""
    from radia.peec_matrices import PEECBuilder

    proc = psutil.Process(os.getpid())

    # ---- Setup: helical geometry ----------------------------------
    t0 = time.perf_counter()

    builder = PEECBuilder()
    nsteps = BASE_TURNS * SEG_PER_TURN
    prev = builder.add_node_at(COIL_RADIUS_M, 0.0, 0.0)
    first_node = prev
    for i in range(1, nsteps + 1):
        theta = 2 * math.pi * i / SEG_PER_TURN
        z = COIL_PITCH_M * (i / SEG_PER_TURN)
        x = COIL_RADIUS_M * math.cos(theta)
        y = COIL_RADIUS_M * math.sin(theta)
        curr = builder.add_node_at(x, y, z)
        builder.add_connected_segment(
            prev, curr, WIRE_W_M, WIRE_H_M,
            sigma=SIGMA_CU, nwinc=nwinc, nhinc=nwinc,
        )
        prev = curr
    last_node = prev
    builder.add_port(first_node, last_node)

    t_setup = time.perf_counter() - t0

    # ---- Build: dense L, R_dc assembly ----------------------------
    t0 = time.perf_counter()
    L, R_dc, P, M_LS = builder.build(include_star=False)
    t_build = time.perf_counter() - t0

    n_fil = L.shape[0]
    L_bytes = L.size * L.itemsize

    # ---- Solve: LAPACK LU port impedance --------------------------
    t0 = time.perf_counter()
    Z_matrix = builder.solve_z_matrix(FREQ_HZ, None)
    t_solve = time.perf_counter() - t0
    Z11 = complex(Z_matrix[0, 0])

    # ---- Peak memory -----------------------------------------------
    mem = proc.memory_info()
    peak_memory_mb = (
        mem.peak_wset / (1024 * 1024)
        if hasattr(mem, "peak_wset")
        else mem.rss / (1024 * 1024)
    )

    converged = np.isfinite(Z11.real) and np.isfinite(Z11.imag)

    return {
        "nwinc": nwinc,
        "nhinc": nwinc,
        "n_connected_segments": nsteps,
        "n_fil": int(n_fil),
        "L_matrix_bytes": int(L_bytes),
        "L_matrix_mb": L_bytes / (1024 * 1024),
        "freq_hz": FREQ_HZ,
        "t_setup": t_setup,
        "t_build": t_build,
        "t_solve": t_solve,
        "iterations": 1,
        "converged": bool(converged),
        "peak_memory_mb": peak_memory_mb,
        "Z11_real": Z11.real,
        "Z11_imag": Z11.imag,
        "abs_Z11": abs(Z11),
    }


# ----------------------------------------------------------------------
# Driver: spawn one subprocess per case for clean memory measurement
# ----------------------------------------------------------------------

def run_all() -> None:
    print(f"[{BENCH_NAME}] running {len(SUBDIVISIONS)} cases (one subprocess each)")
    print(f"[{BENCH_NAME}] base: {BASE_TURNS} turns x {SEG_PER_TURN} seg/turn "
          f"= {BASE_TURNS * SEG_PER_TURN} connected segments")
    print(f"[{BENCH_NAME}] freq: {FREQ_HZ * 1e-3:.1f} kHz")
    print()

    script = os.path.abspath(__file__)
    results = []
    for idx, nwinc in enumerate(SUBDIVISIONS):
        n_fil_expected = BASE_TURNS * SEG_PER_TURN * nwinc * nwinc
        print(f"  [{idx + 1}/{len(SUBDIVISIONS)}] nwinc={nwinc}  "
              f"N_fil~{n_fil_expected}  (L ~ {n_fil_expected**2 * 8 / (1024**3):.2f} GB real64)",
              flush=True)
        t0 = time.perf_counter()
        proc = subprocess.run(
            [sys.executable, script, "--case", str(idx), "--json-stdout"],
            capture_output=True, text=True,
        )
        wall = time.perf_counter() - t0
        if proc.returncode != 0:
            print(f"     FAILED (rc={proc.returncode}, wall={wall:.1f}s)")
            print("     stderr:", proc.stderr[-500:])
            results.append({
                "nwinc": nwinc,
                "error": proc.stderr[-2000:],
                "wall_time": wall,
                "converged": False,
            })
            continue
        # Last JSON line is the result
        last_json = None
        for line in proc.stdout.splitlines():
            line = line.strip()
            if line.startswith("{"):
                last_json = line
        if last_json is None:
            print(f"     FAILED: no JSON in stdout")
            results.append({"nwinc": nwinc, "error": "no JSON", "converged": False})
            continue
        res = json.loads(last_json)
        res["wall_time"] = wall
        results.append(res)
        print(f"     N={res['n_fil']}  t_build={res['t_build']:.2f}s  "
              f"t_solve={res['t_solve']:.2f}s  peak_mem={res['peak_memory_mb']:.0f} MB  "
              f"|Z|={res['abs_Z11']*1e3:.3f} mOhm")

    data = {
        "timestamp": datetime.now().isoformat(),
        "hostname": platform.node(),
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
            "subdivisions": SUBDIVISIONS,
        },
        "results": results,
    }
    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
    with open(OUTFILE, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print()
    print(f"[{BENCH_NAME}] results saved to {OUTFILE}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, default=None,
                        help="Run a single case by index (subprocess mode)")
    parser.add_argument("--json-stdout", action="store_true",
                        help="Emit the result as a single JSON line on stdout")
    parser.add_argument("--ncases", action="store_true",
                        help="Print number of cases and exit")
    args = parser.parse_args()

    if args.ncases:
        print(len(SUBDIVISIONS))
        return

    if args.case is not None:
        nwinc = SUBDIVISIONS[args.case]
        res = run_case(nwinc)
        if args.json_stdout:
            print(json.dumps(res))
        else:
            print(json.dumps(res, indent=2))
        return

    run_all()


if __name__ == "__main__":
    main()
