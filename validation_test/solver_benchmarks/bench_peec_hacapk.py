"""
bench_peec_hacapk.py -- HACApK-PEEC scaling benchmark (Phase 1)

Mirrors bench_peec_dense.py on the same helical solenoid geometry (20 turns
x 25 seg/turn, nwinc = nhinc subdivision) but builds the filament mutual
inductance L as an H-matrix via HACApKPEECManager and solves the frequency-
domain PEEC system

    (diag(R_dc) + j*omega*L) @ I = V

with complex BiCGSTAB (Option A: split real/imag, two real HACApK matvecs
per complex matvec).  No Zs term here -- DC resistance only, so that the
comparison vs bench_peec_dense.py isolates the linear-algebra cost.

Per repository benchmark policy:
  - one subprocess per case (fresh python, clean peak-memory measurement)
  - JSON output: docs/solver_benchmarks/results_bench_peec_hacapk.json

Subdivisions include values beyond the dense-LU feasibility window
(nwinc in [1..7] -> N_fil = 500, 2000, 4500, 8000, 12500, 18000, 24500)
to show the HACApK-only regime.

Baseline for HACApK-PEEC paper (see docs/research/HACAPK_PEEC_PRIMA_PAPER.md).

Usage:
    python bench_peec_hacapk.py                 # run all cases, save JSON
    python bench_peec_hacapk.py --case 0        # run one case (subprocess mode)
    python bench_peec_hacapk.py --ncases        # print number of cases
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
# Benchmark configuration (matches bench_peec_dense.py geometry)
# ----------------------------------------------------------------------

BASE_TURNS = 20
SEG_PER_TURN = 25                      # -> 500 connected segments
COIL_RADIUS_M = 0.050
COIL_PITCH_M = 0.006
WIRE_W_M = 0.003
WIRE_H_M = 0.003
SIGMA_CU = 5.8e7
FREQ_HZ = 10.0e3

# nwinc = nhinc for each case; N_fil = 500 * nwinc^2
# Extended beyond dense (which stops at 5) to exercise the HACApK-only regime.
SUBDIVISIONS = [1, 2, 3, 4, 5, 6, 7]

# BiCGSTAB convergence target
BICG_RTOL = 1.0e-6
BICG_MAXITER = 2000

BENCH_NAME = "bench_peec_hacapk"
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESULT_DIR = os.path.join(REPO_ROOT, "docs", "solver_benchmarks")
OUTFILE = os.path.join(
    RESULT_DIR,
    f"results_{BENCH_NAME}.json",
)


# ----------------------------------------------------------------------
# Geometry: helical solenoid with nwinc x nhinc cross-section subdivision
# ----------------------------------------------------------------------

def _build_helical_filaments(nwinc: int) -> dict:
    """
    Build the raw filament arrays for a helical solenoid with nwinc x nhinc
    cross-section subdivision.  Subdivision logic mirrors PEECMatrixBuilder:

      - parent segment centers lie on a helix
      - each parent segment is replaced by nwinc*nhinc parallel filaments
        spaced by (width/nwinc, height/nhinc) in the cross-section frame
      - local frame: e_w = (ref x dir)/|...|, e_h = dir x e_w
                     ref = (1,0,0) unless dir.x close to 1
      - offset_w = (width/nwinc) * (iw - (nwinc-1)/2)
      - offset_h = (height/nhinc) * (ih - (nhinc-1)/2)
    """
    nw = nh = int(nwinc)
    nsteps = BASE_TURNS * SEG_PER_TURN

    # Parent segment centers / directions
    parents = []
    for i in range(nsteps):
        t0 = 2.0 * math.pi * i / SEG_PER_TURN
        t1 = 2.0 * math.pi * (i + 1) / SEG_PER_TURN
        z0 = COIL_PITCH_M * (i / SEG_PER_TURN)
        z1 = COIL_PITCH_M * ((i + 1) / SEG_PER_TURN)
        p0 = np.array([COIL_RADIUS_M * math.cos(t0),
                       COIL_RADIUS_M * math.sin(t0), z0])
        p1 = np.array([COIL_RADIUS_M * math.cos(t1),
                       COIL_RADIUS_M * math.sin(t1), z1])
        d = p1 - p0
        L = float(np.linalg.norm(d))
        center = 0.5 * (p0 + p1)
        direction = d / L
        parents.append((center, direction, L))

    # Expand to sub-filaments
    N = nsteps * nw * nh
    centers = np.zeros((N, 3), dtype=np.float64)
    directions = np.zeros((N, 3), dtype=np.float64)
    lengths = np.zeros(N, dtype=np.float64)
    widths = np.full(N, WIRE_W_M / nw, dtype=np.float64)
    heights = np.full(N, WIRE_H_M / nh, dtype=np.float64)
    sigmas = np.full(N, SIGMA_CU, dtype=np.float64)

    sub_w = WIRE_W_M / nw
    sub_h = WIRE_H_M / nh

    k = 0
    for (c, d, L) in parents:
        if abs(d[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])
        e_w = np.cross(ref, d)
        e_w /= np.linalg.norm(e_w)
        e_h = np.cross(d, e_w)

        if nw <= 1 and nh <= 1:
            centers[k] = c
            directions[k] = d
            lengths[k] = L
            widths[k] = WIRE_W_M
            heights[k] = WIRE_H_M
            k += 1
            continue

        for iw in range(nw):
            for ih in range(nh):
                off_w = sub_w * (iw - (nw - 1) * 0.5)
                off_h = sub_h * (ih - (nh - 1) * 0.5)
                centers[k] = c + off_w * e_w + off_h * e_h
                directions[k] = d
                lengths[k] = L
                k += 1

    assert k == N
    return dict(centers=centers, directions=directions,
                lengths=lengths, widths=widths,
                heights=heights, sigmas=sigmas)


# ----------------------------------------------------------------------
# Single-case runner
# ----------------------------------------------------------------------

def run_case(nwinc: int) -> dict:
    """Run one HACApK benchmark case in the current (fresh) process."""
    from radia.peec_hacapk_solver import PEECHACApKSolver, _default_R_dc
    from scipy.sparse.linalg import LinearOperator

    proc = psutil.Process(os.getpid())

    # ---- Setup: build filament geometry ----------------------------
    t0 = time.perf_counter()
    geom = _build_helical_filaments(nwinc)
    t_setup = time.perf_counter() - t0

    centers = geom["centers"]
    directions = geom["directions"]
    lengths = geom["lengths"]
    widths = geom["widths"]
    heights = geom["heights"]
    sigmas = geom["sigmas"]
    N = centers.shape[0]

    # ---- Build: HACApK H-matrix for L ------------------------------
    t0 = time.perf_counter()
    solver = PEECHACApKSolver(
        centers, directions, lengths, widths, heights, sigmas,
        print_level=0,  # use PEEC-tuned defaults (aca_eps=1e-4, leaf=128, eta=3, max_rank=400)
    )
    t_build = time.perf_counter() - t0

    stats = solver.hacapk_stats

    # ---- RHS: uniform voltage drop per filament --------------------
    #   V_k = length_k  (unit electric field along each filament)
    # Complex RHS makes BiCGSTAB exercise the full complex matvec path.
    rhs = (lengths.astype(np.complex128)
           + 1j * np.zeros(N, dtype=np.complex128))

    # Jacobi preconditioner M = 1 / diag(A)
    # diag(A) = R_dc + j*omega * L_self_dc_diag
    # We don't have the L diagonal without an extra pass; use the cheap
    # approximation: L_ii ~ (mu0/2*pi) * length * (log(2*length / r_eq) - 1)
    # but accurately extracting diag(L) requires N HACApK matvecs on e_j.
    # Skip for small N (cheap); use only R_dc for large N.
    omega = 2.0 * math.pi * FREQ_HZ
    R_dc = _default_R_dc(lengths, widths, heights, sigmas)

    if N <= 2000:
        # Extract diag(L) by N matvecs (expensive but cheap at small N)
        diag_L = np.zeros(N, dtype=np.float64)
        e = np.zeros(N, dtype=np.float64)
        for j in range(N):
            e[:] = 0.0
            e[j] = 1.0
            diag_L[j] = solver.apply_L(e)[j]
        diag_A = R_dc.astype(np.complex128) + 1j * omega * diag_L
    else:
        # Approximate diag(L) via Grover self-inductance (rectangular cross-
        # section, closed form) so we still get a sensible preconditioner.
        #   L_self = (mu0/(2*pi)) * length * (log(2*length/(w+h)) + 0.5 + 0.2235*(w+h)/length)
        MU0 = 4.0e-7 * math.pi
        wh = widths + heights
        diag_L = (MU0 / (2.0 * math.pi)) * lengths * (
            np.log(2.0 * lengths / wh) + 0.5 + 0.2235 * wh / lengths
        )
        diag_A = R_dc.astype(np.complex128) + 1j * omega * diag_L

    M = LinearOperator(
        shape=(N, N),
        matvec=lambda r: np.asarray(r) / diag_A,
        dtype=np.complex128,
    )

    # ---- Solve: complex BiCGSTAB via HACApK matvec ----------------
    t0 = time.perf_counter()
    res = solver.solve(
        FREQ_HZ, rhs,
        rtol=BICG_RTOL, maxiter=BICG_MAXITER, M=M,
        track_residual=False,
    )
    t_solve = time.perf_counter() - t0

    # ---- Port impedance sanity check ------------------------------
    # For uniform voltage drop V_k = length_k * E, total voltage across
    # the series coil = sum(V_k) and total current flowing into port =
    # mean filament current * (nwinc*nhinc) (parallel filaments carry
    # equal share in series). So an approximate port impedance estimate:
    total_V = float(np.sum(lengths))           # V = E * total length, E = 1 V/m
    mean_I = np.sum(res.x) / N                 # average filament current
    # Z_port ~ V_total / I_series
    # I_series = mean_I * N (sum of parallel currents through any cross-section)
    # ... actually each cross-section has nwinc*nhinc parallel filaments
    # carrying (total I_port / (nw*nh)), and they are summed once per
    # connected segment.  Total sum over all filaments = I_port * N_base.
    # So I_port = sum(res.x) / N_base.
    N_base = BASE_TURNS * SEG_PER_TURN
    I_port = complex(np.sum(res.x) / N_base)
    Z11_approx = complex(total_V / I_port) if abs(I_port) > 0 else complex("nan")

    # ---- Peak memory ----------------------------------------------
    mem = proc.memory_info()
    peak_memory_mb = (
        mem.peak_wset / (1024 * 1024)
        if hasattr(mem, "peak_wset")
        else mem.rss / (1024 * 1024)
    )

    # ---- Dense L memory (reference, not allocated) -----------------
    dense_L_mb = (N * N * 8) / (1024 * 1024)

    return {
        "nwinc": nwinc,
        "nhinc": nwinc,
        "n_connected_segments": BASE_TURNS * SEG_PER_TURN,
        "n_fil": int(N),
        "dense_L_mb": dense_L_mb,
        "hmatrix_mb": float(stats["memory_mb"]),
        "compression": float(stats["compression"]),
        "n_leaves": int(stats["n_leaves"]),
        "n_lowrank": int(stats["n_lowrank"]),
        "n_dense_blocks": int(stats["n_dense"]),
        "max_rank": int(stats["max_rank"]),
        "freq_hz": FREQ_HZ,
        "t_setup": t_setup,
        "t_build": t_build,
        "hacapk_build_time": float(stats["build_time"]),
        "t_solve": t_solve,
        "iterations": int(res.iterations),
        "n_matvec_real": int(res.n_matvec_real),
        "n_matvec_complex": int(res.n_matvec_complex),
        "final_residual": float(res.final_residual),
        "converged": bool(res.converged),
        "info": int(res.info),
        "peak_memory_mb": peak_memory_mb,
        "Z11_real": float(Z11_approx.real),
        "Z11_imag": float(Z11_approx.imag),
        "abs_Z11": float(abs(Z11_approx)),
    }


# ----------------------------------------------------------------------
# Driver: spawn one subprocess per case for clean memory measurement
# ----------------------------------------------------------------------

def run_all() -> None:
    print(f"[{BENCH_NAME}] running {len(SUBDIVISIONS)} cases (one subprocess each)")
    print(f"[{BENCH_NAME}] base: {BASE_TURNS} turns x {SEG_PER_TURN} seg/turn "
          f"= {BASE_TURNS * SEG_PER_TURN} connected segments")
    print(f"[{BENCH_NAME}] freq: {FREQ_HZ * 1e-3:.1f} kHz, BiCGSTAB rtol={BICG_RTOL}")
    print()

    script = os.path.abspath(__file__)
    results = []
    for idx, nwinc in enumerate(SUBDIVISIONS):
        n_fil_expected = BASE_TURNS * SEG_PER_TURN * nwinc * nwinc
        dense_gb = n_fil_expected ** 2 * 8 / (1024 ** 3)
        print(f"  [{idx + 1}/{len(SUBDIVISIONS)}] nwinc={nwinc}  "
              f"N_fil~{n_fil_expected}  (dense L ~ {dense_gb:.2f} GB)",
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
        print(f"     N={res['n_fil']}  "
              f"t_build={res['t_build']:.2f}s  "
              f"H-mat={res['hmatrix_mb']:.0f}MB ({res['compression']*100:.1f}% of dense)  "
              f"t_solve={res['t_solve']:.2f}s  "
              f"iters={res['iterations']}  "
              f"res={res['final_residual']:.1e}  "
              f"peak={res['peak_memory_mb']:.0f}MB")

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
            "bicg_rtol": BICG_RTOL,
            "bicg_maxiter": BICG_MAXITER,
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
