#!/usr/bin/env python
"""Per-step timing benchmark for vim.SolveHysteresis (B-input play, hex cube).

The measurement behind the SA-26-0xx section-7 %TODO(hysteresis-timing): the
same 1 m cube, +-200 kA/m, 20-step full loop as the committed goldens
(binput_hdiv_loop*.json), swept over hex mesh sizes.  What the numbers must
show is the CONSTANT-LHS reuse regime: the charge-Gram H-matrix build and the
mass-Riesz PARDISO factor are paid once (t_setup_s / charge_gram_wall_s),
after which every quasi-static step costs only CG applies + the batched
material update (per-step t_step_s ~ flat across the loop).

Benchmark Policy: publication timings run on hibino first, or on mdx only when
hibino is unavailable and the mdx CI queue is idle; LAB runs are correctness
smoke only. Each case runs in its OWN subprocess so peak memory
(psutil peak working set) is per-case accurate.

Usage:
  python bench_hysteresis_step.py --sizes 8,12,16,20        # compute-host driver
  python bench_hysteresis_step.py --sizes 4                 # LAB smoke
  python bench_hysteresis_step.py --case-n 12               # one case, JSON on stdout
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
# In-repo runs resolve radia from the source tree; a copied script (e.g. the
# C:\temp\ mdx deployment) falls back to the installed package.
if len(_HERE.parents) > 2:
    _SRC = _HERE.parents[2] / "src"
    if _SRC.is_dir() and str(_SRC) not in sys.path:
        sys.path.insert(0, str(_SRC))

MU0 = 4.0e-7 * math.pi
H0 = 200.0e3
N_RAMP, N_BRANCH = 4, 8          # the goldens' drive: virgin 4 + descending 8 + ascending 8


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def _synthetic_play(mu_r=200.0, irr1=0.30, irr2=0.15, eta=(0.0, 0.3, 0.7), rmax=4.0, n=81):
    # rmax must exceed the peak per-ELEMENT |B|: cube corner elements concentrate flux
    # under refinement (N=16 reaches ~2.4 T although the volume average peaks at 0.9 T).
    # The tables are exact lines, so widening the range changes nothing in-range -- it
    # makes the linear law genuinely defined where the C++ table lookup would otherwise
    # CLAMP (flat extrapolation) and the SolveHysteresis b_max guard rightly refuses.
    import numpy as np
    r = np.linspace(0.0, rmax, n)
    a0 = 1.0 / (MU0 * mu_r)
    tables = [(r.tolist(), (a0 * r).tolist()),
              (r.tolist(), (-irr1 * a0 * r).tolist()),
              (r.tolist(), (-irr2 * a0 * r).tolist())]
    return len(eta), np.asarray(eta, float), tables


def _drive():
    import numpy as np
    up0 = np.linspace(0.0, H0, N_RAMP + 1)[1:]
    down = np.linspace(H0, -H0, N_BRANCH + 1)[1:]
    up1 = np.linspace(-H0, H0, N_BRANCH + 1)[1:]
    hz = np.concatenate([up0, down, up1])
    h_steps = np.zeros((hz.size, 3))
    h_steps[:, 2] = hz
    return h_steps


def _peak_memory_mb() -> float:
    import psutil
    mem = psutil.Process(os.getpid()).memory_info()
    peak = getattr(mem, "peak_wset", None)
    return float(peak if peak is not None else mem.rss) / (1024.0 * 1024.0)


def run_case(n: int) -> dict:
    import numpy as np
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    import radia.vim as vim

    K, eta, tables = _synthetic_play()
    h_steps = _drive()
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n,
                                    mapping=lambda x, y, z: (x - 0.5, y - 0.5, z - 0.5))
        res = vim.SolveHysteresis(mesh, h_steps, play=(K, eta, tables))

    t_steps = [float(s["t_step_s"]) for s in res["steps"]]
    iters = [int(s["iters"]) for s in res["steps"]]
    hz = [float(s["h_applied"][2]) for s in res["steps"]]
    Hz = np.array([s["H_avg"][2] for s in res["steps"]])
    Bz = np.array([s["B_avg"][2] for s in res["steps"]])
    area = 0.5 * float(abs(np.sum(Hz[N_RAMP:] * np.roll(Bz[N_RAMP:], -1)
                                  - np.roll(Hz[N_RAMP:], -1) * Bz[N_RAMP:])))
    return dict(
        N=int(n), ndof=int(res["ndof"]), n_el=int(res["n_el"]), n_charge=int(res["n_charge"]),
        n_steps=len(t_steps),
        t_setup=float(res["t_setup_s"]),
        charge_gram_wall_s=float(res["charge_gram_wall_s"]),
        t_solve=float(res["t_steps_s"]),
        t_step_s=t_steps,
        t_step_mean_s=float(np.mean(t_steps)),
        t_step_max_s=float(np.max(t_steps)),
        iterations=int(sum(iters)),
        picard_iters=iters,
        cg_iters=[int(s["cg_iters"]) for s in res["steps"]],
        h_applied_z=hz,
        loop_area_J_per_m3=area,
        cpp_solve_timings={k: float(v) for k, v in res["cpp_solve_timings"].items()},
        converged=True,
        peak_memory_mb=_peak_memory_mb(),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sizes", default="8,12,16,20",
                    help="comma-separated hex cube sizes N (driver mode)")
    ap.add_argument("--case-n", type=int, default=None,
                    help="run ONE case in this process and print its JSON (subprocess mode)")
    ap.add_argument("--out", default=str(_HERE.parent / "results_bench_hysteresis_step.json"))
    args = ap.parse_args()

    if args.case_n is not None:
        print(json.dumps(run_case(args.case_n)))
        return 0

    sizes = [int(s) for s in args.sizes.split(",") if s.strip()]
    results = []
    for n in sizes:
        print(f"[bench_hysteresis_step] N={n} ...", flush=True)
        proc = subprocess.run([sys.executable, str(_HERE), "--case-n", str(n)],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"case N={n} failed:\n{proc.stdout}\n{proc.stderr}")
        case = json.loads(proc.stdout.strip().splitlines()[-1])
        results.append(case)
        print(f"  ndof={case['ndof']}  setup={case['t_setup']:.2f}s "
              f"(gram {case['charge_gram_wall_s']:.2f}s)  steps total={case['t_solve']:.2f}s  "
              f"mean/step={case['t_step_mean_s']*1e3:.1f}ms  iters={case['picard_iters']}",
              flush=True)

    data = dict(
        timestamp=_now(),
        hostname=platform.node(),
        benchmark="hysteresis_step_hex_cube",
        problem=dict(cube_size_m=1.0, H0_A_per_m=H0,
                     drive=f"virgin {N_RAMP} + descending {N_BRANCH} + ascending {N_BRANCH}",
                     material="synthetic play K=3 (mu_r 200, thresholds 0/0.3/0.7 T)",
                     mesh="structured hex NxNxN", order="RT1"),
        results=results,
    )
    Path(args.out).write_text(json.dumps(data, indent=1))
    print(f"[bench_hysteresis_step] wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
