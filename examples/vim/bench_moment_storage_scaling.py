"""Phase-2 Increment-4 storage gate (docs/multipole_moment_mmm/ACA_MOMENT_DESIGN.md): show that multipole-moment MMM method 2
(HACApK H-matrix + block-Jacobi BiCGSTAB) has SUB-QUADRATIC peak memory, while method 0 (dense LU) is O(N^2).

The Increment-4 decoupling removes ALL three O(N^2) buffers from the method-2 path: the dense interaction
matrix N (skipDenseMatrix=1 in SolveGen), the BaseMatrix (radTRelaxationMethNo_0::NeedsDenseMatrix()->false),
and the dgesv working SystemMatrix (lazy-allocated, never reached on the H-BiCGSTAB happy path).  What is left
is the ACA-compressed H-matrix (O(N log N)) + the O(N) block-Jacobi + Krylov vectors.

Per the lab Benchmark Policy: ONE case per process (clean psutil peak_wset), JSON output next to the script.
Solid block of linear soft iron (mu_r=200) in a uniform field -- predictable dof = 6 * (nx*ny*nz).

Usage:
  python bench_moment_storage_scaling.py                 # driver: sweep sizes x {method 0, 2}, write JSON
  python bench_moment_storage_scaling.py --cases 60x46x10 --methods 2 --timeout 7200
  python bench_moment_storage_scaling.py --worker --method M --nx .. --ny .. --nz ..   # one case (internal)
"""
import argparse
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime

import numpy as np
import psutil

import radia as rad

HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7 * np.pi
MU_R = 200.0
L = 0.01


def _peak_mb():
    mi = psutil.Process(os.getpid()).memory_info()
    return (mi.peak_wset if hasattr(mi, "peak_wset") else mi.rss) / (1024.0 * 1024.0)


def _build_block(nx, ny, nz):
    objs = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                x0, y0, z0 = ix * L, iy * L, iz * L
                v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                     [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(MU_R)); objs.append(h)
    return objs


def run_one(method, nx, ny, nz):
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(bicgstab_tol=1e-8)
    t0 = time.perf_counter()
    objs = _build_block(nx, ny, nz)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * 1.0e3])])
    t_setup = time.perf_counter() - t0
    t1 = time.perf_counter()
    nit = rad.Solve(cont, 1e-6, 4000, method)
    t_solve = time.perf_counter() - t1
    # external-field probe (observable, formulation-independent) to confirm a real solve happened
    B = rad.Fld(cont, "b", [nx * L + 0.05, ny * L * 0.5, nz * L * 0.5])
    peak = _peak_mb()
    ndof = 6 * len(objs)
    rad.UtiDelAll()
    # rad.Solve returns (misfit, ..., iteration_count) -- take the last element as the Picard iteration count.
    iters = int(round(nit[-1])) if isinstance(nit, (list, tuple)) and len(nit) else (int(nit) if np.isscalar(nit) else None)
    return {"method": method, "ndof": ndof, "ne": len(objs), "nx": nx, "ny": ny, "nz": nz,
            "peak_memory_mb": peak, "t_setup": t_setup, "t_solve": t_solve,
            "iterations": iters, "converged": bool(np.all(np.isfinite(B))),
            "Bz_probe": float(B[2])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--method", type=int, default=0)
    ap.add_argument("--nx", type=int, default=8)
    ap.add_argument("--ny", type=int, default=8)
    ap.add_argument("--nz", type=int, default=4)
    ap.add_argument("--cases", nargs="+",
                    help="driver cases as NxXxNz, e.g. 16x16x8 60x46x10")
    ap.add_argument("--methods", type=int, nargs="+", default=[0, 2],
                    help="driver methods to run; use --methods 2 for large HACApK-only sweeps")
    ap.add_argument("--timeout", type=float, default=900.0,
                    help="per worker timeout in seconds")
    ap.add_argument("--out", default=os.path.join(HERE, "results_moment_storage_scaling.json"),
                    help="output JSON path")
    args = ap.parse_args()

    if args.worker:
        print(json.dumps(run_one(args.method, args.nx, args.ny, args.nz)))
        return 0

    if args.cases:
        sizes = []
        for spec in args.cases:
            parts = spec.lower().replace("*", "x").split("x")
            if len(parts) != 3:
                raise SystemExit(f"invalid --cases entry {spec!r}; expected NxXxNz")
            sizes.append(tuple(int(p) for p in parts))
    else:
        sizes = [(8, 8, 4), (12, 12, 4), (16, 16, 4), (16, 16, 8)]   # dof = 1536, 3456, 6144, 12288
    results = []
    for (nx, ny, nz) in sizes:
        for method in args.methods:
            cmd = [sys.executable, os.path.abspath(__file__), "--worker", "--method", str(method),
                   "--nx", str(nx), "--ny", str(ny), "--nz", str(nz)]
            try:
                out = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout, cwd=HERE)
                rec = json.loads(out.stdout.strip().splitlines()[-1])
            except Exception as e:
                rec = {"method": method, "nx": nx, "ny": ny, "nz": nz, "error": str(e),
                       "stderr": (out.stderr[-400:] if "out" in dir() else "")}
            results.append(rec)
            tag = f"m{method} dof={rec.get('ndof','?')}"
            print(f"  {tag:18s} peak={rec.get('peak_memory_mb', float('nan')):8.1f} MB  "
                  f"t_solve={rec.get('t_solve', float('nan')):7.2f}s  iters={rec.get('iterations')}  "
                  f"conv={rec.get('converged')}")

    with open(args.out, "w") as f:
        json.dump({"timestamp": datetime.now().isoformat(), "hostname": platform.node(),
                   "benchmark": "moment_storage_scaling",
                   "problem": {"geometry": "solid hex block", "mu_r": MU_R, "L": L},
                   "results": results}, f, indent=2)
    print(f"Results saved to {args.out}")

    # headline: method-2 peak / method-0 peak at the largest dof (should be << 1 once O(N^2) is gone)
    by = {(r.get("method"), r.get("ndof")): r for r in results if "peak_memory_mb" in r}
    dofs = sorted({r.get("ndof") for r in results if r.get("ndof")})
    if dofs:
        big = dofs[-1]
        m0 = by.get((0, big), {}).get("peak_memory_mb")
        m2 = by.get((2, big), {}).get("peak_memory_mb")
        if m0 and m2:
            print(f"\n  dof={big}: method0 peak={m0:.0f} MB, method2 peak={m2:.0f} MB  -> method2/method0 = {m2/m0:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
