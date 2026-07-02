"""bench_moment_solvers.py -- 3-way solver benchmark for the collocation-MMMM surface-charge demag:

  method 0 : dense LU            (direct; O(dof^2) memory, O(dof^3) factor)
  method 1 : BiCGSTAB, dense K   (matrix-free driver; chi-free geometry K built DENSE per solve)
  method 2 : HACApK-BiCGSTAB     (coarse tier: K as a RadHACApKMomentSystem H-matrix, O(N log N)
                                  matvec + storage, CROSS-SOLVE cached on radTApplication)

Measures, per (method, geometry, size) case in its OWN subprocess (Benchmark Policy: one case per
process so psutil peak_wset is per-case accurate):

  - t_setup (t_moment_system_build + t_moment_fieldgrad: K/L/diagK + precond build)
  - t_solve (t_linear_solve + t_lu_decomp)
  - t_wall_solve1 / t_wall_solve2: the SECOND solve on the same geometry shows the method-2
    cross-solve cache (K build skipped) vs methods 0/1 which rebuild -- the optimization-inner-loop
    per-iteration cost.
  - iterations, converged, peak_memory_mb, external-B probe values (cross-method correctness).

Geometries: 'cube' n x n x n (compact, near-field dominated -- H-matrix mostly dense blocks),
'bar' 2 x 2 x m (elongated -- admissible far blocks, ACA compression pays off), and 'ctype'
(the validated voxelized hex C-yoke of bench_multipole_moment_scaling.py: outer 0.12 m square,
inner cavity, +x gap; nxy=nside, nz=max(2,nside//3); in-plane drive H=[0,H0,0] -- the realistic
loop-heavy engineering geometry).

mdx-ready: imports `radia` DIRECTLY (PyPI on mdx / editable on LAB), fully self-contained.

  python bench_moment_solvers.py --sweep                      # full 3-method sweep -> JSON
  python bench_moment_solvers.py --method 2 --geom bar --nside 200   # one case, JSON to stdout

Requires radia >= the release carrying the method-2 moment HACApK route (rad.Solve(..., 2) on a
pure-hex moment object; 2026-07-02).
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

import radia as rad   # PyPI (mdx) or editable (LAB); NO sys.path hack -> mdx-clean

HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7 * np.pi
L = 0.01            # hex edge [m]
H0 = 1.0e3          # source field [A/m]


def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return (mem.peak_wset if hasattr(mem, "peak_wset") else mem.rss) / (1024.0 * 1024.0)


# --- C-yoke voxelization (inline copy of bench_multipole_moment_scaling.py's validated builder;
# kept inline so the benchmark has zero import-path dependency = mdx-clean) ---

def _inside_cyoke(cx, cy):
    # outer square [-0.06,0.06]^2, minus inner cavity [-0.035,0.035]^2, minus gap x>0.018 (C opens +x)
    if not (-0.06 <= cx <= 0.06 and -0.06 <= cy <= 0.06):
        return False
    if -0.035 <= cx <= 0.035 and -0.035 <= cy <= 0.035:
        return False
    if cx > 0.018 and -0.035 <= cy <= 0.035:
        return False
    return True


def ctype_grid(nside):
    """(nxy, nz) for a ctype case: nz keeps ~cubic voxels (z extent 0.04 vs xy extent 0.12)."""
    return nside, max(2, nside // 3)


def count_ctype_hexes(nside):
    """Exact hex count of the ctype voxelization (pure python, no radia -- used by the sweep planner)."""
    nxy, nz = ctype_grid(nside)
    xs = np.linspace(-0.06, 0.06, nxy + 1)
    n_inplane = sum(1 for j in range(nxy) for i in range(nxy)
                    if _inside_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])))
    return n_inplane * nz


def build_model(geom, nside, mu_r):
    """Pure-hex model: cube nside^3, bar 2x2xnside, or the voxelized C-yoke.
    Returns (container, nHex, extents, origin_centered)."""
    objs = []
    if geom in ("cube", "bar"):
        if geom == "cube":
            nx = ny = nz = nside
        else:
            nx = ny = 2
            nz = nside
        for iz in range(nz):
            for ix in range(nx):
                for iy in range(ny):
                    x0, y0, z0 = ix * L, iy * L, iz * L
                    v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                         [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                    h = rad.ObjHexahedron(v, [0, 0, 0])
                    rad.MatApl(h, rad.MatLin(mu_r))
                    objs.append(h)
        extents = (nx * L, ny * L, nz * L)
        centered = False
        bckg = [0.0, 0.0, MU0 * H0]                      # z-drive
    elif geom == "ctype":
        nxy, nz = ctype_grid(nside)
        xs = np.linspace(-0.06, 0.06, nxy + 1)
        zs = np.linspace(-0.02, 0.02, nz + 1)
        for k in range(nz):
            for j in range(nxy):
                for i in range(nxy):
                    if not _inside_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])):
                        continue
                    x0, x1, y0, y1 = xs[i], xs[i + 1], xs[j], xs[j + 1]
                    z0, z1 = zs[k], zs[k + 1]
                    v = [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                         [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]]
                    h = rad.ObjHexahedron(v, [0, 0, 0])
                    rad.MatApl(h, rad.MatLin(mu_r))
                    objs.append(h)
        extents = (0.12, 0.12, 0.04)
        centered = True
        bckg = [0.0, MU0 * H0, 0.0]                      # in-plane y-drive (excites the C loop)
    else:
        raise ValueError("geom must be cube|bar|ctype")
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: list(bckg))])
    return cont, len(objs), extents, centered


def probe_points(geom, extents):
    ex, ey, ez = extents
    if geom == "ctype":
        # gap-center (the engineering observable) + two external points; yoke is origin-centered
        return [[0.045, 0.0, 0.0],
                [0.09, 0.03, 0.005],
                [0.0, 0.09, 0.01]]
    return [[1.5 * ex, 0.25 * ey, 0.25 * ez],
            [0.25 * ex, 1.5 * ey, 0.50 * ez],
            [0.50 * ex, 0.50 * ey, 1.4 * ez]]


def run_case(args):
    """One (method, geom, nside) case in THIS process; JSON dict to stdout."""
    rad.UtiDelAll()
    rad.set_demag_backend("collocation_mmmm")
    rad.SolverConfig(bicgstab_tol=args.bicg_tol,
                     hacapk_eps=args.hacapk_eps, hacapk_leaf=args.hacapk_leaf,
                     hacapk_eta=args.hacapk_eta)
    if args.moment_kernel != "default":   # omit -> library default (analytic since 2026-07-02)
        rad.SolverConfig(moment_analytic_kernel=(args.moment_kernel == "analytic"))

    cont, n_hex, extents, _centered = build_model(args.geom, args.nside, args.mu_r)
    ndof = 6 * n_hex
    pts = probe_points(args.geom, extents)

    case = {
        "method": args.method, "geom": args.geom, "nside": args.nside,
        "n_hex": n_hex, "ndof": ndof, "mu_r": args.mu_r,
        "bicg_tol": args.bicg_tol, "hacapk_eps": args.hacapk_eps,
        "hacapk_leaf": args.hacapk_leaf, "hacapk_eta": args.hacapk_eta,
        "moment_kernel": args.moment_kernel,
    }

    def one_solve():
        t0 = time.perf_counter()
        ret = rad.Solve(cont, 1e-8, 3000, args.method)
        t1 = time.perf_counter()
        st = dict(rad.GetSolveStats())
        return t1 - t0, ret, st

    # solve 1: cold (all builds included)
    t_wall1, ret1, st1 = one_solve()
    # solve 2: same geometry again = the optimization-inner-loop per-iteration cost.
    # method 2 hits the cross-solve K/L/diagK cache; methods 0/1 rebuild.
    t_wall2, ret2, st2 = one_solve()

    B = np.asarray([rad.Fld(cont, "b", p) for p in pts], float)

    def split(st):
        t_setup = float(st.get("t_moment_system_build", 0.0)) + float(st.get("t_moment_fieldgrad", 0.0))
        t_solve = float(st.get("t_linear_solve", 0.0)) + float(st.get("t_lu_decomp", 0.0))
        return t_setup, t_solve

    t_setup1, t_solve1 = split(st1)
    t_setup2, t_solve2 = split(st2)
    finite = bool(np.all(np.isfinite(B)) and np.linalg.norm(B) > 0.0)

    case.update({
        "t_wall_solve1": t_wall1, "t_wall_solve2": t_wall2,
        "t_setup": t_setup1, "t_solve": t_solve1,          # Benchmark Policy names (cold solve)
        "t_setup_solve2": t_setup2, "t_solve_solve2": t_solve2,
        "iterations": int(st1.get("linear_iterations", 0)),
        "iterations_solve2": int(st2.get("linear_iterations", 0)),
        "nonl_iterations": int(st1.get("nonl_iterations", 0)),
        "converged": bool(ret1 is not None and ret2 is not None and finite),
        "solve_return_1": ret1, "solve_return_2": ret2,
        "num_threads": int(st1.get("num_threads", 0)),
        "B_probe": B.tolist(),
        "peak_memory_mb": get_peak_memory_mb(),
        "stats_solve1": {k: (float(v) if isinstance(v, (int, float)) else bool(v)) for k, v in st1.items()},
        "stats_solve2": {k: (float(v) if isinstance(v, (int, float)) else bool(v)) for k, v in st2.items()},
    })
    rad.UtiDelAll()
    print(json.dumps(case))
    return case


# ---------------------------------------------------------------- sweep driver

def sweep(args):
    """Spawn one subprocess per case (Benchmark Policy: per-case memory accuracy), aggregate JSON."""
    cube_sides = args.cube_sides if args.cube_sides is not None else [4, 6, 8, 10, 12, 16, 20]
    bar_sides = args.bar_sides if args.bar_sides is not None else [25, 50, 100, 200, 400, 800]
    ctype_sides = args.ctype_sides if args.ctype_sides is not None else [12, 18, 24, 30, 36]
    plan = []
    for geom, sides in (("cube", cube_sides), ("bar", bar_sides), ("ctype", ctype_sides)):
        for ns in sides:
            if geom == "cube":
                n_hex = ns ** 3
            elif geom == "bar":
                n_hex = 4 * ns
            else:
                n_hex = count_ctype_hexes(ns)
            ndof = 6 * n_hex
            for method in (0, 1, 2):
                if method == 0 and ndof > args.max_lu_dof:
                    continue          # dense LU wall
                if method == 1 and ndof > args.max_dense_dof:
                    continue          # dense-K storage wall
                plan.append((method, geom, ns, ndof))

    results = []
    for i, (method, geom, ns, ndof) in enumerate(plan):
        cmd = [sys.executable, os.path.abspath(__file__),
               "--method", str(method), "--geom", geom, "--nside", str(ns),
               "--mu-r", str(args.mu_r), "--bicg-tol", str(args.bicg_tol),
               "--hacapk-eps", str(args.hacapk_eps), "--hacapk-leaf", str(args.hacapk_leaf),
               "--hacapk-eta", str(args.hacapk_eta)]
        if args.moment_kernel != "default":
            cmd += ["--moment-kernel", args.moment_kernel]
        print(f"[{i+1}/{len(plan)}] method={method} {geom} nside={ns} (ndof={ndof}) ...", flush=True)
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=args.case_timeout)
        dt = time.perf_counter() - t0
        if proc.returncode != 0:
            print(f"    CASE FAILED (exit {proc.returncode}, {dt:.1f} s): {proc.stderr.strip()[-400:]}", flush=True)
            results.append({"method": method, "geom": geom, "nside": ns, "ndof": ndof,
                            "converged": False, "error": f"exit {proc.returncode}",
                            "stderr_tail": proc.stderr.strip()[-400:]})
            continue
        case = json.loads(proc.stdout.strip().splitlines()[-1])
        print(f"    ok {dt:6.1f} s  wall1={case['t_wall_solve1']:8.3f}  wall2={case['t_wall_solve2']:8.3f}"
              f"  iters={case['iterations']:4d}  mem={case['peak_memory_mb']:7.1f} MB", flush=True)
        results.append(case)

    # cross-method correctness: rel |B - B_ref| per (geom, nside); ref = lowest available method
    by_key = {}
    for c in results:
        if c.get("converged"):
            by_key.setdefault((c["geom"], c["nside"]), {})[c["method"]] = np.asarray(c["B_probe"], float)
    for c in results:
        key = (c.get("geom"), c.get("nside"))
        if not c.get("converged") or key not in by_key:
            continue
        refs = by_key[key]
        ref_method = min(refs.keys())
        if c["method"] != ref_method:
            Bref = refs[ref_method]
            Bc = np.asarray(c["B_probe"], float)
            c["rel_vs_ref"] = float(np.linalg.norm(Bc - Bref) / max(np.linalg.norm(Bref), 1e-300))
            c["ref_method"] = ref_method

    data = {
        "timestamp": datetime.now().isoformat(),
        "hostname": platform.node(),
        "benchmark": "moment_solvers_lu_bicgstab_hacapk",
        "problem": {
            "hex_edge_m": L, "H0_A_per_m": H0, "mu_r": args.mu_r,
            "bicg_tol": args.bicg_tol, "hacapk_eps": args.hacapk_eps,
            "hacapk_leaf": args.hacapk_leaf, "hacapk_eta": args.hacapk_eta,
            "cube_sides": cube_sides, "bar_sides": bar_sides, "ctype_sides": ctype_sides,
            "max_lu_dof": args.max_lu_dof, "max_dense_dof": args.max_dense_dof,
            "radia_version": getattr(rad, "__version__", "unknown"),
            "python_version": platform.python_version(),
        },
        "results": results,
    }
    out = args.json_out or os.path.join(HERE, "results_moment_solvers.json")
    with open(out, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Results saved to {out}")

    # console summary table
    print(f"\n{'geom':5s} {'nside':>5s} {'ndof':>7s} | {'m0 wall1':>9s} {'m1 wall1':>9s} {'m2 wall1':>9s} |"
          f" {'m0 wall2':>9s} {'m1 wall2':>9s} {'m2 wall2':>9s} | {'m2 rel':>9s}")
    for key in sorted(by_key.keys(), key=lambda k: (k[0], k[1])):
        row = {c["method"]: c for c in results if (c.get("geom"), c.get("nside")) == key and c.get("converged")}
        def w(m, f):
            return f"{row[m][f]:9.3f}" if m in row else "        -"
        rel2 = f"{row[2]['rel_vs_ref']:9.2e}" if 2 in row and "rel_vs_ref" in row[2] else "        -"
        if key[0] == "cube":
            ndof = 6 * key[1] ** 3
        elif key[0] == "bar":
            ndof = 6 * 4 * key[1]
        else:
            ndof = 6 * count_ctype_hexes(key[1])
        print(f"{key[0]:5s} {key[1]:5d} {ndof:7d} | {w(0,'t_wall_solve1')} {w(1,'t_wall_solve1')} {w(2,'t_wall_solve1')} |"
              f" {w(0,'t_wall_solve2')} {w(1,'t_wall_solve2')} {w(2,'t_wall_solve2')} | {rel2}")


def build_argparser():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sweep", action="store_true", help="run the full 3-method sweep (subprocess per case)")
    ap.add_argument("--method", type=int, choices=(0, 1, 2), default=2, help="0=LU, 1=BiCGSTAB dense-K, 2=HACApK-BiCGSTAB")
    ap.add_argument("--geom", choices=("cube", "bar", "ctype"), default="cube")
    ap.add_argument("--nside", type=int, default=6,
                    help="cube: n^3 hex; bar: 2x2xn hex; ctype: nxy=n voxel C-yoke (nz=max(2,n//3))")
    ap.add_argument("--mu-r", type=float, default=200.0)
    ap.add_argument("--bicg-tol", type=float, default=1e-8)
    ap.add_argument("--hacapk-eps", type=float, default=1e-4)
    ap.add_argument("--hacapk-leaf", type=int, default=32)
    ap.add_argument("--hacapk-eta", type=float, default=2.0)
    ap.add_argument("--moment-kernel", choices=("default", "analytic", "gauss"), default="default",
                    help="face kernel: 'default' follows the library (analytic since 2026-07-02); "
                         "'gauss' forces the 64-sample quadrature for cross-checks")
    ap.add_argument("--cube-sides", type=int, nargs="*", default=None)
    ap.add_argument("--bar-sides", type=int, nargs="*", default=None)
    ap.add_argument("--ctype-sides", type=int, nargs="*", default=None)
    ap.add_argument("--max-lu-dof", type=int, default=8000, help="skip method 0 above this dof (LU wall)")
    ap.add_argument("--max-dense-dof", type=int, default=26000, help="skip method 1 above this dof (dense-K wall)")
    ap.add_argument("--case-timeout", type=float, default=3600.0, help="per-case subprocess timeout [s]")
    ap.add_argument("--json-out", default=None)
    return ap


if __name__ == "__main__":
    a = build_argparser().parse_args()
    if a.sweep:
        sweep(a)
    else:
        run_case(a)
