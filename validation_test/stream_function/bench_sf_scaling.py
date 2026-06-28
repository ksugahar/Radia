#!/usr/bin/env python3
"""SF design scaling benchmark: build time + peak memory vs surface-mesh DOF.

Measures ``calc_streamfunction._build_problem`` (mesh load + Biot-Savart design
matrix assembly + sparse seminorm fold + ACA+TSVD factorisation) on a cylinder
lateral-surface Gx coil at a density sweep, plus the dominant sub-cost
(``_assemble_biot_savart``).  One case at a time (Benchmark Policy); records
peak RSS via psutil.

Motivation (2026-06-03): the regularisation seminorm and the DOF-reduction
matrix R were dense N x N, making ``_build_problem`` an O(N^2) time + memory
wall (13.3 s and 3.75 GB at N = 15 260 DOF).  Sparsifying both (the FE
stiffness is ~2 % dense) + factor-once ``splu`` in
``RegularizedTSVD.from_stiffness`` removed it; the ACA fold is now ~0.3 s and
peak memory ~0.2 GB at the same DOF count (44x / 19x).  This script locks the
post-fix curve; results JSON is committed next to it (Data Persistence Policy).

    python bench_sf_scaling.py                 # default sweep, writes JSON + PNG
    python bench_sf_scaling.py --max-ndof 30000
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import time

import numpy as np


def _peak_mb():
    try:
        import psutil
        m = psutil.Process(os.getpid()).memory_info()
        return float(getattr(m, "peak_wset", m.rss)) / 2 ** 20
    except Exception:
        return -1.0


def _gen_cylinder(maxh, out):
    from netgen.occ import Cylinder, Pnt, Z, OCCGeometry
    a, L = 0.15, 0.50
    cyl = Cylinder(Pnt(0, 0, -L / 2), Z, r=a, h=L)
    lateral = max(cyl.faces, key=lambda f: f.mass)
    OCCGeometry(lateral).GenerateMesh(maxh=maxh).Save(out)


def _gen_eval(out):
    from netgen.occ import Sphere, Pnt, OCCGeometry
    OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(maxh=0.025).Save(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--work-dir", default=os.path.join(os.environ.get("TEMP", "/tmp"),
                                                       "sf_bench_scale"))
    ap.add_argument("--max-ndof", type=float, default=16000.0,
                    help="stop refining once ndof exceeds this")
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--eval-max", type=int, default=80)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "src", "radia"))
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "src", "radia", "panels"))
    import calc_streamfunction as C
    from ngsolve import Mesh, H1, specialcf, TaskManager

    os.makedirs(args.work_dir, exist_ok=True)
    eval_vol = os.path.join(args.work_dir, "eval_dsv.vol")
    _gen_eval(eval_vol)

    rows = []
    with TaskManager():
        # warmup: first build pays NGSolve/MKL/ACA-module cold-start (~7 s);
        # discard it so the recorded curve is steady-state, not init-dominated.
        _wu = os.path.join(args.work_dir, "cyl_warmup.vol")
        _gen_cylinder(0.05, _wu)
        C._build_problem(argparse.Namespace(
            coil_vol=_wu, eval_vol=eval_vol, target_cf="x", order=args.order,
            regularize="h1", confine="abe", alpha=0.0,
            eval_max=args.eval_max, method="design"))
        for maxh in [0.04, 0.022, 0.014, 0.009, 0.006, 0.0045, 0.0035]:
            coil_vol = os.path.join(args.work_dir, f"cyl_h{int(maxh*1e4):04d}.vol")
            _gen_cylinder(maxh, coil_vol)
            ns = argparse.Namespace(
                coil_vol=coil_vol, eval_vol=eval_vol, target_cf="x",
                order=args.order, regularize="h1", confine="abe",
                alpha=0.0, eval_max=args.eval_max, method="design")
            # dominant sub-cost: the Biot-Savart design-matrix assembly
            coil = Mesh(coil_vol)
            fes = H1(coil, order=args.order, definedon=coil.Boundaries(".*"))
            nrm = specialcf.normal(3)
            evalm = Mesh(eval_vol)
            cpts = np.array([list(p.point) for p in evalm.vertices])
            cpts = cpts[np.linspace(0, len(cpts) - 1,
                                    min(len(cpts), args.eval_max)).astype(int)]
            t0 = time.perf_counter()
            C._assemble_biot_savart(fes, nrm, cpts, [2])
            t_assemble = time.perf_counter() - t0
            # full build (load + assemble + sparse fold + ACA)
            t0 = time.perf_counter()
            P = C._build_problem(ns)
            t_build = time.perf_counter() - t0
            row = {"maxh": maxh, "nv": int(coil.nv), "ndof": int(P["fes"].ndof),
                   "n_free": int(P["n_free"]), "M": int(P["n_constraint"]),
                   "t_assemble_s": round(t_assemble, 4),
                   "t_build_s": round(t_build, 4),
                   "peak_memory_mb": round(_peak_mb(), 1)}
            rows.append(row)
            print(f"ndof={row['ndof']:6d} M={row['M']:3d}  "
                  f"assemble={row['t_assemble_s']*1e3:7.0f}ms  "
                  f"build={row['t_build_s']*1e3:7.0f}ms  "
                  f"peak={row['peak_memory_mb']:.0f}MB")
            if P["fes"].ndof > args.max_ndof:
                break

    from datetime import datetime
    data = {"timestamp": datetime.now().isoformat(timespec="seconds"),
            "hostname": platform.node(),
            "benchmark": "sf_design_scaling",
            "problem": {"geometry": "cylinder lateral surface r=0.15 L=0.5",
                        "target": "Gx (Bz=x)", "order": args.order,
                        "regularize": "h1", "confine": "abe"},
            "note": "post-sparse-seminorm fix (2026-06-03); t_build is "
                    "mesh-load + Biot-Savart assemble + sparse R^T S R fold + "
                    "ACA+TSVD factorise; the dense seminorm wall (13.3s/3.75GB "
                    "at ndof=15260) is gone -> ACA fold ~0.3s, peak ~0.2GB.",
            "results": rows}
    jpath = os.path.join(args.out_dir, "bench_sf_scaling.json")
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"\nResults -> {jpath}")

    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            # Lab figure convention (radia_mcp.figure): authored AT the 16 cm
            # embed width, 10 pt, Times New Roman, TrueType; NO in-figure title
            # (the docs / LaTeX caption carry it -- save_lab_figure RAISES on a
            # baked-in title).  examples/ may depend on radia-mcp; the shipped
            # panel calc cannot, so its diagnostic plots apply the same
            # conventions inline instead.
            from radia_mcp.figure import lab_figure, save_lab_figure
            nd = [r["ndof"] for r in rows]
            fig, (ax1, ax2) = lab_figure(embed_width_cm=16.0, aspect=0.42, ncols=2)
            ax1.loglog(nd, [r["t_build_s"] for r in rows], "o-",
                       label="build (total)")
            ax1.loglog(nd, [r["t_assemble_s"] for r in rows], "s--",
                       label="assemble A")
            ax1.set_xlabel("surface FE DOF")
            ax1.set_ylabel("build time (s)")
            ax1.legend()
            ax1.grid(True, which="both", alpha=0.3)
            ax2.loglog(nd, [r["peak_memory_mb"] for r in rows], "o-", color="C3")
            ax2.set_xlabel("surface FE DOF")
            ax2.set_ylabel("peak memory (MB)")
            ax2.grid(True, which="both", alpha=0.3)
            info = save_lab_figure(
                fig, os.path.join(args.out_dir, "bench_sf_scaling"),
                embed_width_cm=16.0)
            print(f"Plot    -> {', '.join(info['wrote'])}")
        except Exception as e:
            print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
