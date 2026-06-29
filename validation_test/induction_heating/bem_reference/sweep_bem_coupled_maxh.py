"""wp_maxh sweep for CoupledBEMSolver Delta_L convergence study.

Goal: decide whether Fix A (nodal grad_s J_wp + multi-point A_wp) stabilises
the back-reaction Delta_L. Pass criterion: |Delta_L| variation < 10% across
the sweep. Fail: abandon Delta_L from PEEC+BEM; stick with PEEC+FEM-Kelvin.

Reuses ``experiment_coupled_bem`` geometry (gapped torus + cylinder wp, Cu,
1 kHz). No new physics.

Usage:
    python sweep_bem_coupled_maxh.py                      # default sweep
    python sweep_bem_coupled_maxh.py --maxh 0.004 0.002   # explicit
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src", "radia"))
sys.path.insert(0, os.path.join(REPO, "src", "radia", "panels"))
sys.path.insert(0, os.path.join(REPO, "src"))
# Sibling directory for experiment_coupled_bem imports (impedance_esim etc.)
sys.path.insert(0, os.path.join(REPO, "validation_test", "induction_heating", "cubit_panels_legacy"))
sys.path.insert(0, HERE)

MU_0 = 4e-7 * math.pi

R_COIL, A_COIL, GAP_DEG = 0.030, 0.003, 5
R_WP, H_WP = 0.010, 0.020
SIGMA = 5.8e7
FREQ = 1000.0


def run_one(maxh_wp):
    from ngsolve import Mesh, BND
    from ngsolve import TaskManager
    from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry, Glue
    from impedance_esim import make_gapped_torus_mesh
    from radia.bem_coupled_solver import CoupledBEMSolver

    omega = 2 * math.pi * FREQ
    delta = math.sqrt(2.0 / (omega * MU_0 * SIGMA))
    Z_s = complex(1, 1) / (SIGMA * delta)

    # Coil mesh is fixed (sweep only wp)
    mesh_coil = make_gapped_torus_mesh(R_COIL, A_COIL, gap_deg=GAP_DEG,
                                        maxh=A_COIL)

    wp_cyl = Cylinder(Pnt(0, 0, -H_WP / 2), Dir(0, 0, 1), R_WP, H_WP)
    for f in wp_cyl.faces:
        f.name = "wp"
    mesh_wp = Mesh(OCCGeometry(Glue(wp_cyl.faces)).GenerateMesh(maxh=maxh_wp))
    with TaskManager():
        mesh_wp.Curve(3)

        n_wp = mesh_wp.GetNE(BND)
        n_coil = mesh_coil.GetNE(BND)

        t0 = time.perf_counter()
        solver = CoupledBEMSolver(mesh_coil, mesh_wp)
        t_asm = time.perf_counter() - t0

        t0 = time.perf_counter()
        r = solver.solve(Z_s=Z_s, omega=omega, max_iter=10, tol=1e-3,
                         verbose=False)
        t_solve = time.perf_counter() - t0

        return {
            "maxh_wp": maxh_wp,
            "n_wp_elems": n_wp,
            "n_coil_elems": n_coil,
            "L_air_nH": r["L_air"] * 1e9,
            "L_total_nH": r["L_total"] * 1e9,
            "Delta_L_nH": r["Delta_L"] * 1e9,
            "P_total_W": r["P_total"],
            "iterations": r["iterations"],
            "t_assembly_s": round(t_asm, 1),
            "t_solve_s": round(t_solve, 1),
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxh", type=float, nargs="+",
                        default=[0.004, 0.003, 0.002, 0.0015],
                        help="wp_maxh values to sweep [m]")
    parser.add_argument("--tag", default="baseline",
                        help="Label for this sweep (e.g. 'baseline', 'fixA')")
    parser.add_argument("--output", default=os.path.join(HERE, "sweep_results.json"))
    args = parser.parse_args()

    results = []
    for maxh in args.maxh:
        print(f"\n=== wp_maxh = {maxh*1e3:.2f} mm ===")
        try:
            r = run_one(maxh)
        except Exception as e:
            print(f"  FAILED: {e}")
            r = {"maxh_wp": maxh, "error": str(e)}
        results.append(r)
        if "error" not in r:
            print(f"  n_wp={r['n_wp_elems']}  L_air={r['L_air_nH']:.3f} nH")
            print(f"  Delta_L={r['Delta_L_nH']:+.3f} nH  L_total={r['L_total_nH']:.3f} nH")
            print(f"  P={r['P_total_W']:.3e} W  iters={r['iterations']}  "
                  f"t_asm={r['t_assembly_s']}s  t_solve={r['t_solve_s']}s")

    # Append to output JSON (keep history of runs)
    try:
        with open(args.output) as f:
            existing = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        existing = []
    existing.append({"tag": args.tag, "sweep": results})
    with open(args.output, "w") as f:
        json.dump(existing, f, indent=2)

    print()
    print(f"=== Summary (tag='{args.tag}') ===")
    good = [r for r in results if "error" not in r]
    if len(good) >= 2:
        dLs = [r["Delta_L_nH"] for r in good]
        spread = (max(dLs) - min(dLs)) / abs(np.mean(dLs)) * 100
        print(f"{'maxh[mm]':>10s} {'n_wp':>8s} {'Delta_L[nH]':>13s}")
        for r in good:
            print(f"{r['maxh_wp']*1e3:>10.3f} {r['n_wp_elems']:>8d} "
                  f"{r['Delta_L_nH']:>13.3f}")
        print(f"\nDelta_L spread (max-min)/mean = {spread:.1f}%  "
              f"(pass criterion: < 10%)")
    print(f"Results appended to {args.output}")


if __name__ == "__main__":
    main()
