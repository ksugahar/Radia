"""Benchmark: Compact AMS + COCR  vs  ABMC-ICCG  vs  ABMC-COCR.

Compares auxiliary-space preconditioning (AMS) against incomplete Cholesky (IC)
with ABMC parallel ordering for complex eddy current problems.

EMD paper reference (Hiruma, SA-26-001, 3,670,328 DOFs, 30 kHz):
  IC only:              15,838 iter, 5964.8s  (376.7 ms/iter)
  EMD (IC+AMG V-cycle):  4,069 iter, 1716.8s  (422.0 ms/iter)
  EMD (IC+AMG W-cycle):  2,935 iter, 1552.9s  (529.1 ms/iter)
  EMD (IC+GenEO 24dom):  1,004 iter,  550.8s  (548.6 ms/iter)

Usage:
    python bench_ams_vs_abmc.py [mesh_name]
    python bench_ams_vs_abmc.py mesh1_3.5T
    python bench_ams_vs_abmc.py --all
"""

import json
import os
import platform
import psutil
import sys
import time
from datetime import datetime

import numpy as np
from netgen.read_gmsh import ReadGmsh
from ngsolve import *
import sparsesolv_ngsolve as ssn

mu0 = 4e-7 * np.pi
freq = 30e3
omega = 2 * np.pi * freq
eps = 1e-6
sigma_cu = 5.96e7
mu_r_core = 1000

maxiter = 20000
tol = 1e-10

ALL_MESHES = ["mesh1_2.5T", "mesh1_3.5T", "mesh1_4.5T", "mesh1_5.5T"]


def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return mem.peak_wset / (1024 * 1024) if hasattr(mem, 'peak_wset') else mem.rss / (1024 * 1024)


def setup_problem(mesh_name):
    mesh_file = os.path.join(os.path.dirname(__file__), mesh_name)
    if not mesh_file.endswith(".msh"):
        mesh_file += ".msh"
    mesh = Mesh(ReadGmsh(mesh_file))

    nu_cf = 1.0 / (mu0 * IfPos(mesh.MaterialCF({"core": 1}), mu_r_core, 1.0))
    sigma_cf = mesh.MaterialCF({"cond": sigma_cu}, default=0)

    fes = HCurl(mesh, order=1, nograds=True,
                dirichlet="dirichlet", complex=True)
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx
    a += eps * nu_cf * u * v * dx
    a += 1j * omega * sigma_cf * u * v * dx("cond")
    a.Assemble()

    f = LinearForm(fes)
    f += nu_cf * CF((0, 0, 1)) * v * dx("cond")
    f.Assemble()

    fes_real = HCurl(mesh, order=1, nograds=True,
                     dirichlet="dirichlet", complex=False)
    u_r, v_r = fes_real.TnT()
    a_real = BilinearForm(fes_real)
    a_real += nu_cf * curl(u_r) * curl(v_r) * dx
    a_real += eps * nu_cf * u_r * v_r * dx
    a_real += abs(omega) * sigma_cf * u_r * v_r * dx("cond")
    a_real.Assemble()

    G_mat, h1_fes = fes_real.CreateGradient()

    nv = mesh.nv
    coord_x = [0.0] * nv
    coord_y = [0.0] * nv
    coord_z = [0.0] * nv
    for i in range(nv):
        p = mesh.ngmesh.Points()[i + 1]
        coord_x[i] = p[0]
        coord_y[i] = p[1]
        coord_z[i] = p[2]

    # RHS norm
    fd = fes.FreeDofs()
    rhs_free = f.vec.CreateVector()
    rhs_free.data = f.vec
    for i in range(fes.ndof):
        if not fd[i]:
            rhs_free[i] = 0
    rhs_norm = float(sqrt(abs(InnerProduct(rhs_free, rhs_free))))

    return {
        "mesh": mesh, "fes": fes, "a": a, "f": f,
        "a_real": a_real, "fes_real": fes_real,
        "G_mat": G_mat, "h1_fes": h1_fes,
        "coord_x": coord_x, "coord_y": coord_y, "coord_z": coord_z,
        "rhs_norm": rhs_norm, "ne": mesh.ne, "nv": mesh.nv,
    }


def true_residual(p, gfu):
    r = p["f"].vec.CreateVector()
    r.data = p["f"].vec - p["a"].mat * gfu.vec
    fd = p["fes"].FreeDofs()
    for i in range(p["fes"].ndof):
        if not fd[i]:
            r[i] = 0
    return float(sqrt(abs(InnerProduct(r, r)))) / p["rhs_norm"]


def run_compact_ams_cocr(p):
    """Compact AMS + COCR solver."""
    ndof = p["fes"].ndof

    t0 = time.perf_counter()
    pre = ssn.ComplexCompactAMSPreconditioner(
        a_real_mat=p["a_real"].mat, grad_mat=p["G_mat"],
        freedofs=p["fes_real"].FreeDofs(),
        coord_x=p["coord_x"], coord_y=p["coord_y"], coord_z=p["coord_z"],
        ndof_complex=ndof, cycle_type=1, print_level=0)
    t_setup = time.perf_counter() - t0

    gfu = GridFunction(p["fes"])
    with TaskManager():
        t0 = time.perf_counter()
        inv = ssn.COCRSolver(
            p["a"].mat, pre, freedofs=p["fes"].FreeDofs(),
            maxiter=maxiter, tol=tol, printrates=False)
        gfu.vec.data = inv * p["f"].vec
        t_solve = time.perf_counter() - t0

    iters = inv.iterations
    res = true_residual(p, gfu)
    return {
        "method": "Compact AMS + COCR",
        "iterations": iters,
        "t_setup": round(t_setup, 3),
        "t_solve": round(t_solve, 3),
        "t_total": round(t_setup + t_solve, 3),
        "ms_per_iter": round(t_solve / max(iters, 1) * 1000, 1),
        "true_residual": float(res),
        "converged": iters < maxiter,
    }


def run_abmc_iccg(p):
    """ABMC-ordered ICCG solver (IC preconditioner)."""
    t0 = time.perf_counter()
    solver = ssn.SparseSolvSolver(
        p["a"].mat, method="ICCG",
        freedofs=p["fes"].FreeDofs(),
        tol=tol, maxiter=maxiter, shift=1.05,
        save_best_result=True, printrates=False,
        use_abmc=True, abmc_block_size=4, abmc_num_colors=4)
    solver.auto_shift = True
    t_setup = time.perf_counter() - t0

    gfu = GridFunction(p["fes"])
    t0 = time.perf_counter()
    gfu.vec.data = solver * p["f"].vec
    t_solve = time.perf_counter() - t0

    result = solver.last_result
    iters = result.iterations
    res = true_residual(p, gfu)
    return {
        "method": "ABMC-ICCG (IC, block=4, colors=4)",
        "iterations": iters,
        "t_setup": round(t_setup, 3),
        "t_solve": round(t_solve, 3),
        "t_total": round(t_setup + t_solve, 3),
        "ms_per_iter": round(t_solve / max(iters, 1) * 1000, 1),
        "true_residual": float(res),
        "converged": result.converged,
    }


def run_abmc_cocr(p):
    """ABMC-ordered COCR solver (IC preconditioner)."""
    t0 = time.perf_counter()
    solver = ssn.SparseSolvSolver(
        p["a"].mat, method="COCR",
        freedofs=p["fes"].FreeDofs(),
        tol=tol, maxiter=maxiter, shift=1.05,
        save_best_result=True, printrates=False,
        use_abmc=True, abmc_block_size=4, abmc_num_colors=4)
    solver.auto_shift = True
    t_setup = time.perf_counter() - t0

    gfu = GridFunction(p["fes"])
    t0 = time.perf_counter()
    gfu.vec.data = solver * p["f"].vec
    t_solve = time.perf_counter() - t0

    result = solver.last_result
    iters = result.iterations
    res = true_residual(p, gfu)
    return {
        "method": "ABMC-COCR (IC, block=4, colors=4)",
        "iterations": iters,
        "t_setup": round(t_setup, 3),
        "t_solve": round(t_solve, 3),
        "t_total": round(t_setup + t_solve, 3),
        "ms_per_iter": round(t_solve / max(iters, 1) * 1000, 1),
        "true_residual": float(res),
        "converged": result.converged,
    }


def run_single(mesh_name):
    if not mesh_name.endswith(".msh"):
        mesh_name += ".msh"

    print(f"\n{'='*90}")
    print(f"Mesh: {mesh_name}, f={freq/1e3:.0f} kHz, tol={tol}")
    print(f"{'='*90}")

    p = setup_problem(mesh_name)
    ndof = p["fes"].ndof
    print(f"  ne={p['ne']:,}, HCurl DOFs={ndof:,}, H1 DOFs={p['h1_fes'].ndof:,}")

    results = []

    # 1. Compact AMS + COCR
    print(f"\n  [1] Compact AMS + COCR ...", flush=True)
    r1 = run_compact_ams_cocr(p)
    print(f"      {r1['iterations']} iter, {r1['t_total']:.1f}s total, "
          f"{r1['ms_per_iter']:.1f} ms/iter, res={r1['true_residual']:.2e}")
    results.append(r1)

    # 2. ABMC-ICCG
    print(f"\n  [2] ABMC-ICCG ...", flush=True)
    r2 = run_abmc_iccg(p)
    status = "OK" if r2['converged'] else "FAILED"
    print(f"      {r2['iterations']} iter, {r2['t_total']:.1f}s total, "
          f"{r2['ms_per_iter']:.1f} ms/iter, res={r2['true_residual']:.2e} [{status}]")
    results.append(r2)

    # 3. ABMC-COCR
    print(f"\n  [3] ABMC-COCR ...", flush=True)
    r3 = run_abmc_cocr(p)
    status = "OK" if r3['converged'] else "FAILED"
    print(f"      {r3['iterations']} iter, {r3['t_total']:.1f}s total, "
          f"{r3['ms_per_iter']:.1f} ms/iter, res={r3['true_residual']:.2e} [{status}]")
    results.append(r3)

    # Speedup
    t_ams = r1['t_total']
    for r in results[1:]:
        if r['converged'] and r['t_total'] > 0:
            r['speedup_vs_ams'] = round(r['t_total'] / t_ams, 1)
        else:
            r['speedup_vs_ams'] = None
    r1['speedup_vs_ams'] = 1.0

    return {
        "mesh": mesh_name,
        "ne": p["ne"],
        "ndof_hcurl": ndof,
        "ndof_h1": p["h1_fes"].ndof,
        "results": results,
    }


def main():
    if "--all" in sys.argv:
        mesh_names = ALL_MESHES
    elif len(sys.argv) > 1:
        mesh_names = [a for a in sys.argv[1:] if not a.startswith("-")]
    else:
        mesh_names = ["mesh1_3.5T"]

    all_cases = []
    for name in mesh_names:
        case = run_single(name)
        all_cases.append(case)

    # Summary table
    print(f"\n{'='*110}")
    hdr = (f"{'Mesh':<16s} {'DOFs':>9s} | {'Method':<38s} {'Iters':>6s} "
           f"{'Total[s]':>9s} {'ms/iter':>8s} {'xAMS':>6s} {'Residual':>10s}")
    print(hdr)
    print(f"{'-'*110}")
    for case in all_cases:
        first = True
        for r in case['results']:
            mesh_col = f"{case['mesh']:<16s}" if first else " " * 16
            dof_col = f"{case['ndof_hcurl']:>9,d}" if first else " " * 9
            sp = f"{r['speedup_vs_ams']:.1f}x" if r.get('speedup_vs_ams') else "N/A"
            cv = "" if r['converged'] else " [FAIL]"
            print(f"{mesh_col} {dof_col} | {r['method']:<38s} {r['iterations']:>6d} "
                  f"{r['t_total']:>9.1f} {r['ms_per_iter']:>8.1f} {sp:>6s} "
                  f"{r['true_residual']:>10.2e}{cv}")
            first = False
        print(f"{'-'*110}")
    print(f"{'='*110}")
    print(f"xAMS = how many times SLOWER than Compact AMS + COCR")

    # EMD reference
    print(f"\nEMD paper reference (Hiruma SA-26-001, 3,670,328 DOFs, 30 kHz):")
    print(f"  IC only:              15,838 iter, 5964.8s  (376.7 ms/iter)")
    print(f"  EMD (IC+AMG V-cycle):  4,069 iter, 1716.8s  (422.0 ms/iter)")
    print(f"  EMD (IC+GenEO 24dom):  1,004 iter,  550.8s  (548.6 ms/iter)")

    # Save JSON
    out = {
        "timestamp": datetime.now().isoformat(),
        "hostname": platform.node(),
        "benchmark": "ams_vs_abmc",
        "physics": {
            "frequency_hz": freq, "sigma_cu": sigma_cu,
            "mu_r_core": mu_r_core, "tol": tol, "maxiter": maxiter,
        },
        "environment": {
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
        },
        "emd_reference": {
            "paper": "Hiruma SA-26-001, 2026.3.5",
            "ndof": 3670328,
            "results": [
                {"method": "IC only", "iterations": 15838, "t_solve": 5964.8},
                {"method": "EMD (IC+AMG V)", "iterations": 4069, "t_solve": 1716.8},
                {"method": "EMD (IC+AMG W)", "iterations": 2935, "t_solve": 1552.9},
                {"method": "EMD (IC+GenEO 24)", "iterations": 1004, "t_solve": 550.8},
            ],
        },
        "cases": all_cases,
    }
    json_path = os.path.join(os.path.dirname(__file__),
                             "results_ams_vs_abmc.json")
    with open(json_path, "w") as fp:
        json.dump(out, fp, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {json_path}")


if __name__ == "__main__":
    main()
