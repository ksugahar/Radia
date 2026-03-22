"""Benchmark: Compact AMS + COCR for complex eddy current problems.

Compact AMS + COCR exploits complex symmetry (A^T = A) of eddy current system.
COCR: O(n) memory, no restart parameter, short recurrences.

Usage:
    python bench_compact_ams.py                        # default (mesh1_3.5T)
    python bench_compact_ams.py mesh1_2.5T             # single mesh
    python bench_compact_ams.py --all                  # all available meshes
    python bench_compact_ams.py mesh1_2.5T mesh1_3.5T  # multiple meshes
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

# Physical parameters (match Hiruma SA-26-001)
mu0 = 4e-7 * np.pi
freq = 30e3
omega = 2 * np.pi * freq
eps = 1e-6
sigma_cu = 5.96e7
mu_r_core = 1000

maxiter = 2000
tol = 1e-10

ALL_MESHES = [
    "mesh1_2.5T",
    "mesh1_3.5T",
    "mesh1_4.5T",
    "mesh1_5.5T",
]


def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return mem.peak_wset / (1024 * 1024) if hasattr(mem, 'peak_wset') else mem.rss / (1024 * 1024)


def setup_problem(mesh_name):
    mesh_file = os.path.join(os.path.dirname(__file__), mesh_name)
    print(f"Loading mesh: {mesh_file}", flush=True)
    mesh = Mesh(ReadGmsh(mesh_file))
    ne = mesh.ne
    nv = mesh.nv
    print(f"  ne={ne:,}, nv={nv:,}", flush=True)

    nu_cf = 1.0 / (mu0 * IfPos(mesh.MaterialCF({"core": 1}), mu_r_core, 1.0))
    sigma_cf = mesh.MaterialCF({"cond": sigma_cu}, default=0)

    fes = HCurl(mesh, order=1, nograds=True,
                dirichlet="dirichlet", complex=True)
    u, v = fes.TnT()
    print(f"  HCurl DOFs: {fes.ndof:,}", flush=True)

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
    print(f"  H1 DOFs: {h1_fes.ndof:,}", flush=True)

    coord_x = [0.0] * nv
    coord_y = [0.0] * nv
    coord_z = [0.0] * nv
    for i in range(nv):
        p = mesh.ngmesh.Points()[i + 1]
        coord_x[i] = p[0]
        coord_y[i] = p[1]
        coord_z[i] = p[2]

    # Compute rhs norm on free DOFs
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
        "rhs_norm": rhs_norm,
        "ne": ne, "nv": nv,
    }


def true_residual(p, gfu):
    """Compute ||b - Ax||/||b|| on free DOFs only."""
    r = p["f"].vec.CreateVector()
    r.data = p["f"].vec - p["a"].mat * gfu.vec
    fd = p["fes"].FreeDofs()
    for i in range(p["fes"].ndof):
        if not fd[i]:
            r[i] = 0
    return float(sqrt(abs(InnerProduct(r, r)))) / p["rhs_norm"]


def run_single(mesh_name):
    """Run benchmark for a single mesh, return result dict."""
    if not mesh_name.endswith(".msh"):
        mesh_name += ".msh"

    print("\n" + "=" * 80)
    print(f"Benchmark: Compact AMS + COCR  |  {mesh_name}")
    print(f"f={freq/1e3:.0f} kHz, sigma={sigma_cu:.2e} S/m, mu_r={mu_r_core}, tol={tol}")
    print("=" * 80)

    p = setup_problem(mesh_name)
    ndof_hcurl = p["fes"].ndof
    ndof_h1 = p["h1_fes"].ndof

    # Setup preconditioner
    print("\nPreconditioner setup:", flush=True)
    t0 = time.perf_counter()
    pre = ssn.ComplexCompactAMSPreconditioner(
        a_real_mat=p["a_real"].mat, grad_mat=p["G_mat"],
        freedofs=p["fes_real"].FreeDofs(),
        coord_x=p["coord_x"], coord_y=p["coord_y"], coord_z=p["coord_z"],
        ndof_complex=ndof_hcurl, cycle_type=1, print_level=1)
    t_setup = time.perf_counter() - t0
    print(f"\n  Compact AMS setup: {t_setup:.3f}s", flush=True)

    # Solve
    print(f"\nCompact AMS + COCR:", flush=True)
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
    ms_per_iter = (t_solve / iters * 1000) if iters > 0 else 0
    t_total = t_setup + t_solve
    peak_mem = get_peak_memory_mb()

    print(f"  Iters:     {iters}")
    print(f"  Setup:     {t_setup:.3f}s")
    print(f"  Solve:     {t_solve:.3f}s ({ms_per_iter:.1f} ms/iter)")
    print(f"  Total:     {t_total:.3f}s")
    print(f"  Residual:  {res:.2e}")
    print(f"  Memory:    {peak_mem:.0f} MB")

    return {
        "mesh": mesh_name,
        "ne": p["ne"],
        "nv": p["nv"],
        "ndof_hcurl": ndof_hcurl,
        "ndof_h1": ndof_h1,
        "method": "Compact AMS + COCR",
        "iterations": iters,
        "t_setup": round(t_setup, 3),
        "t_solve": round(t_solve, 3),
        "t_total": round(t_total, 3),
        "ms_per_iter": round(ms_per_iter, 1),
        "true_residual": float(res),
        "converged": iters < maxiter,
        "peak_memory_mb": round(peak_mem, 1),
    }


def main():
    if "--all" in sys.argv:
        mesh_names = ALL_MESHES
    elif len(sys.argv) > 1:
        mesh_names = [a for a in sys.argv[1:] if not a.startswith("-")]
    else:
        mesh_names = ["mesh1_3.5T"]

    results = []
    for name in mesh_names:
        r = run_single(name)
        results.append(r)

    # Summary table
    print("\n" + "=" * 100)
    hdr = (f"{'Mesh':<18s} {'DOFs':>10s} {'Iters':>6s} {'Setup[s]':>9s} "
           f"{'Solve[s]':>9s} {'Total[s]':>9s} {'ms/iter':>8s} {'Residual':>10s} "
           f"{'Mem[MB]':>8s}")
    print(hdr)
    print("-" * 100)
    for r in results:
        row = (f"{r['mesh']:<18s} {r['ndof_hcurl']:>10,d} {r['iterations']:>6d} "
               f"{r['t_setup']:>9.3f} {r['t_solve']:>9.3f} {r['t_total']:>9.3f} "
               f"{r['ms_per_iter']:>8.1f} {r['true_residual']:>10.2e} "
               f"{r['peak_memory_mb']:>8.0f}")
        print(row)
    print("=" * 100)

    # Save JSON
    out = {
        "timestamp": datetime.now().isoformat(),
        "hostname": platform.node(),
        "benchmark": "compact_ams_cocr",
        "solver": {
            "krylov": "COCR",
            "preconditioner": "ComplexCompactAMS (fused Re/Im + DualMult AMG)",
            "cycle_type": 1,
            "amg": "CompactAMG (PMIS coarsening, l1-Jacobi, V-cycle)",
        },
        "physics": {
            "frequency_hz": freq,
            "sigma_cu": sigma_cu,
            "mu_r_core": mu_r_core,
            "epsilon_regularization": eps,
            "tol": tol,
            "maxiter": maxiter,
        },
        "environment": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "results": results,
    }

    json_path = os.path.join(os.path.dirname(__file__),
                             "results_compact_ams.json")
    with open(json_path, "w") as fp:
        json.dump(out, fp, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {json_path}")


if __name__ == "__main__":
    main()
