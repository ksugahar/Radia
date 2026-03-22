"""Benchmark: COCR vs GMRES(40) with Compact AMS preconditioner.

COCR exploits complex symmetry (A^T = A) of eddy current system.
Short recurrences: O(n) memory, no restart needed.
1 preconditioner + 2 SpMV per iter (vs GMRES: 1 preconditioner + 1 SpMV + O(j) GS).

Usage:
    python bench_cocr_vs_gmres.py [mesh_name]
"""

import json
import os
import platform
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

maxiter = 2000
tol = 1e-10


def setup_problem(mesh_name):
    mesh_file = os.path.join(os.path.dirname(__file__), mesh_name)
    print(f"Loading mesh: {mesh_file}", flush=True)
    mesh = Mesh(ReadGmsh(mesh_file))
    print(f"  ne={mesh.ne:,}, nv={mesh.nv:,}", flush=True)

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

    nv = mesh.nv
    coord_x = [0.0] * nv
    coord_y = [0.0] * nv
    coord_z = [0.0] * nv
    for i in range(nv):
        p = mesh.ngmesh.Points()[i + 1]
        coord_x[i] = p[0]
        coord_y[i] = p[1]
        coord_z[i] = p[2]

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
    }


def true_residual(p, gfu):
    r = p["f"].vec.CreateVector()
    r.data = p["f"].vec - p["a"].mat * gfu.vec
    fd = p["fes"].FreeDofs()
    for i in range(p["fes"].ndof):
        if not fd[i]:
            r[i] = 0
    return float(sqrt(abs(InnerProduct(r, r)))) / p["rhs_norm"]


def main():
    mesh_name = sys.argv[1] if len(sys.argv) > 1 else "mesh1_3.5T"
    if not mesh_name.endswith(".msh"):
        mesh_name += ".msh"

    print("=" * 80)
    print("Benchmark: COCR vs GMRES(40) with Compact AMS")
    print(f"Mesh: {mesh_name}, f={freq/1e3:.0f} kHz, tol={tol}")
    print("=" * 80)

    p = setup_problem(mesh_name)
    ndof = p["fes"].ndof
    results = []

    # Build preconditioner
    t0 = time.perf_counter()
    pre = ssn.ComplexCompactAMSPreconditioner(
        a_real_mat=p["a_real"].mat, grad_mat=p["G_mat"],
        freedofs=p["fes_real"].FreeDofs(),
        coord_x=p["coord_x"], coord_y=p["coord_y"], coord_z=p["coord_z"],
        ndof_complex=ndof, cycle_type=1)
    t_setup = time.perf_counter() - t0
    print(f"\nCompact AMS setup: {t_setup:.3f}s", flush=True)

    # Case 1: COCR
    print(f"\nCase 1: COCR (complex symmetric, O(n) memory)", flush=True)
    gfu1 = GridFunction(p["fes"])
    with TaskManager():
        t0 = time.perf_counter()
        inv1 = ssn.COCRSolver(
            p["a"].mat, pre, freedofs=p["fes"].FreeDofs(),
            maxiter=maxiter, tol=tol)
        gfu1.vec.data = inv1 * p["f"].vec
        t_solve1 = time.perf_counter() - t0

    iters1 = inv1.iterations
    res1 = true_residual(p, gfu1)
    print(f"  Solve: {t_solve1:.3f}s, Iters: {iters1}, "
          f"true ||b-Ax||/||b||: {res1:.2e}", flush=True)
    results.append({
        "method": "COCR + Compact AMS",
        "iterations": iters1,
        "t_setup": round(t_setup, 3),
        "t_solve": round(t_solve1, 3),
        "t_total": round(t_setup + t_solve1, 3),
        "true_residual": res1,
    })

    # Case 2: GMRES(40)
    print(f"\nCase 2: GMRES(40)", flush=True)
    gfu2 = GridFunction(p["fes"])
    with TaskManager():
        t0 = time.perf_counter()
        inv2 = ssn.GMRESSolver(
            p["a"].mat, pre, freedofs=p["fes"].FreeDofs(),
            maxiter=maxiter, tol=tol, restart=40)
        gfu2.vec.data = inv2 * p["f"].vec
        t_solve2 = time.perf_counter() - t0

    iters2 = inv2.iterations
    res2 = true_residual(p, gfu2)
    print(f"  Solve: {t_solve2:.3f}s, Iters: {iters2}, "
          f"true ||b-Ax||/||b||: {res2:.2e}", flush=True)
    results.append({
        "method": "GMRES(40) + Compact AMS",
        "iterations": iters2,
        "t_setup": round(t_setup, 3),
        "t_solve": round(t_solve2, 3),
        "t_total": round(t_setup + t_solve2, 3),
        "true_residual": res2,
    })

    # Solution comparison
    print("\nSolution comparison:", flush=True)
    norm1 = sqrt(abs(InnerProduct(gfu1.vec, gfu1.vec)))
    norm2 = sqrt(abs(InnerProduct(gfu2.vec, gfu2.vec)))
    diff = gfu1.vec.CreateVector()
    diff.data = gfu1.vec - gfu2.vec
    err = sqrt(abs(InnerProduct(diff, diff))) / norm1
    print(f"  ||COCR||   = {norm1:.6e}")
    print(f"  ||GMRES||  = {norm2:.6e}")
    print(f"  ||COCR - GMRES|| / ||COCR|| = {err:.2e}")

    # Summary
    print()
    print("=" * 80)
    hdr = f"{'Method':<35s} {'Iters':>6s} {'Setup[s]':>9s} " \
          f"{'Solve[s]':>9s} {'Total[s]':>9s} {'TrueRes':>10s} {'ms/iter':>8s}"
    print(hdr)
    print("-" * 80)
    for r in results:
        ms = r['t_solve'] / r['iterations'] * 1000 if r['iterations'] > 0 else 0
        row = f"{r['method']:<35s} {r['iterations']:>6d} " \
              f"{r['t_setup']:>9.3f} {r['t_solve']:>9.3f} " \
              f"{r['t_total']:>9.3f} {r['true_residual']:>10.2e} {ms:>8.1f}"
        print(row)
    print("=" * 80)

    if t_solve1 < t_solve2:
        print(f"\nCOCR {t_solve2/t_solve1:.2f}x faster solve "
              f"({t_solve1:.2f}s vs {t_solve2:.2f}s)")
    else:
        print(f"\nGMRES {t_solve1/t_solve2:.2f}x faster solve "
              f"({t_solve2:.2f}s vs {t_solve1:.2f}s)")


if __name__ == "__main__":
    main()
