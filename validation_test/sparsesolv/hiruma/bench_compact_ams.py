"""Benchmark: Compact AMS + COCR for complex eddy current problems.

Compact AMS + COCR exploits complex symmetry (A^T = A) of eddy current system.
COCR: O(n) memory, no restart parameter, short recurrences.

Usage:
    python bench_compact_ams.py                        # tracked mesh1_2.5T fixture
    python bench_compact_ams.py mesh1_2.5T             # single mesh
    python bench_compact_ams.py --all                  # all available meshes
    python bench_compact_ams.py mesh1_2.5T mesh1_3.5T  # multiple meshes
    python bench_compact_ams.py --verify-baseline      # validation gate
"""

import argparse
import hashlib
import json
import os
import platform
import psutil
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ngsolve import *
import radia.sparsesolv_ngsolve as ssn
from hiruma_mesh import load_hiruma_mesh, mesh_path

HERE = Path(__file__).resolve().parent
BASELINE_PATH = HERE / "compact_ams_baseline.json"
DEFAULT_OUTPUT = (Path(r"C:\temp") / "radia-validation" /
                  "sparsesolv_hiruma" / "compact_ams_results.json")

# Physical parameters (match Hiruma SA-26-001)
mu0 = 4e-7 * np.pi
freq = 30e3
omega = 2 * np.pi * freq
preconditioner_eps = 1e-6
sigma_cu = 5.96e7
mu_r_core = 1000

maxiter = 2000
tol = 1e-10

ALL_MESHES = [
    "mesh1_2.5T",
    "mesh1_3.5T",
    "mesh1_4.5T",
    "mesh1_5.5T",
    "mesh1_20.5T",
]


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return mem.peak_wset / (1024 * 1024) if hasattr(mem, 'peak_wset') else mem.rss / (1024 * 1024)


def setup_problem(mesh_name):
    print(f"Loading mesh: {mesh_path(mesh_name)}", flush=True)
    mesh = load_hiruma_mesh(mesh_name)
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
    # The shift regularizes only the AMS surrogate.  Adding it to the physical
    # system perturbs the operator, not only its gauge representative. Compare
    # curl(A), conductor losses and the common unshifted residual before
    # claiming physical equivalence; a scalar quadratic-form match is not proof.
    a_real += preconditioner_eps * nu_cf * u_r * v_r * dx
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
    mesh_name = mesh_path(mesh_name).name

    print("\n" + "=" * 80)
    print(f"Benchmark: Compact AMS + COCR  |  {mesh_name}")
    print(f"f={freq/1e3:.0f} kHz, sigma={sigma_cu:.2e} S/m, mu_r={mu_r_core}, tol={tol}")
    print("=" * 80)

    resolved_mesh = mesh_path(mesh_name)
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
        "mesh_sha256": sha256_file(resolved_mesh),
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
        "converged": iters < maxiter and res <= 2.0 * tol,
        "peak_memory_mb": round(peak_mem, 1),
    }


def verify_baseline(result, baseline_path=BASELINE_PATH):
    """Check stable numerical ranges; timings remain reported, not gated."""
    baseline = json.loads(Path(baseline_path).read_text(encoding="utf-8"))
    fixture = baseline["fixture"]
    acceptance = baseline["acceptance"]
    failures = []
    exact = {
        "mesh_sha256": fixture["sha256"],
        "ne": fixture["ne"],
        "nv": fixture["nv"],
        "ndof_hcurl": fixture["ndof_hcurl"],
        "ndof_h1": fixture["ndof_h1"],
    }
    for key, expected in exact.items():
        if result[key] != expected:
            failures.append(f"{key}: got {result[key]!r}, expected {expected!r}")
    if not acceptance["iterations_min"] <= result["iterations"] <= acceptance["iterations_max"]:
        failures.append(
            f"iterations: got {result['iterations']}, expected "
            f"[{acceptance['iterations_min']}, {acceptance['iterations_max']}]")
    if result["true_residual"] > acceptance["true_residual_max"]:
        failures.append(
            f"true_residual: got {result['true_residual']:.3e}, expected <= "
            f"{acceptance['true_residual_max']:.3e}")
    if failures:
        raise AssertionError("Hiruma baseline failed:\n  " + "\n  ".join(failures))
    print("\nHiruma numerical baseline: PASS")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("meshes", nargs="*", default=["mesh1_2.5T"])
    parser.add_argument("--all", action="store_true",
                        help="run the full optional five-mesh scaling sweep")
    parser.add_argument("--verify-baseline", action="store_true",
                        help="gate the tracked mesh1_2.5T result")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="result JSON (default: C:\\temp\\radia-validation)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    mesh_names = ALL_MESHES if args.all else args.meshes
    if args.verify_baseline and mesh_names != ["mesh1_2.5T"]:
        raise ValueError("--verify-baseline requires only mesh1_2.5T")

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
        "schema": "radia.validation.sparsesolv-hiruma.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
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
            "system_epsilon_regularization": 0.0,
            "preconditioner_epsilon_regularization": preconditioner_eps,
            "epsilon_placement": "preconditioner_only",
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

    if args.verify_baseline:
        verify_baseline(results[0])

    json_path = args.output.resolve()
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w", encoding="utf-8") as fp:
        json.dump(out, fp, indent=2, ensure_ascii=False)
        fp.write("\n")
    print(f"\nResults saved to {json_path}")


if __name__ == "__main__":
    main()
