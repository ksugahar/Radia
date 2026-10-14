"""Bounded ICCG/AMS study. Run on an idle compute host, not a usage host.

All prototype extensions are serial and private to the supplied build directory.
The installed Radia is used only for the existing AMS/IC comparison.
"""
import argparse
from datetime import datetime, timezone
import importlib
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time

# One BLAS thread; NGSolve threads are selected explicitly in the AMS study.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
import scipy
from scipy import sparse
import psutil


def canonical(a):
    a = sparse.csr_matrix(a)
    a.sum_duplicates()
    a.sort_indices()
    return a


def solve(module, a, b, conjugate=False, abmc=False, scaling=True,
          profile=False, maxiter=3000, auto_shift=True):
    a = canonical(a)
    function = module.complex if np.iscomplexobj(a.data) else module.real
    return function(a.indptr, a.indices, a.data, b, conjugate, abmc,
                    scaling, profile, 1e-10, maxiter, auto_shift, 1.0)


def laplacian(side, dim):
    line = sparse.diags([-np.ones(side-1), 2*np.ones(side), -np.ones(side-1)], [-1, 0, 1])
    a = line
    for _ in range(dim - 1):
        a = sparse.kronsum(a, line, format="csr")
    return canonical(a)


def gauge(a):
    q = np.exp(1j * np.random.default_rng(12).uniform(-np.pi, np.pi, a.shape[0]))
    return canonical(sparse.diags(q) @ a @ sparse.diags(q.conj()))


def residual(a, b, x):
    return float(np.linalg.norm(a @ x - b) / np.linalg.norm(b))


def correctness(modules):
    rng = np.random.default_rng(20261004)
    c = rng.normal(size=(12, 12)) + 1j*rng.normal(size=(12, 12))
    dense = c @ c.conj().T + 2*np.eye(12)
    a = canonical(dense)
    b = rng.normal(size=12) + 1j*rng.normal(size=12)
    reference = np.linalg.solve(dense, b)
    checks = []
    for name in ("hermitian", "combined"):
        for abmc in (False, True):
            for scaling in (False, True):
                r = solve(modules[name], a, b, True, abmc, scaling, auto_shift=False)
                err = float(np.linalg.norm(r["x"] - reference) / np.linalg.norm(reference))
                checks.append(dict(case="dense_nonreal_HPD", variant=name, abmc=abmc,
                    scaling=scaling, iterations=r["iterations"], relative_error=err,
                    true_residual=residual(a, b, r["x"]),
                    passed=bool(r["converged"] and r["iterations"] == 1 and err < 1e-12)))
    # Counterexample: changing only CG's inner product does not make IC Hermitian.
    wrong = solve(modules["baseline"], a, b, True, maxiter=1, auto_shift=False)
    checks.append(dict(case="transpose_IC_is_not_Hermitian_IC", expected_negative_control=True,
                       true_residual=residual(a, b, wrong["x"]),
                       passed=bool(residual(a, b, wrong["x"]) > 1e-3)))
    for abmc in (False, True):
        for scaling in (False, True):
            a = laplacian(12, 2).astype(complex) + 0.4j*sparse.eye(144)
            b = np.ones(144, dtype=complex)
            old = solve(modules["baseline"], a, b, False, abmc, scaling)
            new = solve(modules["hermitian"], a, b, False, abmc, scaling)
            checks.append(dict(case="complex_symmetric_unchanged", abmc=abmc, scaling=scaling,
                passed=bool(old["converged"] and new["converged"] and
                            np.array_equal(old["x"], new["x"]) and old["history"] == new["history"])))
    # A real matrix cast to complex would miss the conjugation defect.
    a = gauge(laplacian(12, 2))
    exact = rng.normal(size=144) + 1j*rng.normal(size=144)
    for magnitude in (1e-100, 1.0, 1e100):
        b = a @ exact * magnitude
        r = solve(modules["hermitian"], a, b, True)
        rr = residual(a, b, r["x"])
        checks.append(dict(case="nonreal_HPD_RHS_scale", rhs_scale=magnitude,
                           true_residual=rr, passed=bool(r["converged"] and rr < 1e-8)))
    if not all(c["passed"] for c in checks):
        raise AssertionError(json.dumps(checks, indent=2))
    return checks


def benchmarks(modules, repeats):
    rows = []
    for dim, side in ((2, 128), (2, 256), (2, 384), (3, 24), (3, 40), (3, 50)):
        base = laplacian(side, dim)
        for herm in (False, True):
            a = gauge(base) if herm else base
            rng = np.random.default_rng(614)
            exact = rng.normal(size=a.shape[0])
            if herm:
                exact = exact + 1j*rng.normal(size=len(exact))
            b = a @ exact
            names = ("hermitian", "combined") if herm else ("baseline", "natural")
            samples = {name: [] for name in names}
            previous = None
            # Warm both DLLs, then alternate their timing order each repetition.
            for name in names:
                solve(modules[name], a, b, herm)
            for repeat in range(repeats):
                for name in names[::1 if repeat % 2 == 0 else -1]:
                    r = solve(modules[name], a, b, herm)
                    rr = residual(a, b, r["x"])
                    if not r["converged"] or rr > 1e-8:
                        raise AssertionError((name, dim, side, rr))
                    identity = (r["x"], r["history"], r["actual_shift"])
                    if previous is not None:
                        if not (np.array_equal(previous[0], identity[0]) and previous[1:] == identity[1:]):
                            raise AssertionError("natural order changed numerical results")
                    previous = identity
                    samples[name].append({k: r[k] for k in ("setup_s", "krylov_s", "total_s", "spmv_s")})
            for name in names:
                prof = solve(modules[name], a, b, herm, profile=True)
                row = dict(dimension=dim, side=side, n=a.shape[0], nnz=a.nnz,
                    hermitian=herm, variant=name, threads=1, repeats=repeats,
                    iterations=prof["iterations"], true_residual=residual(a, b, prof["x"]),
                    actual_shift=prof["actual_shift"], bit_identical_pair=True,
                    samples=samples[name], median={k: statistics.median(s[k] for s in samples[name])
                                                  for k in samples[name][0]},
                    profile={k: prof[k] for k in ("setup_s", "krylov_s", "apply_and_dot_s", "apply_calls")},
                    process_rss_bytes=psutil.Process().memory_info().rss)
                rows.append(row)
                print(json.dumps({k: row[k] for k in ("variant", "n", "iterations", "true_residual", "median")}), flush=True)
    return rows


def ams_study(repeats):
    import ngsolve as ng
    import radia.sparsesolv_ngsolve as ss
    from netgen.csg import unit_cube
    rows = []
    for h in (0.35, 0.16, 0.08):
        ng.SetNumThreads(1)
        mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=h))
        fes = ng.HCurl(mesh, order=1, nograds=True, dirichlet=".*")
        u, v = fes.TnT()
        grad, _ = fes.CreateGradient()
        coords = [[mesh.ngmesh.Points()[i+1][j] for i in range(mesh.nv)] for j in range(3)]
        free = np.asarray(list(fes.FreeDofs()), dtype=bool)
        for contrast, beta in ((1.0, 1.0), (1000.0, 1.0), (1000.0, 1e-6), (1000.0, 0.0)):
            a = ng.BilinearForm(fes)
            coefficient = ng.IfPos(ng.x-0.5, contrast, 1.0)
            a += coefficient * ng.curl(u)*ng.curl(v)*ng.dx + beta*u*v*ng.dx
            a.Assemble()
            exact = ng.GridFunction(fes)
            exact.vec.FV().NumPy()[:] = np.random.default_rng(1004).normal(size=fes.ndof)*free
            rhs = exact.vec.CreateVector()
            rhs.data = a.mat * exact.vec
            rhs.FV().NumPy()[~free] = 0
            norm_rhs = ng.Norm(rhs)
            for threads in (1, 4):
                ng.SetNumThreads(threads)
                for method in ("IC_fixed_shift_1.05", "AMS"):
                    samples = []
                    acceptance = []
                    result = None
                    for repeat in range(repeats + 1):
                        t = time.perf_counter()
                        if method == "AMS":
                            pre = ss.HypreBasedAMSPreconditioner(a.mat, grad, freedofs=fes.FreeDofs(),
                                coord_x=coords[0], coord_y=coords[1], coord_z=coords[2],
                                beta_zero=(beta == 0), cycle_type=1, subspace_solver=0, print_level=0)
                        else:
                            pre = ss.ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
                        setup_s = time.perf_counter()-t
                        sol = ng.GridFunction(fes)
                        inv = ss.NativePCG(a.mat, pre, fes.FreeDofs())
                        with ng.TaskManager():
                            t = time.perf_counter()
                            iterations, reported, converged = inv.Solve(rhs, sol.vec, 1e-10, 3000)
                            solve_s = time.perf_counter()-t
                        r = rhs.CreateVector()
                        r.data = rhs-a.mat*sol.vec
                        rr = float(np.linalg.norm(r.FV().NumPy()[free])/norm_rhs)
                        field_error = float(np.sqrt(abs(ng.Integrate(
                            (ng.curl(sol)-ng.curl(exact))*(ng.curl(sol)-ng.curl(exact)), mesh)) /
                            ng.Integrate(ng.curl(exact)*ng.curl(exact), mesh)))
                        # Repeat the same application to expose state-dependent output.
                        z1, z2 = rhs.CreateVector(), rhs.CreateVector()
                        with ng.TaskManager():
                            z1.data = pre * rhs
                            z2.data = pre * rhs
                        repeat_error = float(ng.Norm(z1-z2)/max(ng.Norm(z1), 1e-300))
                        delta = ng.GridFunction(fes)
                        delta.vec.data = z1-z2
                        reference = ng.GridFunction(fes)
                        reference.vec.data = z1
                        curl_repeat_error = float(np.sqrt(abs(ng.Integrate(ng.curl(delta)*ng.curl(delta), mesh)) /
                            max(abs(ng.Integrate(ng.curl(reference)*ng.curl(reference), mesh)), 1e-300)))
                        ad, az = rhs.CreateVector(), rhs.CreateVector()
                        ad.data = a.mat * delta.vec
                        az.data = a.mat * z1
                        action_repeat_error = float(np.linalg.norm(ad.FV().NumPy()[free]) /
                                                    max(np.linalg.norm(az.FV().NumPy()[free]), 1e-300))
                        result = dict(iterations=iterations, true_residual=rr, reported_residual=reported,
                            curl_relative_error=field_error, repeated_apply_error=repeat_error,
                            repeated_apply_curl_error=curl_repeat_error,
                            repeated_apply_matrix_action_error=action_repeat_error,
                            passed=bool(converged and rr < 1e-7 and field_error < 1e-5 and repeat_error < 1e-12))
                        acceptance.append(result)
                        if repeat:
                            samples.append(dict(setup_s=setup_s, krylov_s=solve_s, total_s=setup_s+solve_s))
                        del inv, pre
                    row = dict(method=method, maxh=h, n=fes.ndof, free_dofs=int(free.sum()),
                        contrast=contrast, beta=beta, threads=threads, samples=samples,
                        median={k: statistics.median(s[k] for s in samples) for k in samples[0]}, **result)
                    row["acceptance_samples"] = acceptance
                    row["passed"] = all(r["passed"] for r in acceptance)
                    rows.append(row)
                    print(json.dumps(row), flush=True)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--source-sha", required=True)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--part", choices=("core", "ams"), default="core")
    args = p.parse_args()
    if args.repeats < 2:
        p.error("at least two timed repetitions required")
    sys.path.insert(0, str(args.build / "build/Release"))
    start = time.perf_counter()
    import ngsolve
    import radia
    result = dict(schema="cae-ai-lab.solver-run.v1", case="ICCG Hermitian and throughput study: " + args.part,
        solver="radia-ngsolve", created_at_utc=datetime.now(timezone.utc).isoformat(),
        source_sha=args.source_sha, platform_class=platform.system(),
        tool_versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                           ngsolve=ngsolve.__version__, radia=radia.__version__),
        tolerances=dict(core_true_relative=1e-8, ams_true_relative=1e-7, curl_relative=1e-5),
        variant_sha256=json.loads((args.build / "variants.json").read_text()),
        run=dict(command=" ".join(sys.argv), workdir=str(Path.cwd()), exit_code=1))
    try:
        if args.part == "core":
            modules = {v: importlib.import_module("iccg_" + v) for v in ("baseline", "natural", "hermitian", "combined")}
            result["correctness"] = correctness(modules)
            result["benchmarks"] = benchmarks(modules, args.repeats)
            passed = True
        else:
            result["benchmarks"] = ams_study(args.repeats)
            passed = all(row["passed"] for row in result["benchmarks"])
        result["pass"] = passed
        if not passed:
            result["failure"] = dict(stage=args.part, message="One or more numerical acceptance checks failed",
                next_action="Inspect all acceptance samples; do not promote an AMS default from a failed study")
        result["run"]["exit_code"] = 0 if passed else 1
        result["checks"] = dict(ran_to_completion=True, validation_passed=passed, result_files_exist=True)
        result["errors"] = dict(max_true_relative=max(row["true_residual"] for row in result["benchmarks"]))
        result["timing_breakdown_s"] = dict(total=time.perf_counter()-start)
        result["verification"] = dict(method="Independent original-system residual; dense direct oracle; paired bit identity; curl error for singular systems")
    except Exception as error:
        result["pass"] = False
        result["checks"] = dict(ran_to_completion=False, validation_passed=False)
        result["failure"] = dict(stage=args.part, message=repr(error), next_action="Inspect failure and repair the experiment before drawing conclusions")
        raise
    finally:
        result["run"]["duration_s"] = time.perf_counter()-start
        result["result_files"] = [str(args.output)]
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
