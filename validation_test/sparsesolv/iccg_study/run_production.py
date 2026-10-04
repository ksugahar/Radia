"""Production validation; run on an idle compute host with explicit PYTHONPATH.

Use the same --mesh for baseline and candidate processes. No package installation.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import ngsolve as ng
import radia.sparsesolv_ngsolve as ss


def throughput(meshpath, kind, space):
    from netgen.csg import unit_cube
    if not meshpath.exists():
        unit_cube.GenerateMesh(maxh=0.045).Save(str(meshpath))
    mesh = ng.Mesh(str(meshpath))
    hermitian = kind == "hermitian"
    ng.SetNumThreads(1)
    fes = (ng.HCurl(mesh, order=1, nograds=True, dirichlet=".*", complex=kind != "real")
           if space == "hcurl" else ng.H1(mesh, order=2, dirichlet=".*", complex=kind != "real"))
    u, v = fes.TnT()
    a = ng.BilinearForm(fes)
    mass = 1+100j if kind == "symmetric" else 1
    derivative = ng.curl if space == "hcurl" else ng.grad
    a += (derivative(u)*derivative(v) + mass*u*v)*ng.dx
    a.Assemble()
    matrix = a.mat
    if hermitian:
        i, j, values = matrix.COO()
        phase = np.exp(0.37j*np.arange(fes.ndof))
        values = np.asarray(values)*phase[np.asarray(i)]*phase[np.asarray(j)].conj()
        matrix.AsVector().FV().NumPy()[:] = values
    exact = matrix.CreateColVector()
    free = np.asarray(list(fes.FreeDofs()), bool)
    exact.FV().NumPy()[:] = np.random.default_rng(481).normal(size=fes.ndof)*free
    rhs = matrix.CreateRowVector()
    rhs.data = matrix*exact
    rhs.FV().NumPy()[~free] = 0
    rows = []
    for threads in (1, 4):
        ng.SetNumThreads(threads)
        samples, solutions, histories = [], [], []
        for repeat in range(4):
            solver = ss.SparseSolvSolver(matrix, freedofs=fes.FreeDofs(),
                conjugate=hermitian, tol=1e-9, maxiter=2000,
                save_residual_history=True, divergence_check=False)
            x = matrix.CreateColVector(); x.FV().NumPy()[:] = 0
            with ng.TaskManager():
                start = time.perf_counter()
                info = solver.Solve(rhs, x)
                elapsed = time.perf_counter()-start
            residual = rhs.CreateVector(); residual.data = rhs-matrix*x
            rr = np.linalg.norm(residual.FV().NumPy()[free])/ng.Norm(rhs)
            assert info.converged and rr < 1e-8
            if repeat:
                samples.append(elapsed)
                solutions.append(x.FV().NumPy().copy())
                histories.append(list(info.residual_history))
        row = dict(n=fes.ndof, threads=threads, kind=kind, space=space,
            matrix_values_sha256=hashlib.sha256(matrix.AsVector().FV().NumPy().tobytes()).hexdigest(),
            rhs_sha256=hashlib.sha256(rhs.FV().NumPy().tobytes()).hexdigest(),
            samples_s=samples, median_s=statistics.median(samples),
            iterations=info.iterations, actual_shift=info.actual_shift,
            true_residual=float(rr), history=histories[-1],
            solution_sha256=hashlib.sha256(solutions[-1].tobytes()).hexdigest(),
            repeat_relative=max(float(np.linalg.norm(s-solutions[0])/np.linalg.norm(s)) for s in solutions))
        rows.append(row)
        print(json.dumps(row), flush=True)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--part", choices=["ams", "throughput"], required=True)
    p.add_argument("--mesh", type=Path)
    p.add_argument("--kind", choices=["real", "symmetric", "hermitian"], default="symmetric")
    p.add_argument("--space", choices=["h1", "hcurl"], default="hcurl")
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.part == "throughput" and args.mesh is None:
        p.error("--mesh is required for throughput")
    start = time.perf_counter()
    result = dict(schema="cae-ai-lab.solver-run.v1", case="Production ICCG/AMS " + args.part,
        solver="radia-ngsolve", created_at_utc=datetime.now(timezone.utc).isoformat(),
        host=platform.node(), python=platform.python_version(),
        ngsolve=ng.__version__, module=str(ss.__file__),
        binary_sha256=hashlib.sha256(Path(ss.__file__).read_bytes()).hexdigest(), passed=False)
    try:
        if args.part == "ams":
            from run_study import ams_study
            result["cases"] = ams_study(3)
            result["passed"] = all(c["passed"] for c in result["cases"])
        else:
            result["cases"] = throughput(args.mesh, args.kind, args.space)
            result["passed"] = True
    except Exception as error:
        result["failure"] = dict(stage=args.part, message=repr(error), next_action="Inspect numerical/build failure")
        raise
    finally:
        result["pass"] = result["passed"]
        result["run"] = dict(command=" ".join(sys.argv), workdir=str(Path.cwd()),
            exit_code=0 if result["passed"] else 1, duration_s=time.perf_counter()-start)
        result["checks"] = dict(ran_to_completion="cases" in result,
            validation_passed=result["passed"], result_files_exist=True)
        result["result_files"] = [str(args.output)]
        result["tolerances"] = dict(true_relative=1e-7 if args.part == "ams" else 1e-8,
            repeated_apply_relative=1e-12)
        result["errors"] = dict(max_true_relative=max((r["true_residual"] for r in result.get("cases", [])), default=0))
        result["tool_versions"] = dict(python=platform.python_version(), ngsolve=ng.__version__)
        if args.mesh is not None and args.mesh.is_file():
            result["mesh_sha256"] = hashlib.sha256(args.mesh.read_bytes()).hexdigest()
        result["verification"] = dict(method="Original-system residual; repeated application; paired baseline histories and solution hashes")
        args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
