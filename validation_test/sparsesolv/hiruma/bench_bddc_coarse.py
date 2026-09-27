"""Benchmark: the wirebasket (coarse) solver inside NGSolve BDDC at order 2.

The Hiruma 30 kHz eddy-current problem on HCurl(order=2, nograds=True) is
solved with BDDC + COCR.  At order 2 the BDDC wirebasket is exactly the
lowest-order edge space, so its Schur complement is a lowest-order Nedelec
operator and an AMS cycle can stand in for the direct factorization:

    coarse=direct          NGSolve BDDC, sparse direct wirebasket inverse
    coarse=sparsesolv_ams  one Compact AMS cycle on the wirebasket

``sparsesolv_ams`` is registered with NGSolve's preconditioner registry when
the sparsesolv module is imported, so it is selected through BDDC's own
``coarsetype`` flag.  The system carries a small ``eps * nu`` mass: without
it the sigma=0 air leaves the lowest-order gradients in the kernel and the
direct wirebasket factorization breaks down.

Usage:
    python bench_bddc_coarse.py mesh1_2.5T --coarse direct
    python bench_bddc_coarse.py mesh1_5.5T --coarse sparsesolv_ams \
        --mesh-dir <dir with mesh1_5.5T.vol> --ssn-dir <dir with the module>
"""

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (Path(r"C:\temp") / "radia-validation" /
                  "sparsesolv_hiruma" / "bddc_coarse_results.json")

mu0 = 4e-7 * np.pi
freq = 30e3
omega = 2 * np.pi * freq
sigma_cu = 5.96e7
mu_r_core = 1000


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def peak_memory_mb():
    try:
        import psutil
        mem = psutil.Process(os.getpid()).memory_info()
        return getattr(mem, "peak_wset", mem.rss) / 2**20
    except ImportError:
        return float("nan")


def import_sparsesolv(ssn_dir):
    if ssn_dir:
        import radia  # noqa: F401  (registers the NGSolve/MKL DLL directories)
        sys.path.insert(0, str(ssn_dir))
        import sparsesolv_ngsolve as ssn
    else:
        import radia.sparsesolv_ngsolve as ssn
    return ssn


def run(args):
    import ngsolve
    from ngsolve import (BilinearForm, CF, GridFunction, HCurl, IfPos,
                         Integrate, InnerProduct, LinearForm, Mesh,
                         Preconditioner, TaskManager, curl, dx, sqrt)
    ssn = import_sparsesolv(args.ssn_dir)
    ngsolve.SetNumThreads(args.threads)

    mesh_dir = Path(args.mesh_dir) if args.mesh_dir else HERE / "meshes"
    mesh_file = mesh_dir / f"{args.mesh}.vol"
    t_start = time.perf_counter()
    with TaskManager():
        mesh = Mesh(str(mesh_file))
    nu = 1.0 / (mu0 * IfPos(mesh.MaterialCF({"core": 1}), mu_r_core, 1.0))
    sigma = mesh.MaterialCF({"cond": sigma_cu}, default=0)

    fes = HCurl(mesh, order=args.order, nograds=True,
                dirichlet="dirichlet", complex=True)
    # NGSolve's HCurl also marks the face dofs of badly shaped faces as
    # wirebasket; AMS covers only the lowest-order edge block, so this option
    # returns them to the interface (applied to every coarse type compared).
    t0 = time.perf_counter()
    face_wirebasket = 0
    if args.edge_wirebasket:
        from ngsolve import COUPLING_TYPE
        for dof in range(mesh.nedge, fes.ndof):
            if fes.CouplingType(dof) == COUPLING_TYPE.WIREBASKET_DOF:
                fes.SetCouplingType(dof, COUPLING_TYPE.INTERFACE_DOF)
                face_wirebasket += 1
    t_coupling = time.perf_counter() - t0
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu * curl(u) * curl(v) * dx
    a += 1j * omega * sigma * u * v * dx("cond")
    if args.eps > 0:
        a += args.eps * nu * u * v * dx
    flags = {"coarsetype": args.coarse}
    if args.coarse != "direct":
        flags["coarseflags"] = {"print_level": 0,
                                "cycle_type": args.cycle_type,
                                "num_smooth": args.num_smooth,
                                "eps": args.coarse_eps,
                                "cycles": args.cycles,
                                "lean_coarse": args.lean_coarse}
    else:
        flags["inverse"] = args.inverse
    pre = Preconditioner(a, "bddc", **flags)
    f = LinearForm(fes)
    f += nu * CF((0, 0, 1)) * v * dx("cond")

    with TaskManager():
        t0 = time.perf_counter()
        a.Assemble()
        f.Assemble()
        t_assemble = time.perf_counter() - t0

        gfu = GridFunction(fes)
        t0 = time.perf_counter()
        inv = ssn.COCRSolver(a.mat, pre, freedofs=fes.FreeDofs(),
                             maxiter=args.maxiter, tol=args.tol,
                             printrates=False)
        gfu.vec.data = inv * f.vec
        t_solve = time.perf_counter() - t0
        iters = inv.iterations

        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * gfu.vec
        fd = fes.FreeDofs()
        rf = r.FV().NumPy()
        bf = f.vec.FV().NumPy()
        mask = np.array(fd, dtype=bool)
        true_res = float(np.linalg.norm(rf[mask]) / np.linalg.norm(bf[mask]))

        A = gfu
        loss = Integrate(0.5 * sigma * omega**2 * InnerProduct(A, A).real,
                         mesh, definedon=mesh.Materials("cond"))
        energy = Integrate(nu * InnerProduct(curl(A), curl(A)).real, mesh)
    t_total = time.perf_counter() - t_start

    coarse_stats = {}
    if args.coarse == "sparsesolv_ams" and hasattr(ssn, "AMSCoarseStats"):
        coarse_stats = dict(ssn.AMSCoarseStats())

    result = {
        "mesh": args.mesh,
        "mesh_sha256": sha256_file(mesh_file),
        "ne": mesh.ne, "nv": mesh.nv, "nedge": mesh.nedge,
        "order": args.order,
        "ndof": fes.ndof,
        "coarse": args.coarse,
        "inverse": args.inverse if args.coarse == "direct" else None,
        "cycle_type": args.cycle_type, "num_smooth": args.num_smooth,
        "cycles": args.cycles,
        "lean_coarse": args.lean_coarse,
        "coarse_eps": args.coarse_eps,
        "system_eps": args.eps,
        "edge_wirebasket": bool(args.edge_wirebasket),
        "face_dofs_returned_to_interface": face_wirebasket,
        "t_coupling_reset": round(t_coupling, 3),
        "threads": args.threads,
        "iterations": iters,
        "converged": bool(iters < args.maxiter and true_res <= 10 * args.tol),
        "true_residual": true_res,
        "t_assemble_and_bddc_setup": round(t_assemble, 3),
        "t_solve": round(t_solve, 3),
        "ms_per_iter": round(1000 * t_solve / max(iters, 1), 2),
        "t_total_wall": round(t_total, 3),
        "peak_memory_mb": round(peak_memory_mb(), 1),
        "conductor_loss_w": float(loss),
        "magnetic_energy_j": float(energy),
        "coarse_stats": coarse_stats,
    }
    print(json.dumps(result, indent=2))
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mesh")
    p.add_argument("--mesh-dir")
    p.add_argument("--ssn-dir", help="directory holding sparsesolv_ngsolve.pyd "
                   "(default: radia.sparsesolv_ngsolve)")
    p.add_argument("--order", type=int, default=2)
    p.add_argument("--coarse", default="direct")
    p.add_argument("--inverse", default="sparsecholesky")
    p.add_argument("--cycle-type", type=int, default=1)
    p.add_argument("--num-smooth", type=int, default=1)
    p.add_argument("--lean-coarse", type=int, default=1,
                   help="auxiliary AMGs stop where coarsening stalls and solve the coarsest level densely")
    p.add_argument("--cycles", type=int, default=1,
                   help="AMS cycles (stationary steps) per wirebasket solve")
    p.add_argument("--coarse-eps", type=float, default=0.0,
                   help="relative diagonal shift of the AMS surrogate; the system's "
                   "own eps*nu mass already regularizes it")
    p.add_argument("--eps", type=float, default=1e-6)
    p.add_argument("--edge-wirebasket", action="store_true",
                   help="return non-edge wirebasket dofs to the interface")
    p.add_argument("--tol", type=float, default=1e-8)
    p.add_argument("--maxiter", type=int, default=2000)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--label", default="")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = p.parse_args(argv)

    import ngsolve
    result = run(args)
    out = {
        "schema": "radia.validation.sparsesolv-hiruma-bddc-coarse.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "label": args.label,
        "physics": {"frequency_hz": freq, "sigma_cu": sigma_cu,
                    "mu_r_core": mu_r_core, "tol": args.tol,
                    "maxiter": args.maxiter},
        "environment": {"platform": platform.platform(),
                        "processor": platform.processor(),
                        "python": platform.python_version(),
                        "ngsolve": ngsolve.__version__,
                        "cpu_count": os.cpu_count()},
        "argv": sys.argv[1:] if argv is None else list(argv),
        "result": result,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"saved {args.output}")


if __name__ == "__main__":
    main()
