"""Complex eddy-current cases for the BDDC wirebasket solver at HCurl order 2.

Each case is solved with NGSolve BDDC + COCR, once with the direct wirebasket
inverse and once with ``coarsetype="sparsesolv_ams"`` (see
``hiruma/bench_bddc_coarse.py`` for the Hiruma coil).  All cases use the
reduced-field form A = A_s + A_r with homogeneous Dirichlet data for A_r on
the outer boundary and a small ``eps*nu`` mass (the sigma=0 air otherwise
leaves the lowest-order gradients in the kernel):

    sphere  Cu sphere (a = 5 mm) in a uniform AC field, air box.  Reference:
            the analytic (Smythe) induced dipole m_z; loss = (w/2)|Im m_z| B0.
    disk    thin disk, mu_r = 100, sigma = 1e7, in a uniform axial H0, air
            box.  Reference: stored axisymmetric volume-averaged Bz/(mu0 H0)
            (validation_test/vim_coupled/
            results_magnetic_conductor_disk_adjudication.json).
    plate   Cu plate (17.72 x 17.72 x 2 mm) on a conductor-only mesh with
            A_r = 0 on its surface (the A1 sweep setting).  Reference: the
            direct run of the same system.

Usage:
    python bddc_coarse_cases.py sphere --freq 700 --coarse direct
    python bddc_coarse_cases.py disk --freq 1e4 --coarse sparsesolv_ams --cycles 4
"""

import argparse
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

mu0 = 4e-7 * np.pi
DISK_REFERENCE = {100.0: 1.09161 - 0.00417j, 10000.0: 0.95172 - 0.11536j}


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


def sphere_case(args):
    from ngsolve import CF, Mesh, x, y
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, Sphere
    a, half = 5e-3, 0.030
    cond = Sphere(Pnt(0, 0, 0), a)
    cond.mat("cond")
    cond.maxh = args.hcond or 6e-4
    box = Box(Pnt(-half, -half, -half), Pnt(half, half, half))
    box.faces.name = "outer"
    air = box - cond
    air.mat("air")
    mesh = Mesh(OCCGeometry(Glue([cond, air])).GenerateMesh(maxh=args.hair or 4e-3))
    mesh.Curve(2)
    sigma_c = 5.8e7
    omega = 2 * np.pi * args.freq
    b0 = 1.0
    delta = np.sqrt(2 / (omega * mu0 * sigma_c))
    xa = (1 + 1j) * a / delta
    m_z = -2 * np.pi * a**3 * (b0 / mu0) * (1 - 3 / xa / np.tanh(xa) + 3 / xa**2)
    return {
        "mesh": mesh, "sigma": {"cond": sigma_c}, "mu_r": {},
        "a_s": CF((-0.5 * b0 * y, 0.5 * b0 * x, 0)), "curl_a_s": CF((0, 0, b0)),
        "omega": omega, "conductor": "cond",
        "reference": {"loss_w": 0.5 * omega * abs(m_z.imag) * b0,
                      "a_over_delta": a / delta, "source": "Smythe dipole moment"},
    }


def disk_case(args):
    from ngsolve import CF, Mesh, x, y
    from netgen.occ import Box, Cylinder, Glue, OCCGeometry, Pnt, Z
    disk = Cylinder(Pnt(0, 0, -0.25e-3), Z, r=0.010, h=0.0005)
    disk.mat("cond")
    disk.maxh = args.hcond or 1e-3
    box = Box(Pnt(-0.04, -0.04, -0.02), Pnt(0.04, 0.04, 0.02))
    box.faces.name = "outer"
    air = box - disk
    air.mat("air")
    mesh = Mesh(OCCGeometry(Glue([disk, air])).GenerateMesh(maxh=args.hair or 6e-3))
    h0 = 1.0
    b = mu0 * h0
    return {
        "mesh": mesh, "sigma": {"cond": 1e7}, "mu_r": {"cond": 100.0},
        "a_s": CF((-0.5 * b * y, 0.5 * b * x, 0)), "curl_a_s": CF((0, 0, b)),
        "omega": 2 * np.pi * args.freq, "conductor": "cond", "bz_normalize": b,
        "reference": {"bz_average": DISK_REFERENCE.get(float(args.freq)),
                      "source": "axisymmetric Q2 (stored)"},
    }


def plate_case(args):
    from ngsolve import CF, Mesh, x, y
    from netgen.occ import Box, OCCGeometry, Pnt
    half = 17.72e-3 / 2
    plate = Box(Pnt(-half, -half, -1e-3), Pnt(half, half, 1e-3))
    plate.mat("cond")
    plate.faces.name = "outer"
    mesh = Mesh(OCCGeometry(plate).GenerateMesh(maxh=args.hcond or 5e-4, grading=0.5))
    b0 = 1.0
    return {
        "mesh": mesh, "sigma": {"cond": 5.8e7}, "mu_r": {},
        "a_s": CF((-0.5 * b0 * y, 0.5 * b0 * x, 0)), "curl_a_s": CF((0, 0, b0)),
        "omega": 2 * np.pi * args.freq, "conductor": "cond",
        "reference": {"source": "direct wirebasket run of the same system"},
    }


CASES = {"sphere": sphere_case, "disk": disk_case, "plate": plate_case}


def run(args):
    import ngsolve
    from ngsolve import (COUPLING_TYPE, BilinearForm, GridFunction, HCurl,
                         InnerProduct, Integrate, LinearForm, Preconditioner,
                         TaskManager, curl, dx)
    ssn = import_sparsesolv(args.ssn_dir)
    ngsolve.SetNumThreads(args.threads)
    t_start = time.perf_counter()
    with TaskManager():
        case = CASES[args.case](args)
    t_mesh = time.perf_counter() - t_start
    mesh, omega = case["mesh"], case["omega"]
    nu0 = 1 / mu0
    nu = mesh.MaterialCF({k: nu0 / v for k, v in case["mu_r"].items()}, default=nu0)
    sigma = mesh.MaterialCF(case["sigma"], default=0)
    cond = case["conductor"]

    fes = HCurl(mesh, order=args.order, nograds=True, dirichlet="outer", complex=True)
    face_wirebasket = 0
    if args.edge_wirebasket:
        for dof in range(mesh.nedge, fes.ndof):
            if fes.CouplingType(dof) == COUPLING_TYPE.WIREBASKET_DOF:
                fes.SetCouplingType(dof, COUPLING_TYPE.INTERFACE_DOF)
                face_wirebasket += 1
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu * curl(u) * curl(v) * dx
    a += 1j * omega * sigma * u * v * dx(cond)
    if args.eps > 0:
        a += args.eps * nu * u * v * dx
    flags = {"coarsetype": args.coarse}
    if args.coarse == "direct":
        flags["inverse"] = args.inverse
    else:
        flags["coarseflags"] = {"cycles": args.cycles, "lean_coarse": args.lean_coarse}
    pre = Preconditioner(a, "bddc", **flags)
    # Reduced field: a(A_r, v) = -a(A_s, v); curl-curl of the uniform-curl A_s
    # only acts where nu differs from nu0 (A_r vanishes tangentially outside).
    f = LinearForm(fes)
    f += -1j * omega * sigma * case["a_s"] * v * dx(cond)
    if case["mu_r"]:
        f += -(nu - nu0) * case["curl_a_s"] * curl(v) * dx(cond)

    with TaskManager():
        t0 = time.perf_counter()
        a.Assemble()
        f.Assemble()
        t_assemble = time.perf_counter() - t0
        gfu = GridFunction(fes)
        t0 = time.perf_counter()
        inv = ssn.COCRSolver(a.mat, pre, freedofs=fes.FreeDofs(), maxiter=args.maxiter,
                             tol=args.tol, printrates=False)
        gfu.vec.data = inv * f.vec
        t_solve = time.perf_counter() - t0
        iters = inv.iterations
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * gfu.vec
        mask = np.array(fes.FreeDofs(), dtype=bool)
        true_res = float(np.linalg.norm(r.FV().NumPy()[mask]) / np.linalg.norm(f.vec.FV().NumPy()[mask]))
        a_tot = gfu + case["a_s"]
        loss = Integrate(0.5 * sigma * omega**2 * InnerProduct(a_tot, a_tot).real, mesh,
                         definedon=mesh.Materials(cond))
        energy = Integrate(nu * InnerProduct(curl(gfu), curl(gfu)).real, mesh)
        quantities = {"conductor_loss_w": float(loss), "reduced_field_energy": float(energy)}
        if "bz_normalize" in case:
            vol = Integrate(1.0, mesh, definedon=mesh.Materials(cond))
            bz = Integrate((curl(gfu)[2] + case["curl_a_s"][2]), mesh,
                           definedon=mesh.Materials(cond)) / (vol * case["bz_normalize"])
            quantities["bz_average"] = [float(bz.real), float(bz.imag)]
    stats = dict(ssn.AMSCoarseStats()) if args.coarse != "direct" else {}
    stats.pop("amg_levels", None)
    ref = dict(case["reference"])
    if ref.get("loss_w"):
        ref["loss_relative_error"] = abs(quantities["conductor_loss_w"] - ref["loss_w"]) / ref["loss_w"]
    if ref.get("bz_average") is not None:
        rb = ref["bz_average"]
        mb = complex(*quantities["bz_average"])
        ref["bz_average"] = [rb.real, rb.imag]
        ref["bz_relative_error"] = abs(mb - rb) / abs(rb)
    return {
        "case": args.case, "freq_hz": args.freq, "ne": mesh.ne, "nv": mesh.nv, "nedge": mesh.nedge,
        "order": args.order, "ndof": fes.ndof, "coarse": args.coarse,
        "inverse": args.inverse if args.coarse == "direct" else None,
        "cycles": args.cycles, "lean_coarse": args.lean_coarse, "system_eps": args.eps,
        "edge_wirebasket": bool(args.edge_wirebasket),
        "face_dofs_returned_to_interface": face_wirebasket, "threads": args.threads,
        "iterations": iters, "converged": bool(iters < args.maxiter and true_res <= 10 * args.tol),
        "true_residual": true_res, "t_mesh": round(t_mesh, 3),
        "t_assemble_and_bddc_setup": round(t_assemble, 3), "t_solve": round(t_solve, 3),
        "t_total_wall": round(time.perf_counter() - t_start, 3),
        "peak_memory_mb": round(peak_memory_mb(), 1), **quantities,
        "reference": ref, "coarse_stats": stats,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("case", choices=sorted(CASES))
    p.add_argument("--freq", type=float, required=True)
    p.add_argument("--hcond", type=float, default=None)
    p.add_argument("--hair", type=float, default=None)
    p.add_argument("--ssn-dir")
    p.add_argument("--order", type=int, default=2)
    p.add_argument("--coarse", default="direct")
    p.add_argument("--inverse", default="sparsecholesky")
    p.add_argument("--cycles", type=int, default=4)
    p.add_argument("--lean-coarse", type=int, default=1)
    p.add_argument("--eps", type=float, default=1e-6)
    p.add_argument("--edge-wirebasket", action="store_true")
    p.add_argument("--tol", type=float, default=1e-8)
    p.add_argument("--maxiter", type=int, default=3000)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--label", default="")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(argv)
    import ngsolve
    result = run(args)
    out = {
        "schema": "radia.validation.sparsesolv-bddc-coarse-cases.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(), "label": args.label,
        "environment": {"platform": platform.platform(), "python": platform.python_version(),
                        "ngsolve": ngsolve.__version__, "cpu_count": os.cpu_count()},
        "argv": sys.argv[1:] if argv is None else list(argv), "result": result,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("case", "freq_hz", "ndof", "coarse", "iterations",
                                             "converged", "t_assemble_and_bddc_setup", "t_solve",
                                             "peak_memory_mb", "conductor_loss_w")}
                     | {"reference": result["reference"], "bz": result.get("bz_average")}))


if __name__ == "__main__":
    main()
