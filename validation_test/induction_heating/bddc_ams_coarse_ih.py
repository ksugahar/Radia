"""BDDC with the AMS wirebasket solver on the IH FEM-SIBC system.

The production FEM-SIBC solver (calc_fem_kelvin) uses, at HCurl order 2 and
3, BDDC with ``coarsetype="sparsesolv_ams"`` and COCR.  The AMS coarse
solver was tuned on a volumetric conductor (``hiruma/bench_bddc_coarse.py``:
cycles >= 2, an edge-only wirebasket); the SIBC system is different -- the
workpiece is a hole, air only, with an impedance term on its surface -- so
the coarse settings are measured here on that system.

Each run solves the production problem through ``calc_fem_kelvin.solve_fem``
in its own process.  The coarse settings are applied from outside by
wrapping ``ngsolve.Preconditioner`` (``coarseflags={"cycles": k}``, and for
``--edge-wirebasket`` the higher-order HCurl wirebasket dofs returned to the
interface), so the production code carries no study switch.  Every setting
is compared with the SparseCholesky solve of the same system.

Cases (150 kHz, copper):
    tube    tube with an axial bore in a coaxial loop
    cross   cylinder with a radial through hole in front of a side loop

Output: results/bddc_ams_coarse_ih.json (one record per run).

    python validation_test/induction_heating/bddc_ams_coarse_ih.py
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))
sys.path.insert(0, HERE)

R, H, BORE, LOOP_R = 0.025, 0.025, 0.008, 0.035
CROSS, SIDE_R, SIDE_X = 0.005, 0.012, 0.036
R_AIR, RING_R, RING_H = 0.15, 0.052, 0.06
MESHES = {"coarse": (0.004, 0.006), "fine": (0.0025, 0.004)}


def build_vol(case, mesh_label, path, order):
    import ngsolve
    import fem_sibc_geometry as G
    mp, mr = MESHES[mesh_label]
    if case == "tube":
        make = lambda: G.cylinder_part(R, H, bore=BORE)  # noqa: E731
    else:
        make = lambda: G.cylinder_part(R, H, cross_hole=CROSS,  # noqa: E731
                                       hole_maxh=0.001)
    with ngsolve.TaskManager():
        G.air_mesh(make, r_air=R_AIR, ring_r=RING_R, ring_h=RING_H, maxh_part=mp,
                   maxh_ring=mr, maxh_far=0.04, order=order).ngmesh.Save(path)


def coil(case):
    import fem_sibc_geometry as G
    if case == "tube":
        return G.loop_paths((0, 0, 0), (0, 0, 1), LOOP_R)
    return G.loop_paths((SIDE_X, 0, 0), (1, 0, 0), SIDE_R)


def install_coarse_settings(cycles, edge_wirebasket):
    """Wrap ngsolve.Preconditioner so the AMS wirebasket solver gets the
    requested settings; returns a record of what was applied."""
    import ngsolve
    from ngsolve import COUPLING_TYPE
    real = ngsolve.Preconditioner
    applied = {"calls": 0, "face_dofs_returned": 0}

    def wrapped(bf, name, **kw):
        if name == "bddc" and kw.get("coarsetype") == "sparsesolv_ams":
            if "coarseflags" in kw:
                raise RuntimeError("the production call now sets coarseflags; "
                                   "the study wrapper would override them")
            if edge_wirebasket:
                fes = bf.space
                nedge = fes.mesh.nedge
                for dof in range(nedge, fes.ndof):
                    if fes.CouplingType(dof) == COUPLING_TYPE.WIREBASKET_DOF:
                        fes.SetCouplingType(dof, COUPLING_TYPE.INTERFACE_DOF)
                        applied["face_dofs_returned"] += 1
            kw["coarseflags"] = {"cycles": int(cycles)}
            applied["calls"] += 1
        return real(bf, name, **kw)

    ngsolve.Preconditioner = wrapped
    return applied


def peak_memory_mb():
    import psutil
    mem = psutil.Process(os.getpid()).memory_info()
    return getattr(mem, "peak_wset", mem.rss) / 2**20


def single(a):
    """One solve in this process; prints one JSON record."""
    import ngsolve
    import calc_fem_kelvin as K
    from radia.em_material import EMMaterial
    ngsolve.SetNumThreads(a.threads)
    applied = None
    if a.solver == "bddc":
        applied = install_coarse_settings(a.cycles, a.edge_wirebasket)
    t0 = time.perf_counter()
    res = K.solve_fem(vol_file=a.vol, fes_order=a.order, frequency=a.frequency,
                      mat=EMMaterial(name="copper", sigma=5.8e7, mu_r=1.0),
                      half_thickness=R, solver=a.solver, nthreads=a.threads,
                      filaments=([coil(a.case)], [100.0]))
    wall = time.perf_counter() - t0
    keys = ("P_total", "L", "error", "linear_solver", "linear_true_relative_residual",
            "linear_krylov_iterations", "linear_residual_corrections",
            "linear_residual_before_corrections", "ndof")
    rec = {k: res.get(k) for k in keys}
    rec.update({"case": a.case, "mesh": a.mesh, "order": a.order,
                "solver": a.solver, "cycles": a.cycles if a.solver == "bddc" else None,
                "edge_wirebasket": bool(a.edge_wirebasket) if a.solver == "bddc" else None,
                "coarse_applied": applied, "wall_s": round(wall, 3),
                "peak_memory_mb": round(peak_memory_mb(), 1),
                "threads": a.threads})
    print("RECORD " + json.dumps(rec, default=str), flush=True)


def compare_results(records):
    """Check convergence and same-mesh physical outputs against direct solves."""
    import math

    def solved(record):
        residual = record.get("linear_true_relative_residual")
        return (record.get("process_exit", 0) == 0 and not record.get("error")
                and "timeout_s" not in record and residual is not None
                and math.isfinite(residual) and residual <= 1e-7
                and all(record.get(key) is not None and math.isfinite(record[key])
                        for key in ("P_total", "L")))

    reference = {(r["case"], r["mesh"], r["order"]): r for r in records
                 if r["solver"] == "sparsecholesky" and solved(r)}
    checks = []
    for record in records:
        key = (record["case"], record["mesh"], record["order"])
        direct = reference.get(key)
        errors = {}
        passed = solved(record) and direct is not None
        if passed:
            for name in ("P_total", "L"):
                errors[name] = abs(record[name]-direct[name]) / max(abs(direct[name]), 1e-300)
            passed = all(value <= 1e-6 for value in errors.values())
        checks.append(dict(case=key[0], mesh=key[1], order=key[2],
                           solver=record["solver"], cycles=record.get("cycles"),
                           edge_wirebasket=record.get("edge_wirebasket"),
                           passed=bool(passed), relative_errors=errors))
    return dict(all_passed=bool(checks) and all(c["passed"] for c in checks),
                true_residual_limit=1e-7, parity_limit=1e-6, checks=checks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--vol")
    ap.add_argument("--case", default="tube")
    ap.add_argument("--mesh", default="coarse")
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--solver", default="bddc")
    ap.add_argument("--cycles", type=int, default=1)
    ap.add_argument("--edge-wirebasket", action="store_true")
    ap.add_argument("--frequency", type=float, default=150000.0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cases", default="tube,cross")
    ap.add_argument("--meshes", default="coarse,fine")
    ap.add_argument("--orders", default="2,3")
    ap.add_argument("--cycles-list", default="1,2,3,4,6")
    ap.add_argument("--timeout", type=float, default=7200.0)
    ap.add_argument("--work", default=os.path.join(
        os.environ.get("TEMP", "C:\\temp"), "bddc_ams_coarse_ih"))
    ap.add_argument("--out", default=os.path.join(
        HERE, "results", "bddc_ams_coarse_ih.json"))
    ap.add_argument("--source-commit", default="")
    a = ap.parse_args()
    if a.single:
        single(a)
        return
    from _provenance import provenance
    prov = provenance(sys.argv, a.source_commit, [
        os.path.abspath(__file__), os.path.join(HERE, "fem_sibc_geometry.py"),
        os.path.join(ROOT, "src", "radia", "panels", "calc_fem_kelvin.py")])
    os.makedirs(a.work, exist_ok=True)
    records = []
    out = {"case_settings": {"frequency_Hz": a.frequency, "meshes": MESHES,
                             "threads": a.threads},
           "records": records, "provenance": prov}

    def save():
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1)

    for case in a.cases.split(","):
        for mesh in a.meshes.split(","):
            vol = os.path.join(a.work, f"{case}_{mesh}.vol")
            # A work directory may contain a mesh from another source/order.
            # Generate this campaign's mesh and use it for every setting.
            build_vol(case, mesh, vol, max(int(v) for v in a.orders.split(",")))
            for order in [int(v) for v in a.orders.split(",")]:
                runs = [("sparsecholesky", 0, False)] + [
                    ("bddc", k, e) for e in (False, True)
                    for k in [int(v) for v in a.cycles_list.split(",")]]
                for solver, k, edge in runs:
                    cmd = [sys.executable, "-u", os.path.abspath(__file__),
                           "--single", "--vol", vol, "--case", case, "--mesh", mesh,
                           "--order", str(order), "--solver", solver,
                           "--cycles", str(k), "--frequency", str(a.frequency),
                           "--threads", str(a.threads)]
                    if edge:
                        cmd.append("--edge-wirebasket")
                    t0 = time.perf_counter()
                    log_stem = f"{case}_{mesh}_p{order}_{solver}_c{k}_e{int(edge)}"
                    stdout_path = Path(a.work) / (log_stem + ".stdout.log")
                    stderr_path = Path(a.work) / (log_stem + ".stderr.log")
                    try:
                        cp = subprocess.run(cmd, capture_output=True, text=True,
                                            timeout=a.timeout)
                        stdout_path.write_text(cp.stdout, encoding="utf-8")
                        stderr_path.write_text(cp.stderr, encoding="utf-8")
                        line = [ln for ln in cp.stdout.splitlines()
                                if ln.startswith("RECORD ")]
                        if line:
                            rec = json.loads(line[-1][7:])
                            rec["process_exit"] = cp.returncode
                        else:
                            rec = {"case": case, "mesh": mesh, "order": order,
                                   "solver": solver, "cycles": k,
                                   "edge_wirebasket": edge,
                                   "process_exit": cp.returncode,
                                   "stderr_tail": cp.stderr[-2000:]}
                    except subprocess.TimeoutExpired as error:
                        for path, data in ((stdout_path, error.stdout), (stderr_path, error.stderr)):
                            if isinstance(data, bytes):
                                data = data.decode("utf-8", errors="replace")
                            path.write_text(data or "", encoding="utf-8")
                        rec = {"case": case, "mesh": mesh, "order": order,
                               "solver": solver, "cycles": k,
                               "edge_wirebasket": edge, "timeout_s": a.timeout}
                    rec["stdout_log"] = stdout_path.name
                    rec["stderr_log"] = stderr_path.name
                    rec["process_wall_s"] = round(time.perf_counter() - t0, 3)
                    records.append(rec)
                    print(json.dumps({k2: rec.get(k2) for k2 in (
                        "case", "mesh", "order", "solver", "cycles",
                        "edge_wirebasket", "error", "linear_krylov_iterations",
                        "linear_residual_corrections", "wall_s",
                        "peak_memory_mb", "process_exit")}), flush=True)
                    save()
    out["acceptance"] = compare_results(records)
    save()
    print("wrote", a.out)
    if not out["acceptance"]["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
