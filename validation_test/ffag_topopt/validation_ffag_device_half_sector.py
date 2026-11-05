#!/usr/bin/env python
"""Device half sector against the device ring: curved Q2 annular HEX, midplane mirror plus rotations, linear and saturating.

One sector of the upper half stands for the whole ring through `image="-z"`, `image_cyclic` and the identified cut faces.
The mesh is deformed by a map that keeps z = 0 and the rotation group, so both models hold the same curved geometry. The field
runs along z, which is the field of an FFAG dipole and the one the mirror `-z` assumes. The reference is the device route
with every image removed, on the closed ring with the same cells above and below the midplane, so this checks the reduction
and not the accuracy of the discretization. A half sector without images shows that a no-op image route could not pass.
Limits are PROPOSED, not adopted. Timings are records, not gates.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parents[0]/"hmatrix_cuda"))

PROPOSED_MEAN_REL = 1e-6
SEAM = ("periodic_min", "periodic_max")
BH_TABLE = ((0.0, 0.0), (200.0, 0.8), (1000.0, 1.4), (5000.0, 1.7), (50000.0, 1.95), (500000.0, 2.6))


def _contract(boundaries):
    return dict(schema="radia.vol-label-contract.v1", application="ffag_half_sector",
                metadata=dict(length_unit="m", unit_source="Caller declaration; the VOL file has no unit tag"),
                strict_labels=True, required=dict(materials=["yoke"], boundaries=list(boundaries)),
                allowed=dict(materials=["yoke"], boundaries=list(boundaries)), quality=dict(min_curve_order=1))


def build_yoke_mesh(*, fold, cells, radial, levels, ring, inner_radius, outer_radius):
    """Annular HEX yoke with `cells` per sector, `radial` cells across and its height split at `levels`."""
    import ngsolve as ng
    from netgen.csg import Pnt
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, IdentificationType, Mesh, MeshPoint

    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "yoke")
    skin = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    mesh.SetBCName(0, "skin")
    if not ring:
        low, high = mesh.Add(FaceDescriptor(surfnr=2, domin=1, domout=0, bc=2)), mesh.Add(FaceDescriptor(surfnr=3, domin=1, domout=0, bc=3))
        mesh.SetBCName(1, SEAM[0])
        mesh.SetBCName(2, SEAM[1])
    total = fold*cells if ring else cells
    stations = total if ring else total + 1
    step = 2.0*math.pi/(fold*cells)
    points = {}
    for it, ir, iz in np.ndindex(stations, radial + 1, len(levels)):
        radius = inner_radius + (outer_radius - inner_radius)*ir/radial
        points[it, ir, iz] = mesh.Add(MeshPoint(Pnt(radius*math.cos(it*step), radius*math.sin(it*step), levels[iz])))
    faces = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6), (3, 0, 4, 7), (1, 2, 6, 5))
    for it, ir, iz in np.ndindex(total, radial, len(levels) - 1):
        following = (it + 1) % stations
        corners = [points[it, ir, iz], points[it, ir + 1, iz], points[following, ir + 1, iz], points[following, ir, iz],
                   points[it, ir, iz + 1], points[it, ir + 1, iz + 1], points[following, ir + 1, iz + 1], points[following, ir, iz + 1]]
        mesh.Add(Element3D(1, corners))
        for index, face in enumerate(faces):
            outer = ((index == 0 and iz == 0) or (index == 1 and iz == len(levels) - 2)
                     or (index == 4 and ir == 0) or (index == 5 and ir == radial - 1))
            if outer:
                mesh.Add(Element2D(skin, [corners[k] for k in face]))
            elif not ring and index == 2 and it == 0:
                mesh.Add(Element2D(low, [corners[k] for k in face]))
            elif not ring and index == 3 and it == total - 1:
                mesh.Add(Element2D(high, [corners[k] for k in face]))
    if not ring:
        for ir, iz in np.ndindex(radial + 1, len(levels)):
            mesh.AddPointIdentification(points[0, ir, iz], points[total, ir, iz], identnr=1, type=IdentificationType.PERIODIC)
    return ng.Mesh(mesh)


def _mean_z(field, mesh):
    import ngsolve as ng
    return float(ng.Integrate(field[2], mesh))/float(ng.Integrate(ng.CF(1.0), mesh))


def _solve(path, checked, args, *, compression, nonlinear, **images):
    import ngsolve as ng
    from radia.hdiv_gpu import PreparedHDivLinear, PreparedHDivNonlinear
    from validation_ffag_connected_cyclic_yoke import apply_cyclic_q2_geometry

    source = ng.CF((0.0, 0.0, args.applied_field_am))
    geometry = dict(curve_order=2, q2_deformation=lambda mesh: apply_cyclic_q2_geometry(mesh, args))
    started = time.perf_counter()
    if nonlinear:
        device = PreparedHDivNonlinear.from_bh_table(path, np.array(BH_TABLE), backend="gpu", order=args.hdiv_order, **geometry, **checked, **images)
    else:
        device = PreparedHDivLinear(path, backend="gpu", order=args.hdiv_order, gram_compression=compression, **geometry, **checked, **images)
    with device:
        setup_s = time.perf_counter() - started
        started = time.perf_counter()
        if nonlinear:
            solved = device.solve(source, tolerance=args.nonlinear_tol, linear_tolerance=args.device_tol)
        else:
            solved = device.solve(source, 1.0/(args.mu_r - 1.0), tolerance=args.device_tol, max_iterations=args.max_iterations)
        solve_s = time.perf_counter() - started
        magnetization = ng.GridFunction(device.fespace)
        magnetization.vec.FV().NumPy()[:] = solved["x"]
        mean = _mean_z(magnetization, device.fespace.mesh)
        ledger = device.stats()["assembly"]["transfer_ledger"]
    return dict(mean_z=mean, dofs=int(device.dofs), setup_s=setup_s, solve_s=solve_s, transfer_ledger=ledger,
                iterations=int(solved.get("newton_iterations", solved.get("iterations", 0))), relative_residual=float(solved["relative_residual"]))


def _case(case, args, contracts):
    from device_scale import checked_inputs

    cells, radial, height = case
    root = args.work/f"c{cells}-r{radial}-h{height}"
    root.mkdir(parents=True, exist_ok=False)
    upper = [args.half_height*k/height for k in range(height + 1)]
    levels = [-z for z in reversed(upper[1:])] + upper
    geometry = dict(fold=args.fold, cells=cells, radial=radial, inner_radius=args.inner_radius, outer_radius=args.outer_radius)
    paths = dict(ring=root/"ring.vol", half=root/"half.vol")
    build_yoke_mesh(levels=levels, ring=True, **geometry).ngmesh.Save(str(paths["ring"]))
    build_yoke_mesh(levels=upper, ring=False, **geometry).ngmesh.Save(str(paths["half"]))
    checked = {name: checked_inputs(paths[name], contracts[name], args.receipt_python) for name in paths}
    images = dict(image="-z", image_cyclic=args.fold, cyclic_periodic_boundaries=SEAM)
    rows, gates = [], {}
    for label, compression, nonlinear in args.routes:
        half = _solve(paths["half"], checked["half"], args, compression=compression, nonlinear=nonlinear, **images)
        ring = _solve(paths["ring"], checked["ring"], args, compression=compression, nonlinear=nonlinear) if args.fold*(cells*radial*height*2) <= args.ring_cell_limit else None
        lone = None if nonlinear else _solve(paths["half"], checked["half"], args, compression=compression, nonlinear=nonlinear)
        row = dict(route=label, cells=cells, radial=radial, height=height, half=half, ring=ring, lone=lone)
        if ring is not None:
            row["half_vs_ring"] = abs(half["mean_z"] - ring["mean_z"])/abs(ring["mean_z"])
            gates[f"c{cells}-r{radial}-h{height}-{label}-half-matches-ring"] = bool(row["half_vs_ring"] <= PROPOSED_MEAN_REL)
            if lone is not None:
                row["lone_vs_ring"] = abs(lone["mean_z"] - ring["mean_z"])/abs(ring["mean_z"])
                gates[f"c{cells}-r{radial}-h{height}-{label}-lone-discriminates"] = bool(row["lone_vs_ring"] > 1e3*PROPOSED_MEAN_REL)
        print(json.dumps(dict(route=label, cells=cells, radial=radial, height=height, half_dofs=half["dofs"], half=half["mean_z"], half_s=half["setup_s"] + half["solve_s"],
                              ring_dofs=None if ring is None else ring["dofs"], ring=None if ring is None else ring["mean_z"], ring_s=None if ring is None else ring["setup_s"] + ring["solve_s"],
                              error=row.get("half_vs_ring"), lone_error=row.get("lone_vs_ring"))), flush=True)
        rows.append(row)
    return rows, gates


def run_validation(args):
    import ngsolve as ng

    args.work.mkdir(parents=True, exist_ok=False)
    contracts = dict(ring=args.work/"labels-ring.json", half=args.work/"labels-half.json")
    contracts["ring"].write_text(json.dumps(_contract(["skin"]), indent=2) + "\n", encoding="utf-8")
    contracts["half"].write_text(json.dumps(_contract(["skin", *SEAM]), indent=2) + "\n", encoding="utf-8")
    ng.SetNumThreads(args.threads)
    args.routes = [(name, name.split("-")[0], name.startswith("nonlinear")) for name in args.route]
    started = time.perf_counter()
    rows, gates = [], {}
    for case in args.cases:
        case_rows, case_gates = _case(tuple(int(value) for value in case.split(",")), args, contracts)
        rows += case_rows
        gates.update(case_gates)
    return dict(schema="radia.ffag-device-half-sector/v1", status="pass" if gates and all(gates.values()) else "fail", adoption="proposed", accepted=False,
                proposed_limits=dict(mean=PROPOSED_MEAN_REL), gates=gates, rows=rows,
                settings={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
                source_sha=os.environ.get("RADIA_SOURCE_SHA"), ngsolve_version=ng.__version__, wall_s=time.perf_counter() - started,
                runtime=dict(platform_class=platform.system(), python=platform.python_version()))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True, help="New folder for the meshes and their receipts")
    parser.add_argument("--receipt-python", type=Path, required=True, help="Interpreter of the mesh checker")
    parser.add_argument("--cases", nargs="+", default=["1,1,1", "2,2,1"], help="cells per sector,radial cells,height cells above the midplane")
    parser.add_argument("--route", nargs="+", choices=("dense", "aca", "nonlinear"), default=["dense", "aca", "nonlinear"])
    parser.add_argument("--ring-cell-limit", type=int, default=2000, help="Skip the closed ring above this cell count and record the half alone")
    parser.add_argument("--fold", type=int, default=12)
    parser.add_argument("--hdiv-order", type=int, choices=(1, 2), default=2)
    parser.add_argument("--geometry-order", type=int, choices=(2,), default=2)
    parser.add_argument("--inner-radius", type=float, default=0.030)
    parser.add_argument("--outer-radius", type=float, default=0.050)
    parser.add_argument("--half-height", type=float, default=0.004)
    parser.add_argument("--radial-deformation-scale", type=float, default=2.5)
    parser.add_argument("--axial-deformation-scale", type=float, default=1.5)
    parser.add_argument("--mu-r", type=float, default=100.0)
    parser.add_argument("--applied-field-am", type=float, default=1.0e5)
    parser.add_argument("--device-tol", type=float, default=1.0e-10)
    parser.add_argument("--nonlinear-tol", type=float, default=1.0e-9)
    parser.add_argument("--max-iterations", type=int, default=4000)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    report = run_validation(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report["settings"]["routes"] = [list(row) for row in args.routes]
    args.output.write_text(json.dumps(report, indent=1, allow_nan=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps(dict(status=report["status"], gates=report["gates"]), indent=1))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
