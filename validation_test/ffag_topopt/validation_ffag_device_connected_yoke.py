#!/usr/bin/env python
"""Device connected-yoke sector against the explicit ring: curved Q2 BDM2 HEX, periodic trace plus cyclic images, JSON record.

The contract of ``validation_ffag_connected_cyclic_yoke.py`` on the device route. The yoke crosses both azimuthal cuts
with nonzero normal magnetization, so the sector needs the identified trace, the paired faces out of the charge skin and
the rotated images. NGSolve supplies only non-matrix inputs here: the ``Compress(Periodic(HDiv))`` numbering, the periodic
wrapper's orientation and the deformed Q2 lattices. Every operator of a device solve is built and kept on the device, and
the transfer ledgers are recorded. The linear reference is the CPU equation of the explicit ring: the CPU charge map and
Gram with the mass on the device's material rule, which differs from the rule of the CPU linear solve on non-affine hosts.
The gates chain sector to device ring to that equation, the device field on the CPU ring's own state closes the chain to the
CPU field, and the rows against the CPU linear solve stay recorded beside the size of the rule difference, ungated. The CPU
nonlinear solve takes the device's rule, so the nonlinear sector is held to the CPU ring directly. On several cells the CPU
sector is recorded only, since the CPU image dispatch selects a far rule for touching seam pairs of a multi-cell sector. A
device sector without images shows that a no-op image route could not pass. The limits repeat the CPU route's numbers but
are PROPOSED for the device route, not adopted.
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

PROPOSED_MAXIMUM_MEAN_RELATIVE_ERROR = 1.0e-8
PROPOSED_MAXIMUM_FIELD_RELATIVE_TO_DRIVE = 2.0e-7
PROPOSED_MAXIMUM_STATE_RELATIVE_ERROR = 1.0e-8
SEAM = ("periodic_min", "periodic_max")
BH_TABLE = ((0.0, 0.0), (200.0, 0.8), (1000.0, 1.4), (5000.0, 1.7), (50000.0, 1.95), (500000.0, 2.6))
_PHYSICAL_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (1, 2, 6, 5), (3, 0, 4, 7))
_MIN_FACE, _MAX_FACE = (0, 1, 5, 4), (2, 3, 7, 6)


def _contract(boundaries):
    return dict(schema="radia.vol-label-contract.v1", application="ffag_connected_cyclic_yoke",
                metadata=dict(length_unit="m", unit_source="Caller declaration; the VOL file has no unit tag"),
                strict_labels=True, required=dict(materials=["yoke"], boundaries=list(boundaries)),
                allowed=dict(materials=["yoke"], boundaries=list(boundaries)), quality=dict(min_curve_order=1))


def build_yoke_mesh(*, fold, cells, full_ring, inner_radius, outer_radius, half_height):
    """Annular HEX yoke: `cells` per sector, either the closed ring or one sector with identified cut faces."""
    import ngsolve as ng
    from netgen.csg import Pnt
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, IdentificationType, Mesh, MeshPoint

    netmesh = Mesh(dim=3)
    netmesh.SetMaterial(1, "yoke")
    skin = netmesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    netmesh.SetBCName(0, "skin")
    if not full_ring:
        low = netmesh.Add(FaceDescriptor(surfnr=2, domin=1, domout=0, bc=2))
        high = netmesh.Add(FaceDescriptor(surfnr=3, domin=1, domout=0, bc=3))
        netmesh.SetBCName(1, SEAM[0])
        netmesh.SetBCName(2, SEAM[1])
    total = fold*cells if full_ring else cells
    stations = total if full_ring else total + 1
    step = 2.0*math.pi/(fold*cells)
    points = {}
    for it, ir, iz in np.ndindex(stations, 2, 2):
        radius = (inner_radius, outer_radius)[ir]
        points[it, ir, iz] = netmesh.Add(MeshPoint(Pnt(
            radius*math.cos(it*step), radius*math.sin(it*step), (-half_height, half_height)[iz])))
    for it in range(total):
        following = (it + 1) % stations
        element = [points[it, 0, 0], points[it, 1, 0], points[following, 1, 0], points[following, 0, 0],
                   points[it, 0, 1], points[it, 1, 1], points[following, 1, 1], points[following, 0, 1]]
        netmesh.Add(Element3D(1, element))
        for face in _PHYSICAL_FACES:
            netmesh.Add(Element2D(skin, [element[index] for index in face]))
        if not full_ring and it == 0:
            netmesh.Add(Element2D(low, [element[index] for index in _MIN_FACE]))
        if not full_ring and it == total - 1:
            netmesh.Add(Element2D(high, [element[index] for index in _MAX_FACE]))
    if not full_ring:
        for ir, iz in np.ndindex(2, 2):
            netmesh.AddPointIdentification(points[0, ir, iz], points[total, ir, iz], identnr=1, type=IdentificationType.PERIODIC)
    return ng.Mesh(netmesh)


def _applied(args):
    import ngsolve as ng
    radius = ng.sqrt(ng.x*ng.x + ng.y*ng.y)
    return ng.CF((-args.applied_field_am*ng.y/radius, args.applied_field_am*ng.x/radius, 0.0))


def _device_solve(path, checked, probes, args, compression, reference_state=None, **images):
    import ngsolve as ng
    from radia.hdiv_gpu import PreparedHDivLinear
    from validation_ffag_connected_cyclic_yoke import _mean_tangential_magnetization, apply_cyclic_q2_geometry

    started = time.perf_counter()
    with PreparedHDivLinear(path, backend="gpu", order=args.hdiv_order, gram_compression=compression, curve_order=2,
                            q2_deformation=lambda mesh: apply_cyclic_q2_geometry(mesh, args), **checked, **images) as device:
        setup_s = time.perf_counter() - started
        started = time.perf_counter()
        solved = device.solve(_applied(args), 1.0/(args.mu_r - 1.0), tolerance=args.device_tol, max_iterations=args.max_iterations)
        solve_s = time.perf_counter() - started
        field = device.field(solved["x"], probes)
        reference_field = None if reference_state is None else np.asarray(device.field(reference_state, probes)["field"])
        stats = device.stats()
        magnetization = ng.GridFunction(device.fespace)
        magnetization.vec.FV().NumPy()[:] = solved["x"]
        mean = _mean_tangential_magnetization(dict(gfM=magnetization), device.fespace.mesh)
    return dict(x=np.array(solved["x"]), mean_tangential=mean, field=np.asarray(field["field"]), dofs=int(device.dofs),
                field_of_reference_state=reference_field,
                iterations=int(solved["iterations"]), relative_residual=float(solved["relative_residual"]),
                gram_route=stats["gram_route"], field_route=field["route"], setup_s=setup_s, solve_s=solve_s,
                transfer_ledger=stats["assembly"]["transfer_ledger"], field_ledger=field["transfer_ledger"])


def _device_nonlinear(path, checked, probes, args, **images):
    import ngsolve as ng
    from radia.hdiv_gpu import PreparedHDivNonlinear
    from validation_ffag_connected_cyclic_yoke import _mean_tangential_magnetization, apply_cyclic_q2_geometry

    started = time.perf_counter()
    with PreparedHDivNonlinear.from_bh_table(path, np.array(BH_TABLE), backend="gpu", order=args.hdiv_order, curve_order=2,
                                             q2_deformation=lambda mesh: apply_cyclic_q2_geometry(mesh, args), **checked, **images) as device:
        setup_s = time.perf_counter() - started
        started = time.perf_counter()
        solved = device.solve(_applied(args), tolerance=args.nonlinear_tol, linear_tolerance=args.device_tol)
        solve_s = time.perf_counter() - started
        field = device.field(solved["x"], probes)
        stats = device.stats()
        magnetization = ng.GridFunction(device.fespace)
        magnetization.vec.FV().NumPy()[:] = solved["x"]
        mean = _mean_tangential_magnetization(dict(gfM=magnetization), device.fespace.mesh)
    return dict(mean_tangential=mean, field=np.asarray(field["field"]), dofs=int(device.dofs), newton_iterations=int(solved["newton_iterations"]),
                relative_residual=float(solved["relative_residual"]), field_route=field["route"], setup_s=setup_s, solve_s=solve_s,
                transfer_ledger=stats["assembly"]["transfer_ledger"], field_ledger=field["transfer_ledger"])


def _cpu_ring_equation(path, args):
    """The ring's CPU charge map and Gram with the mass on the device's material rule: the equation the device solves."""
    import ngsolve as ng
    from scipy.sparse import coo_matrix, csr_matrix
    from radia.vim._vim import _build_charge_gram_hex
    from validation_ffag_connected_cyclic_yoke import _mean_tangential_magnetization, apply_cyclic_q2_geometry

    mesh = ng.Mesh(str(path))
    keep = apply_cyclic_q2_geometry(mesh, args)
    space = ng.HDiv(mesh, order=args.hdiv_order)
    with ng.TaskManager():
        charge_map, gram, _, _ = _build_charge_gram_hex(space, eps=args.gram_eps, build_hmatrix=False, materialize_mass=False)
        u, v = space.TnT()
        mass = ng.BilinearForm(space)
        mass += u*v*ng.dx(intrules={ng.HEX: ng.IntegrationRule(ng.HEX, 2*max(6, 2*args.hdiv_order + 2))})
        mass.Assemble()
        source = ng.LinearForm(space)
        source += _applied(args)*v*ng.dx
        source.Assemble()
    gram.build_exact_dense_normalized_gram(maximum_memory_mb=8192)
    rows, cols, values = mass.mat.COO()
    mass_matrix = coo_matrix((values, (rows, cols)), shape=(space.ndof, space.ndof)).toarray()
    charge_map = csr_matrix(charge_map)
    energy = np.column_stack([charge_map.T @ np.asarray(gram.matvec_sym(np.ascontiguousarray(charge_map[:, k].toarray().ravel())))
                              for k in range(space.ndof)])
    x = np.linalg.solve(energy + mass_matrix/(args.mu_r - 1.0), np.array(source.vec.FV().NumPy()))
    magnetization = ng.GridFunction(space)
    magnetization.vec.FV().NumPy()[:] = x
    mean = _mean_tangential_magnetization(dict(gfM=magnetization), mesh)
    del keep
    return x, mean


def _public(run):
    return {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in run.items() if key != "x"}


def _nonlinear_case(paths, checked, probes, args):
    """The saturating yoke: the device sector against the device ring and against the CPU ring, on the compressed route."""
    import ngsolve as ng
    import radia as rad
    import radia.vim as vim
    from validation_ffag_connected_cyclic_yoke import _mean_tangential_magnetization, apply_cyclic_q2_geometry

    ring = ng.Mesh(str(paths["ring"]))
    keep = apply_cyclic_q2_geometry(ring, args)
    with ng.TaskManager():
        rad.UtiDelAll()
        cpu_ring = vim.Solve(ring, bh_table=np.array(BH_TABLE), H_ext=_applied(args), order=args.hdiv_order, gram_eps=args.gram_eps, leaf=args.leaf,
                             curve_order=2, tol=args.solve_tol, nl_tol=args.nonlinear_tol)
        cpu_mean = _mean_tangential_magnetization(cpu_ring, ring)
        cpu_h = np.asarray(vim.FieldFromSolution(cpu_ring, probes, algorithm="direct"), float)
    del keep
    drive = args.applied_field_am
    sector = _device_nonlinear(paths["sector"], checked["sector"], probes, args, image_cyclic=args.fold, cyclic_periodic_boundaries=SEAM)
    device_ring = _device_nonlinear(paths["ring"], checked["ring"], probes, args)
    measured = dict(
        cpu_ring_mean_tangential=cpu_mean, sector_mean_tangential=sector["mean_tangential"], device_ring_mean_tangential=device_ring["mean_tangential"],
        sector_vs_device_ring_mean=abs(sector["mean_tangential"] - device_ring["mean_tangential"])/abs(device_ring["mean_tangential"]),
        sector_vs_device_ring_field_relative_to_drive=float(np.max(np.linalg.norm(sector["field"] - device_ring["field"], axis=1))/drive),
        sector_vs_cpu_ring_mean=abs(sector["mean_tangential"] - cpu_mean)/abs(cpu_mean),
        device_ring_vs_cpu_ring_mean=abs(device_ring["mean_tangential"] - cpu_mean)/abs(cpu_mean),
        sector_vs_cpu_ring_field_relative_to_drive=float(np.max(np.linalg.norm(sector["field"] - cpu_h, axis=1))/drive))
    gates = dict(
        nonlinear_sector_mean_matches_device_ring=bool(measured["sector_vs_device_ring_mean"] <= PROPOSED_MAXIMUM_MEAN_RELATIVE_ERROR),
        nonlinear_sector_field_matches_device_ring=bool(measured["sector_vs_device_ring_field_relative_to_drive"] <= PROPOSED_MAXIMUM_FIELD_RELATIVE_TO_DRIVE),
        nonlinear_sector_mean_matches_cpu_ring=bool(measured["sector_vs_cpu_ring_mean"] <= PROPOSED_MAXIMUM_MEAN_RELATIVE_ERROR),
        nonlinear_sector_field_matches_cpu_ring=bool(measured["sector_vs_cpu_ring_field_relative_to_drive"] <= PROPOSED_MAXIMUM_FIELD_RELATIVE_TO_DRIVE))
    return dict(measured=measured, gates=gates, sector=_public(sector), device_ring=_public(device_ring))


def _case(cells, args, contracts):
    import ngsolve as ng
    import radia as rad
    import radia.vim as vim
    from device_scale import checked_inputs
    from validation_ffag_connected_cyclic_yoke import _mean_tangential_magnetization, apply_cyclic_q2_geometry

    geometry = dict(inner_radius=args.inner_radius, outer_radius=args.outer_radius, half_height=args.half_height)
    root = args.work/f"cells{cells}"
    root.mkdir(parents=True, exist_ok=False)
    paths = dict(ring=root/"ring.vol", sector=root/"sector.vol")
    build_yoke_mesh(fold=args.fold, cells=cells, full_ring=True, **geometry).ngmesh.Save(str(paths["ring"]))
    build_yoke_mesh(fold=args.fold, cells=cells, full_ring=False, **geometry).ngmesh.Save(str(paths["sector"]))
    checked = {name: checked_inputs(paths[name], contracts[name], args.receipt_python) for name in paths}
    probe_radius = np.linspace(0.8*args.inner_radius, 1.2*args.outer_radius, args.probe_count)
    probes = np.column_stack((probe_radius, np.zeros(args.probe_count), np.full(args.probe_count, 1.5*args.half_height)))

    ring, sector = ng.Mesh(str(paths["ring"])), ng.Mesh(str(paths["sector"]))
    keep = apply_cyclic_q2_geometry(ring, args), apply_cyclic_q2_geometry(sector, args)
    with ng.TaskManager():
        rad.UtiDelAll()
        cpu_ring = vim.Solve(ring, mu_r=args.mu_r, H_ext=_applied(args), order=args.hdiv_order, gram_eps=args.gram_eps,
                             leaf=args.leaf, curve_order=2, tol=args.solve_tol)
        rad.UtiDelAll()
        cpu_sector = vim.Solve(sector, mu_r=args.mu_r, H_ext=_applied(args), order=args.hdiv_order, image_cyclic=args.fold,
                               cyclic_periodic_boundaries=SEAM, curve_order=2, gram_eps=args.gram_eps, leaf=args.leaf,
                               tol=args.solve_tol)
        cpu_ring_x = np.array(cpu_ring["_m_coefficients"], dtype=float)
        cpu_ring_mean = _mean_tangential_magnetization(cpu_ring, ring)
        cpu_sector_mean = _mean_tangential_magnetization(cpu_sector, sector)
        cpu_ring_h = np.asarray(vim.FieldFromSolution(cpu_ring, probes, algorithm="direct"), float)
        cpu_sector_h = np.asarray(vim.FieldFromSolution(cpu_sector, probes, algorithm="direct"), float)
        cpu_sector_x = np.asarray(cpu_sector["_m_coefficients"], dtype=float)
    del keep
    equation_x, equation_mean = _cpu_ring_equation(paths["ring"], args)
    drive = args.applied_field_am
    runs, measured, gates = {}, {}, {}
    for compression in args.compressions:
        runs[compression] = dict(
            lone=_device_solve(paths["sector"], checked["sector"], probes, args, compression),
            sector=_device_solve(paths["sector"], checked["sector"], probes, args, compression, image_cyclic=args.fold,
                                 cyclic_periodic_boundaries=SEAM),
            ring=_device_solve(paths["ring"], checked["ring"], probes, args, compression, reference_state=cpu_ring_x))
        run = runs[compression]
        mean_error = abs(run["sector"]["mean_tangential"] - cpu_ring_mean)/abs(cpu_ring_mean)
        lone_error = abs(run["lone"]["mean_tangential"] - cpu_ring_mean)/abs(cpu_ring_mean)
        field_error = float(np.max(np.linalg.norm(run["sector"]["field"] - cpu_ring_h, axis=1))/drive)
        measured[compression] = dict(
            sector_mean_tangential=run["sector"]["mean_tangential"], ring_mean_tangential=run["ring"]["mean_tangential"],
            lone_mean_tangential=run["lone"]["mean_tangential"], sector_vs_cpu_ring_mean=mean_error, lone_vs_cpu_ring_mean=lone_error,
            sector_vs_cpu_ring_field_relative_to_drive=field_error,
            device_ring_vs_cpu_ring_mean=abs(run["ring"]["mean_tangential"] - cpu_ring_mean)/abs(cpu_ring_mean),
            device_ring_vs_cpu_ring_field_relative_to_drive=float(np.max(np.linalg.norm(run["ring"]["field"] - cpu_ring_h, axis=1))/drive),
            sector_x_vs_cpu_sector=float(np.linalg.norm(run["sector"]["x"] - cpu_sector_x)/np.linalg.norm(cpu_sector_x)),
            sector_vs_device_ring_mean=abs(run["sector"]["mean_tangential"] - run["ring"]["mean_tangential"])/abs(run["ring"]["mean_tangential"]),
            sector_vs_device_ring_field_relative_to_drive=float(np.max(np.linalg.norm(run["sector"]["field"] - run["ring"]["field"], axis=1))/drive),
            device_ring_x_vs_cpu_equation=float(np.linalg.norm(run["ring"]["x"] - equation_x)/np.linalg.norm(equation_x)),
            device_ring_vs_cpu_equation_mean=abs(run["ring"]["mean_tangential"] - equation_mean)/abs(equation_mean),
            cpu_equation_vs_cpu_ring_mean=abs(equation_mean - cpu_ring_mean)/abs(cpu_ring_mean),
            device_field_of_cpu_ring_state_vs_cpu_ring_field_relative_to_drive=float(
                np.max(np.linalg.norm(run["ring"]["field_of_reference_state"] - cpu_ring_h, axis=1))/drive))
        gates.update({
            f"{compression}_sector_dofs_scale_as_one_fold": bool(run["sector"]["dofs"]*args.fold == run["ring"]["dofs"]),
            f"{compression}_periodic_trace_is_compressed": bool(run["sector"]["dofs"] < run["lone"]["dofs"]),
            f"{compression}_sector_without_images_is_discriminating": bool(lone_error > 1.0e3*PROPOSED_MAXIMUM_MEAN_RELATIVE_ERROR),
            f"{compression}_sector_mean_matches_device_ring": bool(measured[compression]["sector_vs_device_ring_mean"] <= PROPOSED_MAXIMUM_MEAN_RELATIVE_ERROR),
            f"{compression}_sector_field_matches_device_ring": bool(
                measured[compression]["sector_vs_device_ring_field_relative_to_drive"] <= PROPOSED_MAXIMUM_FIELD_RELATIVE_TO_DRIVE),
            f"{compression}_device_ring_state_matches_cpu_equation": bool(measured[compression]["device_ring_x_vs_cpu_equation"] <= PROPOSED_MAXIMUM_STATE_RELATIVE_ERROR),
            f"{compression}_device_ring_mean_matches_cpu_equation": bool(measured[compression]["device_ring_vs_cpu_equation_mean"] <= PROPOSED_MAXIMUM_MEAN_RELATIVE_ERROR),
            f"{compression}_device_field_of_cpu_ring_state_matches_cpu_field": bool(
                measured[compression]["device_field_of_cpu_ring_state_vs_cpu_ring_field_relative_to_drive"] <= PROPOSED_MAXIMUM_FIELD_RELATIVE_TO_DRIVE),
            f"{compression}_images_crossed_only_as_rows": bool(
                next(row["bytes"] for row in run["sector"]["transfer_ledger"]
                     if row["payload"] == "image_transforms" and row["direction"] == "H2D") == 72*(args.fold - 1))})
    nonlinear = _nonlinear_case(paths, checked, probes, args)
    gates.update(nonlinear.pop("gates"))
    return dict(
        cells_per_sector=cells,
        nonlinear=nonlinear,
        cpu=dict(ring_mean_tangential=cpu_ring_mean, sector_mean_tangential=cpu_sector_mean, ring_dofs=int(cpu_ring["ndof"]),
                 sector_dofs=int(cpu_sector["ndof"]), periodic_slave_dofs=int(cpu_sector["periodic_slave_dofs"]),
                 sector_vs_ring_mean=abs(cpu_sector_mean - cpu_ring_mean)/abs(cpu_ring_mean),
                 sector_vs_ring_field_relative_to_drive=float(np.max(np.linalg.norm(cpu_sector_h - cpu_ring_h, axis=1))/drive),
                 ring_field_am=cpu_ring_h.tolist(), probe_points_m=probes.tolist()),
        measured=measured, gates=gates,
        runs={compression: {name: _public(run) for name, run in cases.items()} for compression, cases in runs.items()})


def run_validation(args):
    import ngsolve as ng
    from radia.ffag_topopt import validate_ffag_cyclic_sector_contract

    contract = validate_ffag_cyclic_sector_contract(args.fold, body_crosses_periodic_planes=True, periodic_trace_identified=True)
    args.work.mkdir(parents=True, exist_ok=False)
    contracts = dict(ring=args.work/"labels-ring.json", sector=args.work/"labels-sector.json")
    contracts["ring"].write_text(json.dumps(_contract(["skin"]), indent=2) + "\n", encoding="utf-8")
    contracts["sector"].write_text(json.dumps(_contract(["skin", *SEAM]), indent=2) + "\n", encoding="utf-8")
    started = time.perf_counter()
    ng.SetNumThreads(args.threads)
    cases = [_case(cells, args, contracts) for cells in args.cells_per_sector]
    gates = {f"cells{case['cells_per_sector']}_{name}": value for case in cases for name, value in case["gates"].items()}
    gates["connected_periodic_contract"] = bool(contract.reduction_mode == "connected-periodic-fem-sector")
    return dict(
        schema="radia.ffag-device-connected-yoke/v1", status="pass" if all(gates.values()) else "fail", adoption="proposed",
        accepted=False, scope="Connected curved-Q2 annular HEX yoke on the device route; the CPU equation of the explicit ring is the reference.",
        proposed_limits=dict(maximum_mean_relative_error=PROPOSED_MAXIMUM_MEAN_RELATIVE_ERROR,
                             maximum_state_relative_error=PROPOSED_MAXIMUM_STATE_RELATIVE_ERROR,
                             maximum_field_relative_to_drive=PROPOSED_MAXIMUM_FIELD_RELATIVE_TO_DRIVE),
        settings={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        source_sha=os.environ.get("RADIA_SOURCE_SHA"), ngsolve_version=ng.__version__, gates=gates, cases=cases,
        validation_wall_s=time.perf_counter() - started, runtime=dict(platform_class=platform.system(), python=platform.python_version()))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True, help="New folder for the meshes and their receipts")
    parser.add_argument("--receipt-python", type=Path, required=True, help="Interpreter of the mesh checker")
    parser.add_argument("--compressions", nargs="+", choices=("dense", "aca"), default=["dense", "aca"])
    parser.add_argument("--cells-per-sector", nargs="+", type=int, default=[1, 2])
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
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--leaf", type=int, default=32)
    parser.add_argument("--solve-tol", type=float, default=1.0e-10)
    parser.add_argument("--device-tol", type=float, default=1.0e-10)
    parser.add_argument("--nonlinear-tol", type=float, default=1.0e-9)
    parser.add_argument("--max-iterations", type=int, default=4000)
    parser.add_argument("--probe-count", type=int, default=5)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    report = run_validation(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(status=report["status"], gates=report["gates"],
                          measured={case["cells_per_sector"]: dict(cpu_sector_vs_ring=case["cpu"]["sector_vs_ring_mean"], device=case["measured"])
                                    for case in report["cases"]}), indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
