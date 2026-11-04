#!/usr/bin/env python
"""Device ``image_cyclic=N`` cell against the explicitly meshed N-cell ring: disjoint poles, JSON record.

The contract of ``validation_ffag_cyclic_bdm2_ring.py`` on the device route. One two-pole HEX cell solved on the
device with its rotated images must reproduce the ring of N exact rotations of that cell, in the mean magnetization and
in gap-field samples. The ring is solved on the CPU route, the validated reference, and on the device route; the
isolated cell shows that a no-op image route could not pass. Every operator of a device solve is built and kept on the
device, and the transfer ledgers are recorded. The limits repeat the CPU route's numbers but are PROPOSED for the device
route, not adopted.
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

PROPOSED_MINIMUM_GAP_FRACTION_CLOSED = 0.999
PROPOSED_MAXIMUM_FIELD_RELATIVE_ERROR = 2.0e-8
_CONTRACT = dict(schema="radia.vol-label-contract.v1", application="ffag_cyclic_pole_ring",
                 metadata=dict(length_unit="m", unit_source="Caller declaration; the VOL file has no unit tag"),
                 strict_labels=True, required=dict(materials=["pole"], boundaries=["pole_boundary"]),
                 allowed=dict(materials=["pole"], boundaries=["pole_boundary"]), quality=dict(min_curve_order=1))


def _mean_mz(result, mesh):
    import ngsolve as ng
    return float(ng.Integrate(result["gfM"][2], mesh))/float(ng.Integrate(ng.CF(1.0), mesh))


def _device_solve(path, checked, probes, args, compression, **images):
    from radia.hdiv_gpu import PreparedHDivLinear
    started = time.perf_counter()
    with PreparedHDivLinear(path, backend="gpu", order=args.hdiv_order, gram_compression=compression,
                            **checked, **images) as device:
        setup_s = time.perf_counter() - started
        started = time.perf_counter()
        solved = device.solve((0.0, 0.0, args.applied_field_am), 1.0/(args.mu_r - 1.0), tolerance=args.device_tol,
                              max_iterations=args.max_iterations)
        solve_s = time.perf_counter() - started
        field = device.field(solved["x"], probes)
        stats = device.stats()
    return dict(x=solved["x"], mean_mz=float(solved["mean_magnetization"][2]), field=np.asarray(field["field"]),
                dofs=int(device.dofs), iterations=int(solved["iterations"]), relative_residual=float(solved["relative_residual"]),
                gram_route=stats["gram_route"], field_route=field["route"], setup_s=setup_s, solve_s=solve_s,
                transfer_ledger=stats["assembly"]["transfer_ledger"], field_ledger=field["transfer_ledger"])


def _public(run):
    return {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in run.items() if key != "x"}


def run_validation(args):
    import ngsolve as ng
    import radia as rad
    import radia.vim as vim
    from device_scale import checked_inputs
    from radia.ffag_topopt import validate_ffag_cyclic_sector_contract
    from validation_ffag_cyclic_bdm2_ring import build_ffag_pole_ring_mesh

    validate_ffag_cyclic_sector_contract(args.fold)
    common = dict(fold=args.fold, radius=args.radius, radial_half_width=args.radial_half_width,
                  tangential_half_width=args.tangential_half_width, half_gap=args.half_gap,
                  pole_thickness=args.pole_thickness)
    args.work.mkdir(parents=True, exist_ok=False)
    contract = args.work/"labels.json"
    contract.write_text(json.dumps(_CONTRACT, indent=2) + "\n", encoding="utf-8")
    started = time.perf_counter()
    ng.SetNumThreads(args.threads)
    offsets = np.linspace(-0.75*args.radial_half_width, 0.75*args.radial_half_width, args.probe_count)
    probes = np.column_stack((args.radius + offsets, np.zeros(args.probe_count), np.zeros(args.probe_count)))
    paths, checked = {}, {}
    for name, sectors in (("ring", range(args.fold)), ("cell", (0,))):
        paths[name] = args.work/(name + ".vol")
        build_ffag_pole_ring_mesh(sectors, **common).ngmesh.Save(str(paths[name]))
        checked[name] = checked_inputs(paths[name], contract, args.receipt_python)
    with ng.TaskManager():
        ring, cell = ng.Mesh(str(paths["ring"])), ng.Mesh(str(paths["cell"]))
        external = ng.CF((0.0, 0.0, args.applied_field_am))
        rad.UtiDelAll()
        cpu_ring = vim.Solve(ring, mu_r=args.mu_r, H_ext=external, order=args.hdiv_order, gram_eps=args.gram_eps, tol=args.solve_tol)
        cpu_ring_mz = _mean_mz(cpu_ring, ring)
        cpu_ring_h = np.asarray(vim.FieldFromSolution(cpu_ring, probes, algorithm="direct"), float)
        rad.UtiDelAll()
        cpu_cell = vim.Solve(cell, mu_r=args.mu_r, H_ext=external, order=args.hdiv_order, image_cyclic=args.fold,
                             gram_eps=args.gram_eps, tol=args.solve_tol)
        cpu_cell_mz = _mean_mz(cpu_cell, cell)
        cpu_cell_x = np.asarray(cpu_cell["_m_coefficients"], dtype=float)
    runs = {}
    for compression in args.compressions:
        runs[compression] = dict(
            lone=_device_solve(paths["cell"], checked["cell"], probes, args, compression),
            cyclic=_device_solve(paths["cell"], checked["cell"], probes, args, compression, image_cyclic=args.fold),
            ring=_device_solve(paths["ring"], checked["ring"], probes, args, compression))

    field_scale = max(float(np.max(np.linalg.norm(cpu_ring_h, axis=1))), 1.0)
    gates, measured = {}, {}
    for compression, run in runs.items():
        gap = abs(run["lone"]["mean_mz"] - cpu_ring_mz)
        residual = abs(run["cyclic"]["mean_mz"] - cpu_ring_mz)
        closed = 1.0 - residual/gap if gap > 0.0 else -math.inf
        field_error = float(np.max(np.linalg.norm(run["cyclic"]["field"] - cpu_ring_h, axis=1))/field_scale)
        lone_error = float(np.max(np.linalg.norm(run["lone"]["field"] - cpu_ring_h, axis=1))/field_scale)
        measured[compression] = dict(
            lone_cell_mean_mz=run["lone"]["mean_mz"], cyclic_cell_mean_mz=run["cyclic"]["mean_mz"], ring_mean_mz=run["ring"]["mean_mz"],
            lone_vs_cpu_ring_gap=gap, cyclic_vs_cpu_ring_residual=residual, gap_fraction_closed=closed,
            cyclic_vs_cpu_ring_field=field_error, lone_vs_cpu_ring_field=lone_error,
            cyclic_vs_device_ring_field=float(np.max(np.linalg.norm(run["cyclic"]["field"] - run["ring"]["field"], axis=1))/field_scale),
            device_ring_vs_cpu_ring_field=float(np.max(np.linalg.norm(run["ring"]["field"] - cpu_ring_h, axis=1))/field_scale),
            cyclic_x_vs_cpu_cell=float(np.linalg.norm(run["cyclic"]["x"] - cpu_cell_x)/np.linalg.norm(cpu_cell_x)))
        gates.update({
            f"{compression}_lone_cell_is_discriminating": bool(gap/max(abs(cpu_ring_mz), 1.0) > 1.0e-3),
            f"{compression}_lone_cell_field_differs": bool(lone_error > 1.0e3*PROPOSED_MAXIMUM_FIELD_RELATIVE_ERROR),
            f"{compression}_cyclic_closes_magnetization_gap": bool(closed >= PROPOSED_MINIMUM_GAP_FRACTION_CLOSED),
            f"{compression}_cyclic_gap_field_matches_ring": bool(field_error <= PROPOSED_MAXIMUM_FIELD_RELATIVE_ERROR),
            f"{compression}_one_cell_reduces_dofs": bool(run["cyclic"]["dofs"]*(args.fold - 1) < run["ring"]["dofs"]),
            f"{compression}_images_crossed_only_as_rows": bool(
                next(row["bytes"] for row in run["cyclic"]["transfer_ledger"]
                     if row["payload"] == "image_transforms" and row["direction"] == "H2D") == 72*(args.fold - 1))})
    gates["bdm2_hex_is_used"] = bool(args.hdiv_order == 2 and {len(e.vertices) for e in cell.Elements(ng.VOL)} == {8})
    return dict(
        schema="radia.ffag-device-cyclic-ring/v1", status="pass" if all(gates.values()) else "fail", adoption="proposed",
        accepted=False, scope="Disjoint two-pole FFAG cells on the device route; the CPU ring is the reference.",
        proposed_limits=dict(minimum_gap_fraction_closed=PROPOSED_MINIMUM_GAP_FRACTION_CLOSED,
                             maximum_field_relative_error=PROPOSED_MAXIMUM_FIELD_RELATIVE_ERROR),
        settings={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        source_sha=os.environ.get("RADIA_SOURCE_SHA"), ngsolve_version=ng.__version__,
        cpu=dict(ring_mean_mz=cpu_ring_mz, cell_mean_mz=cpu_cell_mz, ring_dofs=int(cpu_ring["ndof"]), cell_dofs=int(cpu_cell["ndof"]),
                 ring_field_am=cpu_ring_h.tolist(), probe_points_m=probes.tolist()),
        measured=measured, gates=gates, runs={compression: {name: _public(run) for name, run in cases.items()} for compression, cases in runs.items()},
        validation_wall_s=time.perf_counter() - started, runtime=dict(platform_class=platform.system(), python=platform.python_version()))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True, help="New folder for the meshes and their receipts")
    parser.add_argument("--receipt-python", type=Path, required=True, help="Interpreter of the mesh checker")
    parser.add_argument("--compressions", nargs="+", choices=("dense", "aca"), default=["dense", "aca"])
    parser.add_argument("--fold", type=int, default=6)
    parser.add_argument("--hdiv-order", type=int, default=2, choices=(1, 2))
    parser.add_argument("--radius", type=float, default=0.040)
    parser.add_argument("--radial-half-width", type=float, default=0.008)
    parser.add_argument("--tangential-half-width", type=float, default=0.006)
    parser.add_argument("--half-gap", type=float, default=0.004)
    parser.add_argument("--pole-thickness", type=float, default=0.006)
    parser.add_argument("--mu-r", type=float, default=500.0)
    parser.add_argument("--applied-field-am", type=float, default=1.0e5)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--solve-tol", type=float, default=1.0e-10)
    parser.add_argument("--device-tol", type=float, default=1.0e-10)
    parser.add_argument("--max-iterations", type=int, default=4000)
    parser.add_argument("--probe-count", type=int, default=5)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.fold < 2 or args.probe_count < 2 or args.threads < 1 or args.mu_r <= 1.0:
        raise ValueError("invalid FFAG device cyclic settings")
    report = run_validation(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(status=report["status"], measured=report["measured"], gates=report["gates"]), indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
