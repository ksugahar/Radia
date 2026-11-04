"""Ideal sphere on connected tet, pyramid and hex meshes: mesh and order convergence of the device route; JSON record.

A uniform field magnetizes an ideal sphere uniformly: M = 3 chi/(3 + chi) H0 for a linear
material and H = H0 - M/3 on the BH curve for a nonlinear one. Straight cells converge to the
sphere as they refine, so the error also holds the geometry error, which the volume ratio
records. The CPU route runs on the coarse levels as a second implementation of one
discretization. The limit is proposed, not adopted. Timings are diagnostic.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import time

import numpy as np

from device_scale import carried_inputs, checked_inputs, sha256
from hdiv_mixed_meshes import mixed_grid, to_ball
from validate_device_nonlinear import bh_table

PROPOSED_MEAN_LIMIT = 0.01


def analytic_nonlinear(table, field):
    from radia.vim._nonlinear import _bh_inverse_funcs
    from scipy.optimize import brentq
    curves, _, _ = _bh_inverse_funcs(table[:, 0], table[:, 1])
    return brentq(lambda value: float(curves(np.array([value]))[0][0])*value + value/3 - field, 0, 3*field, xtol=1e-6)


def error(mean, expected):
    return float(np.linalg.norm(np.asarray(mean) - (0., 0., expected))/abs(expected))


def core_error(mesh, order, coefficients, expected):
    import ngsolve as ng
    grid = np.linspace(-.45, .45, 4)
    points = [(x + .013, y - .009, z + .017) for x in grid for y in grid for z in grid]
    field = ng.GridFunction(ng.HDiv(mesh, order=order))
    field.vec.FV().NumPy()[:] = coefficients
    values = np.array([tuple(field(mesh(*point))) for point in points])
    return float(np.sqrt(np.mean(np.sum((values - (0., 0., expected))**2, axis=1)))/abs(expected))


def device_level(path, order, args, checked, table, linear_expected, nonlinear_expected, mesh):
    from radia.hdiv_gpu import PreparedHDivLinear, PreparedHDivNonlinear
    start = time.perf_counter()
    with PreparedHDivLinear(path, backend="gpu", order=order, gram_compression="aca", **checked) as prepared:
        setup = time.perf_counter() - start
        start = time.perf_counter()
        linear = prepared.solve((0., 0., args.field), 1/args.chi, tolerance=1e-10, max_iterations=args.max_iterations,
                                preconditioner="mass-diagonal")
        solve = time.perf_counter() - start
        stats = prepared.stats()
    start = time.perf_counter()
    with PreparedHDivNonlinear.from_bh_table(path, table, backend="gpu", order=order, **checked) as prepared:
        nonlinear_setup = time.perf_counter() - start
        start = time.perf_counter()
        nonlinear = prepared.solve((0., 0., args.field), tolerance=1e-10, max_iterations=args.max_newton, linear_tolerance=args.linear_tolerance,
                                   max_linear_iterations=args.max_linear_iterations, preconditioner=args.preconditioner,
                                   degree=args.degree, lower_fraction=args.lower_fraction)
        nonlinear_solve = time.perf_counter() - start
    coefficients = dict(linear=linear["x"], nonlinear=nonlinear["x"])
    return coefficients, dict(dofs=stats["dofs"], family=stats["family"], hmatrix=stats["assembly"]["hmatrix"],
                setup_peak_device_bytes=stats["assembly"]["setup_peak_device_bytes"],
                operator_owned_device_bytes=stats["operator"]["owned_device_bytes"],
                setup_seconds_diagnostic=setup, solve_seconds_diagnostic=solve,
                nonlinear_setup_seconds_diagnostic=nonlinear_setup, nonlinear_solve_seconds_diagnostic=nonlinear_solve,
                linear=dict(mean_error=error(linear["mean_magnetization"], linear_expected),
                            core_error=core_error(mesh, order, linear["x"], linear_expected),
                            relative_residual=float(linear["relative_residual"]), iterations=int(linear["iterations"]),
                            mean=[float(v) for v in linear["mean_magnetization"]]),
                nonlinear=dict(mean_error=error(nonlinear["mean_magnetization"], nonlinear_expected),
                               core_error=core_error(mesh, order, nonlinear["x"], nonlinear_expected),
                               relative_residual=float(nonlinear["relative_residual"]),
                               newton_iterations=int(nonlinear["newton_iterations"]),
                               linear_iterations=[int(row["linear_iterations"]) for row in nonlinear["history"]],
                               mean=[float(v) for v in nonlinear["mean_magnetization"]]))


def cpu_level(mesh, order, args, table, device, coefficients):
    import ngsolve as ng
    from radia.vim import Solve
    applied = ng.CoefficientFunction((0, 0, args.field))
    with ng.TaskManager():
        start = time.perf_counter()
        linear = Solve(mesh, mu_r=1 + args.chi, H_ext=applied, order=order, tol=1e-10)
        linear_seconds = time.perf_counter() - start
        start = time.perf_counter()
        nonlinear = Solve(mesh, bh_table=table, H_ext=applied, order=order, nl_tol=1e-9)
        nonlinear_seconds = time.perf_counter() - start
    def relative(result, x):
        reference = result["gfM"].vec.FV().NumPy()
        return float(np.linalg.norm(reference - x)/np.linalg.norm(reference))
    return dict(linear_seconds_diagnostic=linear_seconds, nonlinear_seconds_diagnostic=nonlinear_seconds,
                linear_iterations=int(linear["iters"]), newton_iterations=int(nonlinear["iters"]),
                linear_coefficients_vs_device=relative(linear, coefficients["linear"]),
                nonlinear_coefficients_vs_device=relative(nonlinear, coefficients["nonlinear"]),
                linear_mean_vs_device=float(np.linalg.norm(np.asarray(linear["M_avg"]) - device["linear"]["mean"])/np.linalg.norm(linear["M_avg"])),
                nonlinear_mean_vs_device=float(np.linalg.norm(np.asarray(nonlinear["M_avg"]) - device["nonlinear"]["mean"])/np.linalg.norm(nonlinear["M_avg"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    meshes = parser.add_mutually_exclusive_group(required=True)
    meshes.add_argument("--receipt-python", type=Path, help="Generate each ball and check it with this interpreter")
    meshes.add_argument("--prepared", type=Path, help="Folder of ball-<n>.vol already checked, with their receipts")
    parser.add_argument("--levels", type=int, nargs="+", default=[2, 4, 6, 8])
    parser.add_argument("--orders", type=int, nargs="+", choices=[1, 2], default=[1, 2])
    parser.add_argument("--cpu-levels", type=int, nargs="*", default=[2, 4])
    parser.add_argument("--field", type=float, default=1e5)
    parser.add_argument("--chi", type=float, default=10.)
    parser.add_argument("--max-iterations", type=int, default=4000)
    parser.add_argument("--max-newton", type=int, default=100)
    parser.add_argument("--linear-tolerance", type=float, default=1e-6)
    parser.add_argument("--max-linear-iterations", type=int, default=30000)
    parser.add_argument("--preconditioner", default="reference-diagonal")
    parser.add_argument("--degree", type=int, default=32)
    parser.add_argument("--lower-fraction", type=float, default=.01)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    import ngsolve as ng
    table = bh_table()
    linear_expected = 3*args.chi/(3 + args.chi)*args.field
    nonlinear_expected = analytic_nonlinear(table, args.field)
    result = dict(schema="radia.hdiv-mixed-sphere.v1", accepted=False, units="SI", adoption="proposed",
                  proposed_mean_limit=PROPOSED_MEAN_LIMIT, field_A_per_m=args.field, chi=args.chi,
                  newton_controls=dict(max_newton=args.max_newton, linear_tolerance=args.linear_tolerance, max_linear_iterations=args.max_linear_iterations,
                                       preconditioner=args.preconditioner, degree=args.degree, lower_fraction=args.lower_fraction),
                  linear_expected_M=linear_expected, nonlinear_expected_M=nonlinear_expected,
                  ngsolve_version=ng.__version__, machine=platform.machine(), source_sha=os.environ.get("RADIA_SOURCE_SHA"), cases=[])
    try:
        for n in args.levels:
            if args.prepared:
                path = args.prepared/f"ball-{n}.vol"
                checked = carried_inputs(path, args.contract)
            else:
                path = args.out/f"ball-{n}.vol"
                mesh = mixed_grid((n, n, n), lengths=(2., 2., 2.), origin=(-1., -1., -1.), mapping=to_ball(1.))
                mesh.ngmesh.Save(str(path))
                checked = checked_inputs(path, args.contract, args.receipt_python)
            mesh = ng.Mesh(str(path))
            counts = {}
            for element in mesh.Elements(ng.VOL):
                counts[len(element.vertices)] = counts.get(len(element.vertices), 0) + 1
            volume_ratio = float(ng.Integrate(ng.CF(1), mesh, order=6)/(4*np.pi/3))
            for order in args.orders:
                coefficients, record = device_level(path, order, args, checked, table, linear_expected, nonlinear_expected, mesh)
                case = dict(n=n, order=order, cells=counts, volume_ratio=volume_ratio, mesh_sha256=sha256(path), **record)
                if n in args.cpu_levels:
                    case["cpu"] = cpu_level(mesh, order, args, table, case, coefficients)
                result["cases"].append(case)
                (args.out/"result.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
        finest = [case for case in result["cases"] if case["n"] == max(args.levels)]
        result["accepted"] = bool(finest) and all(
            case["linear"]["relative_residual"] <= 1e-10 and case["nonlinear"]["relative_residual"] <= 1e-10
            and case["linear"]["mean_error"] <= PROPOSED_MEAN_LIMIT and case["nonlinear"]["mean_error"] <= PROPOSED_MEAN_LIMIT
            for case in finest)
    except BaseException as exc:
        result["error"] = repr(exc)
        raise
    finally:
        (args.out/"result.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
