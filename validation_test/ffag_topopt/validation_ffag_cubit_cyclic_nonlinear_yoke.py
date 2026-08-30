"""End-to-end Cubit curved-Q2 FFAG sector validation with nonlinear BH.

This gate uses the meshes produced by ``build_ffag_cyclic_yoke_hex.py``.
It compares an explicit Cubit full ring with one Cubit sector carrying
Netgen periodic identifications, a conforming BDM2 trace, and cyclic nonlocal
images.  The constitutive law is a single-valued monotone BH table; hysteresis
is deliberately outside this validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _mean_tangential_magnetization(result, mesh):
    import ngsolve as ng

    radius = ng.sqrt(ng.x*ng.x + ng.y*ng.y)
    tangent = ng.CF((-ng.y/radius, ng.x/radius, 0.0))
    volume = float(ng.Integrate(ng.CF(1.0), mesh))
    return float(ng.Integrate(result["gfM"]*tangent, mesh))/volume


def _field_comparison_metrics(full_field, sector_field, full_mean, sector_mean):
    full_field = np.asarray(full_field, dtype=float)
    sector_field = np.asarray(sector_field, dtype=float)
    if full_field.shape != sector_field.shape or full_field.ndim != 2:
        raise ValueError("full and sector field arrays must have equal shape")
    field_absolute = float(np.max(np.linalg.norm(
        sector_field-full_field, axis=1)))
    residual_field_scale = max(
        float(np.max(np.linalg.norm(full_field, axis=1))),
        np.finfo(float).tiny)
    # The demagnetizing field above a closed ring is a small residual of much
    # larger magnetization sources.  Keep the cancellation-sensitive relative
    # number as a diagnostic, but gate the integral against its source scale.
    field_source_scale = max(abs(float(full_mean)), abs(float(sector_mean)),
                             residual_field_scale)
    return {
        "field_maximum_absolute_error_am": field_absolute,
        "field_residual_relative_error": (
            field_absolute/residual_field_scale),
        "field_source_scale_am": field_source_scale,
        "field_source_relative_error": field_absolute/field_source_scale,
    }


def _run_builder(args):
    script = Path(__file__).with_name("build_ffag_cyclic_yoke_hex.py")
    command = [
        sys.executable, str(script),
        "--sector-output", str(args.sector_mesh.resolve()),
        "--full-output", str(args.full_mesh.resolve()),
        "--report", str(args.build_report.resolve()),
        "--fold", str(args.fold),
        "--inner-radius", str(args.inner_radius),
        "--outer-radius", str(args.outer_radius),
        "--half-height", str(args.half_height),
        "--intervals", str(args.intervals),
    ]
    completed = subprocess.run(command, check=False)
    if completed.returncode:
        raise RuntimeError(
            "canonical Cubit FFAG mesh build failed with exit code %d"
            % completed.returncode)


def run_validation(args):
    import ngsolve as ng
    from cubit_mesh_export.check import check_consistency

    import radia as rad
    from radia import vim
    from radia.ffag_topopt import (
        build_ffag_cyclic_density_map,
        validate_ffag_cyclic_sector_contract,
    )

    if args.rebuild_cubit or not (
            args.sector_mesh.is_file() and args.full_mesh.is_file()
            and args.build_report.is_file()):
        _run_builder(args)
    build_report = json.loads(args.build_report.read_text(encoding="utf-8"))
    if not build_report.get("passed", False):
        raise RuntimeError("canonical Cubit FFAG build report is not passing")
    sector_check = check_consistency(
        args.sector_mesh, min_curve_order=2,
        required_materials=("yoke",),
        required_boundaries=("skin", "periodic_min", "periodic_max"))
    full_check = check_consistency(
        args.full_mesh, min_curve_order=2,
        required_materials=("yoke",), required_boundaries=("skin",))
    if not sector_check["passed"] or not full_check["passed"]:
        raise RuntimeError("check-vol rejected the nonlinear FFAG mesh pair")

    sector = ng.Mesh(str(args.sector_mesh))
    full = ng.Mesh(str(args.full_mesh))
    contract = validate_ffag_cyclic_sector_contract(
        args.fold, body_crosses_periodic_planes=True,
        periodic_trace_identified=True)
    density_map = build_ffag_cyclic_density_map(sector)
    reduced_trial = np.linspace(0.2, 0.8, density_map.variable_count)
    element_trial = density_map.expand(reduced_trial)
    roundtrip = density_map.reduce(element_trial)

    mu0 = 4.0e-7*math.pi
    bh_table = np.asarray([
        [0.0, 0.0],
        [100.0, 0.25],
        [300.0, 0.75],
        [1.0e3, 1.30],
        [3.0e3, 1.55],
        [1.0e4, 1.75],
        [1.0e5, 2.00],
    ], dtype=float)
    radius = ng.sqrt(ng.x*ng.x + ng.y*ng.y)
    applied = ng.CF((
        -args.applied_field_am*ng.y/radius,
        args.applied_field_am*ng.x/radius, 0.0))
    options = dict(
        bh_table=bh_table, H_ext=applied, order=2, curve_order=2,
        gram_eps=args.gram_eps, leaf=args.leaf, tol=args.solve_tol,
        maxit=args.linear_maxit, nl_tol=args.nonlinear_tol,
        nl_maxit=args.nonlinear_maxit,
        nonlinear_solver="picard-mass-riesz",
        preconditioner="mass-riesz")
    ng.SetNumThreads(args.threads)
    started = time.perf_counter()
    with ng.TaskManager():
        rad.UtiDelAll()
        full_result = vim.Solve(full, **options)
        rad.UtiDelAll()
        sector_result = vim.Solve(
            sector, image_cyclic=args.fold,
            cyclic_periodic_boundaries=("periodic_min", "periodic_max"),
            **options)
        full_mean = _mean_tangential_magnetization(full_result, full)
        sector_mean = _mean_tangential_magnetization(sector_result, sector)
        angle = math.pi/args.fold
        radii = np.linspace(
            0.8*args.inner_radius, 1.2*args.outer_radius,
            args.probe_count)
        probes = np.column_stack((
            radii*np.cos(angle), radii*np.sin(angle),
            np.full(args.probe_count, 1.5*args.half_height)))
        full_field = np.asarray(vim.FieldFromSolution(
            full_result, probes, algorithm="direct"), dtype=float)
        sector_field = np.asarray(vim.FieldFromSolution(
            sector_result, probes, algorithm="direct"), dtype=float)
    wall = time.perf_counter()-started

    mean_scale = max(abs(full_mean), np.finfo(float).tiny)
    mean_relative = abs(sector_mean-full_mean)/mean_scale
    field_metrics = _field_comparison_metrics(
        full_field, sector_field, full_mean, sector_mean)
    checks = {
        "cubit_build_report_passed": bool(build_report["passed"]),
        "cubit_vol_checks_passed": bool(
            sector_check["passed"] and full_check["passed"]),
        "connected_periodic_contract": (
            contract.reduction_mode == "connected-periodic-fem-sector"),
        "periodic_material_density_roundtrip": bool(np.array_equal(
            roundtrip, reduced_trial)),
        "periodic_material_variables_are_reduced": bool(
            density_map.variable_count < density_map.element_count),
        "full_ring_nonlinear_cpp_picard_converged": bool(
            full_result["last_solve_converged"]
            and full_result["linear_solver"] == "picard-mass-riesz-cpp"
            and full_result["last_solve_final_relative_residual"]
            <= args.maximum_nonlinear_residual),
        "sector_nonlinear_cpp_picard_converged": bool(
            sector_result["last_solve_converged"]
            and sector_result["linear_solver"] == "picard-mass-riesz-cpp"
            and sector_result["last_solve_final_relative_residual"]
            <= args.maximum_nonlinear_residual),
        "nonlinear_iteration_is_exercised": bool(
            int(full_result["iters"]) >= 2 and int(sector_result["iters"]) >= 2),
        "mean_magnetization_matches_full_ring": bool(
            mean_relative <= args.maximum_relative_error),
        "field_matches_full_ring": bool(
            field_metrics["field_source_relative_error"]
            <= args.maximum_relative_error),
    }
    result = {
        "schema": "radia.ffag-cubit-cyclic-nonlinear-yoke/v1",
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "machine": platform.node(),
        "platform": platform.platform(),
        "model": {
            "geometry_source": "Coreform Cubit 2025.12 ACIS sweep-map",
            "fold": int(args.fold),
            "hdiv_family": "BDM2",
            "curve_order": 2,
            "constitutive_model": "single-valued monotone BH",
            "nonlinear_solver": "C++ Picard + mass-Riesz",
            "hysteresis": False,
            "mu0_henry_per_m": mu0,
            "bh_table_h_am_b_t": bh_table.tolist(),
            "applied_tangential_h_am": float(args.applied_field_am),
        },
        "meshes": {
            "sector": str(args.sector_mesh.resolve()),
            "sector_sha256": _sha256(args.sector_mesh),
            "sector_elements": int(sector.ne),
            "sector_dofs": int(sector_result["ndof"]),
            "full_ring": str(args.full_mesh.resolve()),
            "full_ring_sha256": _sha256(args.full_mesh),
            "full_ring_elements": int(full.ne),
            "full_ring_dofs": int(full_result["ndof"]),
        },
        "periodic_topology": {
            "boundary_facet_pairs": int(density_map.boundary_pair_count),
            "element_count": int(density_map.element_count),
            "independent_density_variable_count": int(
                density_map.variable_count),
            "element_to_variable": density_map.element_to_variable.tolist(),
        },
        "nonlinear_solve": {
            "full_ring_iterations": int(full_result["iters"]),
            "sector_iterations": int(sector_result["iters"]),
            "full_ring_final_relative_residual": float(
                full_result["last_solve_final_relative_residual"]),
            "sector_final_relative_residual": float(
                sector_result["last_solve_final_relative_residual"]),
            "full_ring_timing": dict(full_result.get("timing", {})),
            "sector_timing": dict(sector_result.get("timing", {})),
        },
        "comparison": {
            "full_ring_mean_tangential_m_am": full_mean,
            "sector_mean_tangential_m_am": sector_mean,
            "mean_magnetization_relative_error": mean_relative,
            "probe_points_m": probes.tolist(),
            "full_ring_demag_h_am": full_field.tolist(),
            "sector_demag_h_am": sector_field.tolist(),
            **field_metrics,
            "maximum_allowed_relative_error": float(
                args.maximum_relative_error),
        },
        "wall_s": wall,
        "checks": checks,
        "passed": all(checks.values()),
    }
    args.output.resolve().write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise RuntimeError("nonlinear Cubit FFAG cyclic validation failed")
    return result


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    root = Path(r"C:\temp\radia_ffag_cyclic")
    parser.add_argument("--sector-mesh", type=Path,
                        default=root/"ffag_cyclic_sector_q2.vol")
    parser.add_argument("--full-mesh", type=Path,
                        default=root/"ffag_cyclic_full_q2.vol")
    parser.add_argument("--build-report", type=Path,
                        default=root/"ffag_cyclic_mesh_build.json")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name(
        "results_ffag_cubit_cyclic_nonlinear_yoke.json"))
    parser.add_argument("--rebuild-cubit", action="store_true")
    parser.add_argument("--fold", type=int, default=12)
    parser.add_argument("--inner-radius", type=float, default=1.0)
    parser.add_argument("--outer-radius", type=float, default=2.0)
    parser.add_argument("--half-height", type=float, default=0.25)
    parser.add_argument("--intervals", type=int, default=2)
    parser.add_argument("--applied-field-am", type=float, default=300.0)
    parser.add_argument("--probe-count", type=int, default=5)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--gram-eps", type=float, default=1.0e-13)
    parser.add_argument("--leaf", type=int, default=4096)
    parser.add_argument("--solve-tol", type=float, default=1.0e-10)
    parser.add_argument("--linear-maxit", type=int, default=4000)
    parser.add_argument("--nonlinear-tol", type=float, default=1.0e-9)
    parser.add_argument("--nonlinear-maxit", type=int, default=100)
    parser.add_argument("--maximum-nonlinear-residual", type=float,
                        default=2.0e-9)
    parser.add_argument("--maximum-relative-error", type=float, default=2.0e-8)
    return parser.parse_args(argv)


def main():
    args = parse_args()
    if (args.fold < 3 or args.intervals < 2 or args.probe_count < 2
            or args.maximum_relative_error <= 0.0
            or args.maximum_nonlinear_residual <= 0.0):
        raise ValueError("invalid nonlinear FFAG cyclic validation settings")
    args.output.resolve().parent.mkdir(parents=True, exist_ok=True)
    run_validation(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
