"""Compare HDiv-MMM and Omega-reduced-Omega on one C-type magnet.

HCurl reduced-A remains an independent third route unless ``--primary-only``
is selected. All formulations share the exact Cubit/ACIS iron authority, one
solid ``CoilBuilder`` excitation, one B-H table, and one set of physical
observation points. HDiv-MMM intentionally uses the iron-only mesh because its
Coulomb Gram is the open-boundary operator; the FEM formulations use the same
periodic Kelvin mesh. No finite outer air box is part of this comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import sys
import time
from pathlib import Path

import ngsolve as ng
import numpy as np
import radia as rad

from radia import vim
from radia import _radia_pybind as radia_native
from radia.coil_builder import CoilBuilder
from radia.kelvin_identify_ngsolve import (
    detect_kelvin_offset,
    has_kelvin_identification,
)
from radia.scalar_potential_solver import ScalarPotentialSolver
from radia.vector_potential_solver import VectorPotentialSolver


MU0 = 4.0e-7 * math.pi
HERE = Path(__file__).resolve().parent
DEFAULT_BH = (
    Path(rad.__file__).resolve().parent
    / "panels"
    / "samples"
    / "em_sample_bh.txt"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def progress(event: str, **values) -> None:
    print(json.dumps({"event": event, **values}, sort_keys=True), flush=True)


def write_checkpoint(path: Path, contract: dict[str, object], field: np.ndarray,
                     diagnostics: dict[str, object]) -> None:
    payload = {
        "schema": "radia.validation.c-type-three-engine-checkpoint.v1",
        "contract": contract,
        "field_T": np.asarray(field, dtype=float).tolist(),
        "diagnostics": diagnostics,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def read_checkpoint(path: Path, contract: dict[str, object]):
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (payload.get("schema") != "radia.validation.c-type-three-engine-checkpoint.v1"
            or payload.get("contract") != contract):
        raise RuntimeError(
            f"checkpoint contract differs from this run; remove or rename {path}"
        )
    return (
        np.asarray(payload["field_T"], dtype=float).reshape(-1, 3),
        dict(payload["diagnostics"]),
    )


def build_coil() -> tuple[int, dict[str, object]]:
    radius = 0.0225
    straight_x = 0.050
    straight_y = 0.0625
    centre = np.array([0.0, 0.13125, 0.0])
    start = centre + np.array([0.5 * straight_x + radius, -0.5 * straight_y, 0.0])
    builder = (
        CoilBuilder(-2000.0)
        .set_start(start)
        .set_cross_section(0.035, 0.105)
        .add_straight(straight_y)
        .add_arc(radius, 90.0)
        .add_straight(straight_x)
        .add_arc(radius, 90.0)
        .add_straight(straight_y)
        .add_arc(radius, 90.0)
        .add_straight(straight_x)
        .add_arc(radius, 90.0)
    )
    if not builder.is_closed or builder.gap > 1e-12:
        raise RuntimeError(f"C-type CoilBuilder path is open by {builder.gap:.6e} m")
    objects = builder.to_radia(arc_max_segment_length=0.004)
    if not objects:
        raise RuntimeError("C-type CoilBuilder produced no Radia current objects")
    return rad.ObjCnt(objects), {
        "authority": "ESRF Example 5 rounded rectangular solid current",
        "current_A": -2000.0,
        "radius_m": radius,
        "straight_x_m": straight_x,
        "straight_y_m": straight_y,
        "cross_section_m": [0.035, 0.105],
        "centre_m": centre.tolist(),
        "closed": True,
        "closure_gap_m": float(builder.gap),
        "radia_object_count": len(objects),
        "arc_max_segment_length_m": 0.004,
    }


def observation_points() -> np.ndarray:
    stations = np.linspace(-0.040, 0.040, 9)
    offsets_y = (-0.003, 0.0, 0.003)
    offsets_z = (-0.002, 0.0, 0.002)
    return np.asarray(
        [[x, y, z] for x in stations for y in offsets_y for z in offsets_z],
        dtype=float,
    )


def evaluate_cf(field, mesh: ng.Mesh, points: np.ndarray) -> np.ndarray:
    values = []
    for point in points:
        values.append(np.asarray(field(mesh(*map(float, point))), dtype=float))
    return np.asarray(values, dtype=float).reshape(-1, 3)


def relative_rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    denominator = float(np.sqrt(np.mean(np.sum(reference * reference, axis=1))))
    if denominator <= 0.0:
        raise RuntimeError("zero reference field in C-type comparison")
    return float(np.sqrt(np.mean(np.sum((candidate - reference) ** 2, axis=1))) / denominator)


def median_plane_projection(points: np.ndarray, field: np.ndarray) -> np.ndarray:
    """Project axial B onto the exact z-reflection symmetry.

    For z reflection, an axial vector transforms as ``(-Bx, -By, Bz)``.
    At ``z=0`` this is the parity-consistent two-sided trace used for comparing
    the three formulations.  In particular, ``curl(HCurl)`` has a continuous
    normal trace but its one-sided tangential trace is not a pointwise field
    contract on an element face.
    """

    keys = {tuple(np.round(point, 12)): index for index, point in enumerate(points)}
    projected = np.empty_like(field)
    axial_transform = np.asarray([-1.0, -1.0, 1.0])
    for index, point in enumerate(points):
        reflected = (float(point[0]), float(point[1]), float(-point[2]))
        other = keys.get(tuple(np.round(reflected, 12)))
        if other is None:
            raise RuntimeError(f"observation grid lacks reflected point {reflected}")
        projected[index] = 0.5 * (field[index] + axial_transform * field[other])
    return projected


def reflection_diagnostics(
    points: np.ndarray, field: np.ndarray, *, plane_tolerance: float = 1e-14
) -> dict[str, object]:
    """Separate physical off-plane symmetry from a one-sided face trace."""

    keys = {tuple(np.round(point, 12)): index for index, point in enumerate(points)}
    axial_transform = np.asarray([-1.0, -1.0, 1.0])
    off_plane = np.abs(points[:, 2]) > plane_tolerance
    median_plane = ~off_plane
    reflected_values = np.empty_like(field)
    for index, point in enumerate(points):
        reflected = (float(point[0]), float(point[1]), float(-point[2]))
        other = keys.get(tuple(np.round(reflected, 12)))
        if other is None:
            raise RuntimeError(f"observation grid lacks reflected point {reflected}")
        reflected_values[index] = axial_transform * field[other]

    off_reference = field[off_plane]
    off_reflected = reflected_values[off_plane]
    off_denominator = float(
        np.sqrt(np.mean(np.sum(off_reference * off_reference, axis=1)))
    )
    if off_denominator <= 0.0:
        raise RuntimeError("zero off-plane reference field in reflection diagnostic")
    off_difference = off_reference - off_reflected

    plane_values = field[median_plane]
    if plane_values.size == 0:
        raise RuntimeError("observation grid has no median-plane points")
    plane_norm = float(np.sqrt(np.mean(np.sum(plane_values * plane_values, axis=1))))
    tangential = np.linalg.norm(plane_values[:, :2], axis=1)
    return {
        "off_plane_relative_rms": float(
            np.sqrt(np.mean(np.sum(off_difference * off_difference, axis=1)))
            / off_denominator
        ),
        "off_plane_maximum_absolute_difference_T": float(
            np.max(np.linalg.norm(off_difference, axis=1))
        ),
        "off_plane_point_count": int(np.count_nonzero(off_plane)),
        "median_plane_one_sided_tangential_relative_rms": (
            float(np.sqrt(np.mean(tangential * tangential)) / plane_norm)
            if plane_norm > 0.0
            else 0.0
        ),
        "median_plane_one_sided_tangential_maximum_T": float(np.max(tangential)),
        "median_plane_point_count": int(np.count_nonzero(median_plane)),
        "trace_contract": (
            "The median-plane tangential value is a one-sided curl(HCurl) trace "
            "diagnostic, not a reflection-symmetry error."
        ),
    }


def pairwise_metrics(fields: dict[str, np.ndarray], selector: np.ndarray) -> dict[str, object]:
    pairs = {}
    names = tuple(fields)
    for left_index, left in enumerate(names):
        for right in names[left_index + 1 :]:
            left_values = fields[left][selector]
            right_values = fields[right][selector]
            pairs[f"{left}__vs__{right}"] = {
                "relative_rms": relative_rms(left_values, right_values),
                "maximum_absolute_difference_T": float(
                    np.max(np.linalg.norm(left_values - right_values, axis=1))
                ),
            }
    return pairs


def solve_hdiv(
    mesh: ng.Mesh,
    coil: int,
    material,
    *,
    nonlinear: bool,
    order: int,
    gram_eps: float,
    nonlinear_tolerance: float,
    nonlinear_maximum_iterations: int,
    points: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    started = time.perf_counter()
    source_h = rad.RadiaField(coil, "h")
    kwargs = {"bh_table": material} if nonlinear else {"mu_r": float(material)}
    with ng.TaskManager():
        result = vim.Solve(
            mesh,
            H_ext=source_h,
            order=order,
            gram_eps=gram_eps,
            tol=1e-8,
            maxit=12000,
            nl_tol=nonlinear_tolerance,
            nl_maxit=nonlinear_maximum_iterations,
            preconditioner="auto",
            **kwargs,
        )
        demag_h = vim.FieldFromSolution(result, points, algorithm="direct")
    coil_h = np.asarray(rad.Fld(coil, "h", points), dtype=float)
    field = MU0 * (np.asarray(demag_h, dtype=float) + coil_h)
    nonlinear_stats = dict(result.get("nonlinear_solve_stats", {}))
    if nonlinear:
        nonlinear_stats["converged"] = bool(
            nonlinear_stats.get("nonlinear_converged_final_stage", False)
        )
    return field, {
        "formulation": "HDiv-MMM BDM%d" % order,
        "open_boundary": "Coulomb charge Gram; iron-only mesh",
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(result["ndof"]),
        "linear_iterations": int(result["iters"]),
        "nonlinear": bool(result.get("nonlinear", nonlinear)),
        "nonlinear_stats": nonlinear_stats,
        "preconditioner": result.get("preconditioner"),
        "gram_eps": float(gram_eps),
        "runtime_s": float(time.perf_counter() - started),
    }


def solve_reduced_a(
    mesh: ng.Mesh,
    coil: int,
    material,
    *,
    nonlinear: bool,
    order: int,
    linear_solver: str,
    nonlinear_tolerance: float,
    nonlinear_maximum_iterations: int,
    nonlinear_verbose: bool,
    kelvin_center: tuple[float, float, float],
    kelvin_radius: float,
    points: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    started = time.perf_counter()
    solver = VectorPotentialSolver(
        mesh,
        iron_domains="iron",
        mu_r=float(material) if not nonlinear else 1000.0,
        order=order,
        kelvin_region="kelvin",
        kelvin_radius=kelvin_radius,
        kelvin_center=kelvin_center,
    )
    solver.set_source_cf(rad.RadiaField(coil, "b"))
    with ng.TaskManager():
        if nonlinear:
            solution = solver.solve_nonlinear(
                material,
                tol=nonlinear_tolerance,
                maxiter=nonlinear_maximum_iterations,
                dirichlet="GND",
                verbose=nonlinear_verbose,
                solver=linear_solver,
            )
        else:
            solution = solver.solve_linear(
                dirichlet="GND", solver=linear_solver
            )
    field = evaluate_cf(solver.get_B(), mesh, points)
    return field, {
        "formulation": "HCurl reduced-A",
        "open_boundary": "periodic spherical Kelvin transform",
        "source_contract": "vacuum source removed; iron contrast RHS only",
        "kelvin_center_m": list(kelvin_center),
        "kelvin_radius_m": kelvin_radius,
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(solution.space.ndof),
        "nonlinear": nonlinear,
        "linear_solver": linear_solver,
        "nonlinear_stats": getattr(solver, "_last_nonlinear_stats", {}),
        "runtime_s": float(time.perf_counter() - started),
    }


def solve_omega(
    mesh: ng.Mesh,
    coil: int,
    material,
    *,
    nonlinear: bool,
    order: int,
    nonlinear_tolerance: float,
    nonlinear_maximum_iterations: int,
    nonlinear_verbose: bool,
    kelvin_center: tuple[float, float, float],
    kelvin_radius: float,
    points: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    started = time.perf_counter()
    solver = ScalarPotentialSolver(
        mesh,
        iron_domains="iron",
        mu_r=float(material) if not nonlinear else 1000.0,
        order=order,
        kelvin_region="kelvin",
        kelvin_radius=kelvin_radius,
        kelvin_center=kelvin_center,
    )
    solver.set_source_cf(rad.RadiaField(coil, "h"))
    with ng.TaskManager():
        if nonlinear:
            solution = solver.solve_nonlinear(
                material,
                tol=nonlinear_tolerance,
                maxiter=nonlinear_maximum_iterations,
                relax=0.3,
                dirichlet="GND",
                verbose=nonlinear_verbose,
            )
        else:
            solution = solver.solve_single_potential(dirichlet="GND")
    field = evaluate_cf(solver.get_B(), mesh, points)
    return field, {
        "formulation": "H1 Omega-reduced-Omega",
        "open_boundary": "periodic spherical Kelvin transform",
        "source_contract": "vacuum source removed; iron contrast RHS only",
        "kelvin_center_m": list(kelvin_center),
        "kelvin_radius_m": kelvin_radius,
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(solution.space.ndof),
        "nonlinear": nonlinear,
        "nonlinear_stats": getattr(solver, "_last_nonlinear_stats", {}),
        "runtime_s": float(time.perf_counter() - started),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bh-table", type=Path, default=DEFAULT_BH)
    parser.add_argument("--mode", choices=("linear", "nonlinear"), default="nonlinear")
    parser.add_argument("--hdiv-order", type=int, choices=(1, 2), default=2)
    parser.add_argument(
        "--hdiv-gram-eps",
        type=float,
        default=1.0e-14,
        help="HDiv ChargeGram relative tolerance for this accuracy validation",
    )
    parser.add_argument("--fem-order", type=int, default=2)
    parser.add_argument(
        "--reduced-a-solver",
        choices=("direct", "bddc", "ams", "auto"),
        default="direct",
    )
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--nonlinear-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--nonlinear-maximum-iterations", type=int, default=80)
    parser.add_argument("--nonlinear-verbose", action="store_true")
    parser.add_argument(
        "--primary-only",
        action="store_true",
        help="run only the primary HDiv-MMM versus Omega-reduced-Omega comparison",
    )
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--gap-core-half-length", type=float, default=0.010)
    parser.add_argument("--relative-rms-tolerance", type=float, default=0.03)
    options = parser.parse_args()
    if options.threads > 0:
        ng.SetNumThreads(options.threads)
    if options.fem_order < 1:
        raise ValueError("--fem-order must be positive")
    if not 0.0 < options.hdiv_gram_eps < 1.0:
        raise ValueError("--hdiv-gram-eps must lie strictly between 0 and 1")
    if options.nonlinear_tolerance <= 0.0:
        raise ValueError("--nonlinear-tolerance must be positive")
    if options.nonlinear_maximum_iterations < 1:
        raise ValueError("--nonlinear-maximum-iterations must be positive")

    mesh_dir = options.mesh_dir.resolve()
    mesh_report_path = mesh_dir / "mesh_result.json"
    iron_vol = mesh_dir / "iron.vol"
    kelvin_vol = mesh_dir / "kelvin_domain.vol"
    for path in (mesh_report_path, iron_vol, kelvin_vol, options.bh_table):
        if not path.is_file():
            raise FileNotFoundError(path)
    mesh_report = json.loads(mesh_report_path.read_text(encoding="utf-8"))
    if not mesh_report.get("passed", False):
        raise RuntimeError("mesh_result.json is not a passing Cubit mesh contract")
    if mesh_report["artifacts"]["iron_vol_sha256"] != sha256(iron_vol):
        raise RuntimeError("iron.vol hash differs from mesh_result.json")
    if mesh_report["artifacts"]["kelvin_domain_vol_sha256"] != sha256(kelvin_vol):
        raise RuntimeError("kelvin_domain.vol hash differs from mesh_result.json")

    nonlinear = options.mode == "nonlinear"
    material = np.loadtxt(options.bh_table, dtype=float)[:, :2].tolist() if nonlinear else options.mu_r
    rad.UtiDelAll()
    coil, coil_manifest = build_coil()
    points = observation_points()
    iron_mesh = ng.Mesh(str(iron_vol))
    kelvin_mesh = ng.Mesh(str(kelvin_vol))
    if not has_kelvin_identification(kelvin_mesh):
        raise RuntimeError("kelvin_domain.vol has no periodic point identification")
    offset = np.asarray(detect_kelvin_offset(kelvin_mesh), dtype=float)
    physical_center = np.asarray(
        mesh_report["kelvin_physical_center_m"], dtype=float)
    kelvin_center = tuple((physical_center + offset).tolist())
    kelvin_radius = float(mesh_report["kelvin_radius_m"])

    fields = {}
    diagnostics = {}
    hdiv_checkpoint = options.output.with_suffix(".hdiv-checkpoint.json")
    hdiv_contract = {
        "mode": options.mode,
        "hdiv_order": int(options.hdiv_order),
        "hdiv_gram_eps": float(options.hdiv_gram_eps),
        "nonlinear_tolerance": float(options.nonlinear_tolerance),
        "nonlinear_maximum_iterations": int(options.nonlinear_maximum_iterations),
        "iron_vol_sha256": sha256(iron_vol),
        "bh_table_sha256": None if not nonlinear else sha256(options.bh_table),
        "observation_points": points.tolist(),
        "implementation_sha256": {
            "vim_assembly": sha256(Path(vim.__file__).resolve().parent / "_vim.py"),
            "vim_solve": sha256(Path(vim.__file__).resolve().parent / "_solve.py"),
            "radia_pybind": sha256(Path(radia_native.__file__).resolve()),
            "coil_builder": sha256(
                Path(sys.modules[CoilBuilder.__module__].__file__).resolve()
            ),
        },
    }
    resumed = read_checkpoint(hdiv_checkpoint, hdiv_contract) if options.resume else None
    if resumed is None:
        progress("engine_start", engine="hdiv_mmm", mesh_elements=int(iron_mesh.ne))
        fields["hdiv_mmm"], diagnostics["hdiv_mmm"] = solve_hdiv(
            iron_mesh,
            coil,
            material,
            nonlinear=nonlinear,
            order=options.hdiv_order,
            gram_eps=options.hdiv_gram_eps,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
            points=points,
        )
        write_checkpoint(
            hdiv_checkpoint, hdiv_contract, fields["hdiv_mmm"], diagnostics["hdiv_mmm"]
        )
        progress("engine_complete", engine="hdiv_mmm",
                 runtime_s=diagnostics["hdiv_mmm"]["runtime_s"],
                 checkpoint=str(hdiv_checkpoint))
    else:
        fields["hdiv_mmm"], diagnostics["hdiv_mmm"] = resumed
        diagnostics["hdiv_mmm"]["resumed_from_checkpoint"] = True
        progress("engine_resumed", engine="hdiv_mmm", checkpoint=str(hdiv_checkpoint))
    fem_contract = {
        "mode": options.mode,
        "fem_order": int(options.fem_order),
        "nonlinear_tolerance": float(options.nonlinear_tolerance),
        "nonlinear_maximum_iterations": int(options.nonlinear_maximum_iterations),
        "kelvin_domain_vol_sha256": sha256(kelvin_vol),
        "bh_table_sha256": None if not nonlinear else sha256(options.bh_table),
        "observation_points": points.tolist(),
        "implementation_sha256": {
            "radia_pybind": sha256(Path(radia_native.__file__).resolve()),
            "coil_builder": sha256(
                Path(sys.modules[CoilBuilder.__module__].__file__).resolve()
            ),
        },
    }
    reduced_a_checkpoint = options.output.with_suffix(
        ".reduced-a-checkpoint.json")
    reduced_a_contract = {
        **fem_contract,
        "engine": "reduced_a",
        "linear_solver": options.reduced_a_solver,
        "implementation_sha256": sha256(
            Path(sys.modules[VectorPotentialSolver.__module__].__file__).resolve()
        ),
    }
    if not options.primary_only:
        resumed = (
            read_checkpoint(reduced_a_checkpoint, reduced_a_contract)
            if options.resume else None
        )
        if resumed is None:
            progress("engine_start", engine="reduced_a",
                     mesh_elements=int(kelvin_mesh.ne))
            fields["reduced_a"], diagnostics["reduced_a"] = solve_reduced_a(
                kelvin_mesh,
                coil,
                material,
                nonlinear=nonlinear,
                order=options.fem_order,
                linear_solver=options.reduced_a_solver,
                nonlinear_tolerance=options.nonlinear_tolerance,
                nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
                nonlinear_verbose=options.nonlinear_verbose,
                kelvin_center=kelvin_center,
                kelvin_radius=kelvin_radius,
                points=points,
            )
            write_checkpoint(
                reduced_a_checkpoint,
                reduced_a_contract,
                fields["reduced_a"],
                diagnostics["reduced_a"],
            )
            progress("engine_complete", engine="reduced_a",
                     runtime_s=diagnostics["reduced_a"]["runtime_s"],
                     checkpoint=str(reduced_a_checkpoint))
        else:
            fields["reduced_a"], diagnostics["reduced_a"] = resumed
            diagnostics["reduced_a"]["resumed_from_checkpoint"] = True
            progress("engine_resumed", engine="reduced_a",
                     checkpoint=str(reduced_a_checkpoint))

    omega_checkpoint = options.output.with_suffix(".omega-checkpoint.json")
    omega_contract = {
        **fem_contract,
        "engine": "omega_reduced_omega",
        "implementation_sha256": sha256(
            Path(sys.modules[ScalarPotentialSolver.__module__].__file__).resolve()
        ),
    }
    resumed = (
        read_checkpoint(omega_checkpoint, omega_contract)
        if options.resume else None
    )
    if resumed is None:
        progress("engine_start", engine="omega_reduced_omega",
                 mesh_elements=int(kelvin_mesh.ne))
        fields["omega_reduced_omega"], diagnostics["omega_reduced_omega"] = solve_omega(
            kelvin_mesh,
            coil,
            material,
            nonlinear=nonlinear,
            order=options.fem_order,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
            nonlinear_verbose=options.nonlinear_verbose,
            kelvin_center=kelvin_center,
            kelvin_radius=kelvin_radius,
            points=points,
        )
        write_checkpoint(
            omega_checkpoint,
            omega_contract,
            fields["omega_reduced_omega"],
            diagnostics["omega_reduced_omega"],
        )
        progress("engine_complete", engine="omega_reduced_omega",
                 runtime_s=diagnostics["omega_reduced_omega"]["runtime_s"],
                 checkpoint=str(omega_checkpoint))
    else:
        fields["omega_reduced_omega"], diagnostics["omega_reduced_omega"] = resumed
        diagnostics["omega_reduced_omega"]["resumed_from_checkpoint"] = True
        progress("engine_resumed", engine="omega_reduced_omega",
                 checkpoint=str(omega_checkpoint))

    all_points = np.ones(points.shape[0], dtype=bool)
    gap_core = np.abs(points[:, 0]) <= float(options.gap_core_half_length) + 1e-14
    if not np.any(gap_core):
        raise RuntimeError("--gap-core-half-length selected no observation points")
    projected_fields = {
        name: median_plane_projection(points, value)
        for name, value in fields.items()
    }
    raw_pairs = pairwise_metrics(fields, all_points)
    projected_pairs = pairwise_metrics(projected_fields, all_points)
    core_pairs = pairwise_metrics(projected_fields, gap_core)
    reflection = {
        name: reflection_diagnostics(points, value)
        for name, value in fields.items()
    }
    primary_pair_name = "hdiv_mmm__vs__omega_reduced_omega"
    if primary_pair_name not in core_pairs:
        raise RuntimeError("primary HDiv-MMM versus Omega comparison is missing")
    primary_relative_rms = float(core_pairs[primary_pair_name]["relative_rms"])
    maximum_relative_rms = max(row["relative_rms"] for row in core_pairs.values())
    nonlinear_converged = all(
        bool(row.get("nonlinear_stats", {}).get("converged", False))
        for row in diagnostics.values()
    ) if nonlinear else True
    primary_accuracy_passed = primary_relative_rms <= options.relative_rms_tolerance
    all_pairwise_within_tolerance = (
        maximum_relative_rms <= options.relative_rms_tolerance
    )
    passed = (
        primary_accuracy_passed
        and nonlinear_converged
    )
    try:
        radia_version = importlib.metadata.version("radia")
    except importlib.metadata.PackageNotFoundError:
        radia_version = "editable-unversioned"
    output = {
        "schema": "radia.validation.c-type-formulation-comparison.v2",
        "passed": passed,
        "machine": platform.node(),
        "python": sys.version,
        "radia_version": radia_version,
        "radia_module": str(rad.__file__),
        "mode": options.mode,
        "nonlinear_converged": nonlinear_converged,
        "primary_accuracy_passed": primary_accuracy_passed,
        "all_pairwise_within_tolerance": all_pairwise_within_tolerance,
        "comparison_contract": {
            "cad_authority_sha256": mesh_report["cad_sha256"],
            "coilbuilder_source_shared": True,
            "bh_table_shared": nonlinear,
            "bh_interpolation_shared": (
                "monotone PCHIP B(H), vacuum-slope continuation beyond H_max"
                if nonlinear else None
            ),
            "observation_points_shared": True,
            "hdiv_air_mesh_forbidden": True,
            "hdiv_gram_eps": float(options.hdiv_gram_eps),
            "fem_periodic_kelvin_mesh_shared": True,
            "finite_outer_air_box_forbidden": True,
            "quantity": "gauge-invariant magnetic flux density B in tesla",
            "fixed_mesh_equality_claimed": False,
            "primary_pair": ["hdiv_mmm", "omega_reduced_omega"],
            "reduced_a_role": (
                "not run (--primary-only)"
                if options.primary_only
                else "independent third-formulation cross-check"
            ),
            "acceptance": (
                "parity-projected pairwise B convergence in the useful gap core; "
                "off-plane reflection and raw/full-fringe values remain mandatory "
                "diagnostics"
            ),
        },
        "engine_checkpoint_contracts": {
            name: contract
            for name, contract in (
                ("hdiv_mmm", hdiv_contract),
                ("reduced_a", reduced_a_contract),
                ("omega_reduced_omega", omega_contract),
            )
            if name in fields
        },
        "mesh_result": str(mesh_report_path),
        "mesh_result_sha256": sha256(mesh_report_path),
        "bh_table": None if not nonlinear else str(options.bh_table.resolve()),
        "bh_table_sha256": None if not nonlinear else sha256(options.bh_table),
        "coil": coil_manifest,
        "observation_points_m": points.tolist(),
        "engines": diagnostics,
        "fields_T": {name: value.tolist() for name, value in fields.items()},
        "median_plane_projected_fields_T": {
            name: value.tolist() for name, value in projected_fields.items()
        },
        "reflection_diagnostics": reflection,
        "pairwise_raw_full_tube": raw_pairs,
        "pairwise_median_projected_full_tube": projected_pairs,
        "pairwise_median_projected_gap_core": core_pairs,
        "gap_core_half_length_m": float(options.gap_core_half_length),
        "gap_core_point_count": int(np.count_nonzero(gap_core)),
        "primary_gap_core_relative_rms": primary_relative_rms,
        "primary_gap_core_metrics": core_pairs[primary_pair_name],
        "maximum_gap_core_pairwise_relative_rms": maximum_relative_rms,
        "relative_rms_tolerance": options.relative_rms_tolerance,
    }
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    progress("validation_complete", passed=bool(passed), output=str(options.output),
             primary_gap_core_relative_rms=primary_relative_rms,
             maximum_gap_core_pairwise_relative_rms=maximum_relative_rms)
    if not passed:
        raise RuntimeError(
            f"HDiv-MMM versus Omega C-type gate failed: {primary_relative_rms:.6e} > "
            f"{options.relative_rms_tolerance:.6e}; see {options.output}"
        )


if __name__ == "__main__":
    main()
