"""Compare three formulations on one C-type electromagnet.

HDiv-MMM, HCurl reduced-A, and H1 TOSCA-style mixed total/reduced Omega share
the exact Cubit/ACIS iron authority, one solid ``CoilBuilder`` excitation, one
B-H table, and one set of physical observation points. HDiv-MMM intentionally
uses the iron-only mesh because its Coulomb Gram is the open-boundary operator;
the FEM formulations use the same periodic Kelvin mesh. No finite outer air
box is part of this comparison.
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
from datetime import datetime, timezone
from pathlib import Path

# A validation invocation from a source checkout must evaluate the source under
# test, not an unrelated installed wheel that happens to be on ``PATH``.
HERE = Path(__file__).resolve().parent
SOURCE_ROOT = HERE.parents[1] / "src"
SOURCE_PACKAGE = SOURCE_ROOT / "radia"
if (SOURCE_PACKAGE / "__init__.py").is_file():
    sys.path.insert(0, str(SOURCE_ROOT))
    # A result-bearing notebook may already have imported the installed wheel
    # before it loads this validation runner.  Make the source package explicit
    # in that case too, so the runner cannot silently combine its source files
    # with an unrelated wheel-only module set.
    imported_radia = sys.modules.get("radia")
    if imported_radia is not None:
        package_paths = getattr(imported_radia, "__path__", None)
        if package_paths is not None and str(SOURCE_PACKAGE) not in package_paths:
            package_paths.insert(0, str(SOURCE_PACKAGE))

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
from radia.electromagnet_validation import (
    require_static_electromagnet_three_engine_contract,
)
from radia.static_electromagnet import (
    StaticElectromagnetMixedDomain,
    solve_static_electromagnet_mixed_total_reduced_omega,
)
from radia.vector_potential_solver import VectorPotentialSolver


MU0 = 4.0e-7 * math.pi
DEFAULT_BH = (
    Path(rad.__file__).resolve().parent
    / "panels"
    / "samples"
    / "em_sample_bh.txt"
)
MIXED_DOMAIN = StaticElectromagnetMixedDomain(
    reduced_materials=("air",),
    total_materials=("iron", "kelvin"),
    nonlinear_materials=("iron",),
    reduced_total_interface="iron_air_interface",
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
    image: str | None = None,
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
            image=image,
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
        "formulation": "HDiv-MMM",
        "discretization": "BDM%d" % order,
        "open_boundary": "Coulomb charge Gram; iron-only mesh",
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(result["ndof"]),
        "linear_iterations": int(result["iters"]),
        "nonlinear": bool(result.get("nonlinear", nonlinear)),
        "nonlinear_stats": nonlinear_stats,
        "preconditioner": result.get("preconditioner"),
        "gram_eps": float(gram_eps),
        "image": image,
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
    relax: float,
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
                relax=relax,
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
        "relax": float(relax),
        "bh_interpolation": (
            "exact inverse of shared monotone PCHIP B(H), with vacuum-slope "
            "continuation" if nonlinear else None
        ),
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
    source_trace_tolerance: float,
) -> tuple[np.ndarray, dict[str, object]]:
    """Run the TOSCA-style total/reduced Omega route on the Kelvin mesh.

    The physical air contains the CoilBuilder source and is the reduced
    region.  Iron and the Kelvin exterior are total-potential regions, so the
    source field is never numerically cancelled in high-permeability iron.
    """
    started = time.perf_counter()
    source_h = rad.RadiaField(coil, "h")
    with ng.TaskManager():
        result = solve_static_electromagnet_mixed_total_reduced_omega(
            mesh,
            source_h,
            MIXED_DOMAIN,
            kelvin_radius,
            kelvin_center,
            order=order,
            linear_mu_r_by_material=None if nonlinear else {"iron": float(material)},
            bh_table=material if nonlinear else None,
            source_trace_tolerance=source_trace_tolerance,
            source_potential_contract="total_hodge",
            nonlinear_tolerance=nonlinear_tolerance,
            nonlinear_max_iterations=nonlinear_maximum_iterations,
            nonlinear_relaxation=0.3,
        )
    field = evaluate_cf(result["B_cf"], mesh, points)
    source_trace = result["static_electromagnet_contract"]["source_trace"]
    return field, {
        "formulation": "H1 TOSCA mixed total/reduced Omega",
        "open_boundary": "periodic spherical Kelvin transform",
        "source_contract": (
            "exact Radia H is restricted to physical air; the source/iron "
            "interface uses a projected scalar trace; iron and Kelvin use "
            "total Omega"
        ),
        "source_trace": {
            "iron_interface_boundary": "iron_air_interface",
            "projection_order": int(source_trace["projection_order"]),
            "iron_relative_harmonic_norm": float(
                source_trace["iron_relative_harmonic_norm"]
            ),
            "kelvin_interface_boundary": "kelvin_int",
            "kelvin_relative_tangential_residual": float(
                source_trace["kelvin_relative_tangential_residual"]
            ),
            "relative_tolerance": float(source_trace_tolerance),
            "cut_policy": (
                "the iron volume Hodge split retains the linked-source "
                "harmonic/cohomology field; the Kelvin jump carries the "
                "orientation-reversed exact source 0-form"
            ),
        },
        "kelvin_center_m": list(kelvin_center),
        "kelvin_radius_m": kelvin_radius,
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(result["fes"].ndof),
        "nonlinear": nonlinear,
        "nonlinear_stats": result.get("nonlinear_stats", {}),
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
    parser.add_argument("--reduced-a-relax", type=float, default=0.1)
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--nonlinear-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--nonlinear-maximum-iterations", type=int, default=80)
    parser.add_argument("--nonlinear-verbose", action="store_true")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--gap-core-half-length", type=float, default=0.010)
    parser.add_argument("--relative-rms-tolerance", type=float, default=0.03)
    parser.add_argument(
        "--source-trace-tolerance",
        type=float,
        default=0.05,
        help=(
            "maximum relative tangential residual for each projected source "
            "trace; a larger residual requires an explicit cut/cohomology source"
        ),
    )
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
    if not 0.0 < options.reduced_a_relax <= 1.0:
        raise ValueError("--reduced-a-relax must lie in (0, 1]")
    if not 0.0 < options.source_trace_tolerance < 1.0:
        raise ValueError("--source-trace-tolerance must lie in (0, 1)")

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
        "source_trace_tolerance": float(options.source_trace_tolerance),
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
        "relax": float(options.reduced_a_relax),
        "implementation_sha256": sha256(
            Path(sys.modules[VectorPotentialSolver.__module__].__file__).resolve()
        ),
    }
    resumed = (
        read_checkpoint(reduced_a_checkpoint, reduced_a_contract)
        if options.resume else None
    )
    if resumed is None:
        progress("engine_start", engine="reduced_a", mesh_elements=int(kelvin_mesh.ne))
        fields["reduced_a"], diagnostics["reduced_a"] = solve_reduced_a(
            kelvin_mesh,
            coil,
            material,
            nonlinear=nonlinear,
            order=options.fem_order,
            linear_solver=options.reduced_a_solver,
            relax=options.reduced_a_relax,
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
        progress("engine_resumed", engine="reduced_a", checkpoint=str(reduced_a_checkpoint))

    mixed_omega_checkpoint = options.output.with_suffix(
        ".mixed-total-reduced-omega-checkpoint.json"
    )
    mixed_omega_contract = {
        **fem_contract,
        "engine": "mixed_total_reduced_omega",
        "implementation_sha256": sha256(
            Path(sys.modules[
                solve_static_electromagnet_mixed_total_reduced_omega.__module__
            ].__file__).resolve()
        ),
    }
    resumed = (
        read_checkpoint(mixed_omega_checkpoint, mixed_omega_contract)
        if options.resume else None
    )
    if resumed is None:
        progress("engine_start", engine="mixed_total_reduced_omega",
                 mesh_elements=int(kelvin_mesh.ne))
        fields["mixed_total_reduced_omega"], diagnostics[
            "mixed_total_reduced_omega"
        ] = solve_omega(
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
            source_trace_tolerance=options.source_trace_tolerance,
        )
        write_checkpoint(
            mixed_omega_checkpoint,
            mixed_omega_contract,
            fields["mixed_total_reduced_omega"],
            diagnostics["mixed_total_reduced_omega"],
        )
        progress("engine_complete", engine="mixed_total_reduced_omega",
                 runtime_s=diagnostics["mixed_total_reduced_omega"]["runtime_s"],
                 checkpoint=str(mixed_omega_checkpoint))
    else:
        fields["mixed_total_reduced_omega"], diagnostics[
            "mixed_total_reduced_omega"
        ] = resumed
        diagnostics["mixed_total_reduced_omega"]["resumed_from_checkpoint"] = True
        progress("engine_resumed", engine="mixed_total_reduced_omega",
                 checkpoint=str(mixed_omega_checkpoint))

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
    hdiv_mixed_pair_name = "hdiv_mmm__vs__mixed_total_reduced_omega"
    if hdiv_mixed_pair_name not in core_pairs:
        raise RuntimeError("HDiv-MMM versus mixed total/reduced Omega comparison is missing")
    hdiv_mixed_relative_rms = float(core_pairs[hdiv_mixed_pair_name]["relative_rms"])
    maximum_relative_rms = max(row["relative_rms"] for row in core_pairs.values())
    nonlinear_converged = all(
        bool(row.get("nonlinear_stats", {}).get("converged", False))
        for row in diagnostics.values()
    ) if nonlinear else True
    hdiv_mixed_accuracy_passed = (
        hdiv_mixed_relative_rms <= options.relative_rms_tolerance
    )
    all_pairwise_within_tolerance = (
        maximum_relative_rms <= options.relative_rms_tolerance
    )
    require_static_electromagnet_three_engine_contract(diagnostics)
    passed = all_pairwise_within_tolerance and nonlinear_converged
    try:
        radia_version = importlib.metadata.version("radia")
    except importlib.metadata.PackageNotFoundError:
        radia_version = "editable-unversioned"
    output = {
        "schema": "radia.validation.c-type-formulation-comparison.v4",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "passed": passed,
        "machine": platform.node(),
        "python": sys.version,
        "radia_version": radia_version,
        "radia_module": str(rad.__file__),
        "mode": options.mode,
        "nonlinear_converged": nonlinear_converged,
        "hdiv_mixed_accuracy_passed": hdiv_mixed_accuracy_passed,
        "all_pairwise_within_tolerance": all_pairwise_within_tolerance,
        "comparison_contract": {
            "cad_authority_sha256": mesh_report["cad_sha256"],
            "coilbuilder_source_shared": True,
            "bh_table_shared": nonlinear,
            "bh_interpolation_shared": (
                "monotone PCHIP B(H), exact reduced-A inversion, and "
                "vacuum-slope continuation beyond H_max"
                if nonlinear else None
            ),
            "observation_points_shared": True,
            "hdiv_air_mesh_forbidden": True,
            "hdiv_gram_eps": float(options.hdiv_gram_eps),
            "source_trace_tolerance": float(options.source_trace_tolerance),
            "fem_periodic_kelvin_mesh_shared": True,
            "finite_outer_air_box_forbidden": True,
            "quantity": "gauge-invariant magnetic flux density B in tesla",
            "fixed_mesh_equality_claimed": False,
            "static_electromagnet_formulation_contract": require_static_electromagnet_three_engine_contract(
                diagnostics
            ),
            "acceptance": (
                "all three formulations are mandatory and every parity-projected B "
                "comparison in the useful gap core to pass; off-plane reflection "
                "and raw/full-fringe values remain mandatory diagnostics"
            ),
        },
        "engine_checkpoint_contracts": {
            name: contract
            for name, contract in (
                ("hdiv_mmm", hdiv_contract),
                ("reduced_a", reduced_a_contract),
                ("mixed_total_reduced_omega", mixed_omega_contract),
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
        "hdiv_mixed_gap_core_relative_rms": hdiv_mixed_relative_rms,
        "hdiv_mixed_gap_core_metrics": core_pairs[hdiv_mixed_pair_name],
        "maximum_gap_core_pairwise_relative_rms": maximum_relative_rms,
        "relative_rms_tolerance": options.relative_rms_tolerance,
    }
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    progress("validation_complete", passed=bool(passed), output=str(options.output),
             hdiv_mixed_gap_core_relative_rms=hdiv_mixed_relative_rms,
             maximum_gap_core_pairwise_relative_rms=maximum_relative_rms)
    if not passed:
        failures = []
        if not nonlinear_converged:
            failures.append("at least one selected nonlinear formulation did not converge")
        if not all_pairwise_within_tolerance:
            failures.append(
                f"all pairwise relative RMS {maximum_relative_rms:.6e} exceeds "
                f"{options.relative_rms_tolerance:.6e}"
            )
        raise RuntimeError("; ".join(failures) + f"; see {options.output}")


if __name__ == "__main__":
    main()
