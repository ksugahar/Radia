"""TEAM 28 CoilBuilder + HCurl eddy-bubble coupling validation.

The two counter-wound TEAM 28 winding packs are built with the public
``CoilBuilder`` API.  Their incident vector potential and flux density are
first checked against an independent Gauss-integrated circular-filament
reference.  The verified fields then drive the existing p=6 HCurl
eddy-bubble/VIM model and its passive CLN export for MATLAB.

The default case is validation-class work.  It intentionally keeps the coil
source, field cross-check, spatial reduction, solve, force evaluation, and
MATLAB exchange export in one reproducible artifact.
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import radia
import radia.vim as vim
from radia.coil_builder import CoilBuilder

from team28_hcurl_vim_force import (
    FREQUENCY_HZ,
    REFERENCE_COIL_CURRENT_A,
    SIGMA_AL,
    TARGET_PHYSICAL_FORCE_N,
    _build_eddy_basis,
    _external_coil_fields,
    _rectangular_coil_fields,
)

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_OUTPUT = HERE / "team28_coilbuilder_eddy_bubble_results.json"
DEFAULT_EXCHANGE = HERE / "team28_coilbuilder_hcurl_eddy_cln.json"

COIL_SPECS = (
    {
        "name": "inner",
        "mean_radius_m": 0.041,
        "radial_width_m": 0.028,
        "axial_height_m": 0.052,
        "turns": 960.0,
        "direction": 1.0,
    },
    {
        "name": "outer",
        "mean_radius_m": 0.0875,
        "radial_width_m": 0.015,
        "axial_height_m": 0.052,
        "turns": 576.0,
        "direction": -1.0,
    },
)


def build_team28_coils(coil_current_A=REFERENCE_COIL_CURRENT_A):
    """Build the two closed, counter-wound TEAM 28 winding packs."""

    current = float(coil_current_A)
    if not np.isfinite(current):
        raise ValueError("coil_current_A must be finite")
    builders = []
    for spec in COIL_SPECS:
        radius = spec["mean_radius_m"]
        height = spec["axial_height_m"]
        ampere_turns = spec["direction"] * spec["turns"] * current
        builder = (
            CoilBuilder(current=ampere_turns)
            .set_start([radius, 0.0, -0.5 * height])
            .set_cross_section(spec["radial_width_m"], height)
            .add_arc(radius=radius, arc_angle=360.0)
        )
        builders.append(builder)
    return tuple(builders)


def coilbuilder_fields(
    points,
    *,
    coil_current_A=REFERENCE_COIL_CURRENT_A,
    arc_max_segment_length_m=0.002,
    height_offset_m=0.0,
):
    """Evaluate CoilBuilder solid-current A and B, then release Radia objects."""

    sample_points = np.asarray(points, dtype=float)
    if sample_points.ndim != 2 or sample_points.shape[1] != 3:
        raise ValueError("points must have shape (n, 3)")
    evaluation_points = sample_points.copy()
    evaluation_points[:, 2] += float(height_offset_m)
    builders = build_team28_coils(coil_current_A)
    owned_objects = []
    container = None
    try:
        for builder in builders:
            owned_objects.extend(
                builder.to_radia(arc_max_segment_length=float(arc_max_segment_length_m))
            )
        container = radia.ObjCnt(owned_objects)
        vector_potential = np.asarray(
            radia.Fld(container, "a", evaluation_points.tolist()), dtype=float
        )
        flux_density = np.asarray(
            radia.Fld(container, "b", evaluation_points.tolist()), dtype=float
        )
    finally:
        if container is not None:
            radia.UtiDel(container)
        else:
            for obj in owned_objects:
                radia.UtiDel(obj)
    report = {
        "coil_count": len(builders),
        "source_object_count": len(owned_objects),
        "all_paths_closed": all(builder.is_closed for builder in builders),
        "closure_gaps_m": [float(builder.gap) for builder in builders],
        "arc_max_segment_length_m": float(arc_max_segment_length_m),
        "height_offset_m": float(height_offset_m),
        "coils": [
            {
                **spec,
                "ampere_turns_A": float(spec["direction"] * spec["turns"] * coil_current_A),
            }
            for spec in COIL_SPECS
        ],
    }
    return vector_potential, flux_density, report


def compare_coil_fields(
    points,
    *,
    coil_current_A=REFERENCE_COIL_CURRENT_A,
    arc_max_segment_length_m=0.002,
    height_offset_m=0.0,
    cross_check_indices=None,
    reference_radial_order=10,
    reference_axial_order=14,
):
    """Cross-check CoilBuilder fields against the filament-quadrature source."""

    builder_a, builder_b, report = coilbuilder_fields(
        points,
        coil_current_A=coil_current_A,
        arc_max_segment_length_m=arc_max_segment_length_m,
        height_offset_m=height_offset_m,
    )
    if cross_check_indices is None:
        check_points = np.asarray(points, dtype=float)
        check_builder_a = builder_a
        check_builder_b = builder_b
    else:
        indices = np.asarray(cross_check_indices, dtype=int)
        check_points = np.asarray(points, dtype=float)[indices]
        check_builder_a = builder_a[indices]
        check_builder_b = builder_b[indices]
    if reference_radial_order == 10 and reference_axial_order == 14:
        reference_a, reference_b = _external_coil_fields(
            check_points,
            coil_current=coil_current_A,
            height_offset_m=height_offset_m,
        )
    else:
        shifted_points = check_points.copy()
        shifted_points[:, 2] += float(height_offset_m)
        scale = float(coil_current_A) / REFERENCE_COIL_CURRENT_A
        a1, b1 = _rectangular_coil_fields(
            shifted_points,
            0.041,
            0.028,
            0.052,
            960.0 * REFERENCE_COIL_CURRENT_A * scale,
            radial_order=reference_radial_order,
            axial_order=reference_axial_order,
        )
        a2, b2 = _rectangular_coil_fields(
            shifted_points,
            0.0875,
            0.015,
            0.052,
            -576.0 * REFERENCE_COIL_CURRENT_A * scale,
            radial_order=reference_radial_order,
            axial_order=reference_axial_order,
        )
        reference_a = a1 + a2
        reference_b = b1 + b2
    tiny = np.finfo(float).tiny
    a_relative_l2 = float(
        np.linalg.norm(check_builder_a - reference_a)
        / max(np.linalg.norm(reference_a), tiny)
    )
    b_relative_l2 = float(
        np.linalg.norm(check_builder_b - reference_b)
        / max(np.linalg.norm(reference_b), tiny)
    )
    b_pointwise_error = np.linalg.norm(check_builder_b - reference_b, axis=1)
    b_reference_scale = max(float(np.max(np.linalg.norm(reference_b, axis=1))), tiny)
    report["field_cross_check"] = {
        "reference": "Gauss-integrated circular-filament winding pack",
        "sample_count": int(len(check_builder_a)),
        "reference_radial_order": int(reference_radial_order),
        "reference_axial_order": int(reference_axial_order),
        "vector_potential_relative_l2": a_relative_l2,
        "flux_density_relative_l2": b_relative_l2,
        "flux_density_max_pointwise_over_reference_max": float(
            np.max(b_pointwise_error) / b_reference_scale
        ),
    }
    return builder_a, builder_b, report


def _force_operator(current_basis, unit_external_b):
    return np.transpose(
        np.sum(
            current_basis.weights[None, :, None]
            * np.cross(
                current_basis.modes,
                np.conj(unit_external_b)[None, :, :],
            ),
            axis=1,
        ),
        (1, 0),
    )[:, :, None]


def _evaluate_force(force_operator, coefficients, coil_current_A):
    return 0.5 * np.real(
        np.einsum(
            "kab,a,b->k",
            force_operator,
            np.asarray(coefficients),
            np.asarray([coil_current_A], dtype=complex).conj(),
        )
    )


def _repo_path(path):
    return Path(path).resolve().relative_to(REPO_ROOT).as_posix()


def run_case(
    *,
    maxh_m=0.025,
    outer_quad=4,
    arc_max_segment_length_m=0.002,
    exchange_file=DEFAULT_EXCHANGE,
):
    """Run the coupled source, eddy-bubble, CLN, and force validation."""

    timings = {}
    stage = time.perf_counter()
    ng, mesh, fes, basis = _build_eddy_basis(maxh_m)
    timings["mesh_and_eddy_basis"] = time.perf_counter() - stage

    stage = time.perf_counter()
    external_a, external_b, coil_report = compare_coil_fields(
        basis.current_basis.points,
        arc_max_segment_length_m=arc_max_segment_length_m,
    )
    timings["coilbuilder_and_field_cross_check"] = time.perf_counter() - stage

    stage = time.perf_counter()
    with ng.TaskManager():
        interaction = basis.tet_volume_interaction(
            mesh,
            fes,
            degree=5,
            projection_quad=7,
            outer_quad=outer_quad,
            projection_tolerance=1.0e-10,
            materials="Al",
        )
    system = basis.assemble_vim(sigma=SIGMA_AL, interaction=interaction)
    unit_external_a = external_a / REFERENCE_COIL_CURRENT_A
    unit_external_b = external_b / REFERENCE_COIL_CURRENT_A
    rhs = vim.ExternalVectorPotentialRHS(basis.current_basis, unit_external_a)
    cln_model = vim.HCurlEddyCLNFromVIM(system, rhs)
    force_operator = _force_operator(basis.current_basis, unit_external_b)
    timings["vim_and_cln_assembly"] = time.perf_counter() - stage

    stage = time.perf_counter()
    s = 2.0j * np.pi * FREQUENCY_HZ
    currents_A = np.asarray([10.0, 20.0, 30.0])
    force_rows = []
    reference_coefficients = None
    for current_A in currents_A:
        coefficients = cln_model.solve_vector_potential_drive(s, current_A)
        if current_A == REFERENCE_COIL_CURRENT_A:
            reference_coefficients = coefficients
        force_N = _evaluate_force(force_operator, coefficients, current_A)
        force_rows.append(
            {
                "coil_current_A": float(current_A),
                "force_N": force_N.tolist(),
                "upward_force_N": float(force_N[2]),
                "force_per_current_squared_N_per_A2": float(force_N[2] / current_A**2),
            }
        )
    if reference_coefficients is None:
        raise RuntimeError("reference current was not included in the sweep")
    reference_force = np.asarray(force_rows[1]["force_N"])
    force_relative_error = float(
        abs(reference_force[2] - TARGET_PHYSICAL_FORCE_N) / TARGET_PHYSICAL_FORCE_N
    )
    transverse_ratio = float(
        np.linalg.norm(reference_force[:2]) / max(abs(reference_force[2]), np.finfo(float).tiny)
    )
    quadratic_values = np.asarray([row["force_per_current_squared_N_per_A2"] for row in force_rows])
    quadratic_spread = float(
        np.ptp(quadratic_values) / max(abs(np.mean(quadratic_values)), 1.0e-30)
    )
    joule_loss_W = float(
        0.5
        * np.real(
            np.vdot(
                reference_coefficients,
                cln_model.resistance @ reference_coefficients,
            )
        )
    )
    magnetic_energy_J = float(
        0.25
        * np.real(
            np.vdot(
                reference_coefficients,
                cln_model.inductance @ reference_coefficients,
            )
        )
    )
    vim.ExportHCurlEddyCLNJSON(
        cln_model,
        exchange_file,
        force_operator=force_operator,
        metadata={
            "benchmark": "TEAM 28 electrodynamic levitation",
            "frequency_hz": FREQUENCY_HZ,
            "coil_current_reference_A": REFERENCE_COIL_CURRENT_A,
            "conductivity_S_per_m": SIGMA_AL,
            "coil_source": "radia.coil_builder.CoilBuilder",
            "coil_count": 2,
            "winding_directions": "counter-wound",
            "force_operator_units": ("N per reduced coefficient and ampere port current"),
            "parent_space": "HCurl",
            "parent_order": 6,
        },
    )
    timings["harmonic_solve_force_and_export"] = time.perf_counter() - stage

    field_check = coil_report["field_cross_check"]
    interaction_report = interaction.diagnostics()
    cln_report = cln_model.diagnostics()
    checks = {
        "coilbuilder_created_two_closed_counter_wound_coils": bool(
            coil_report["coil_count"] == 2
            and coil_report["all_paths_closed"]
            and coil_report["coils"][0]["ampere_turns_A"] > 0.0
            and coil_report["coils"][1]["ampere_turns_A"] < 0.0
        ),
        "vector_potential_matches_filament_reference_below_0p05_percent": bool(
            field_check["vector_potential_relative_l2"] < 5.0e-4
        ),
        "flux_density_matches_filament_reference_below_0p1_percent": bool(
            field_check["flux_density_relative_l2"] < 1.0e-3
        ),
        "interaction_projection_residual_below_1e-10": bool(
            interaction_report["projection_relative_residual"] < 1.0e-10
        ),
        "reduced_inductance_block_positive": bool(cln_report["min_inductance_eigenvalue"] > 0.0),
        "cln_handoff_passive": bool(cln_report["passive"]),
        "upward_force_matches_reference_below_one_percent": bool(
            reference_force[2] > 0.0 and force_relative_error < 0.01
        ),
        "transverse_force_below_half_percent": bool(transverse_ratio < 0.005),
        "linear_model_force_scales_with_current_squared": bool(quadratic_spread < 1.0e-12),
        "runtime_reduction_below_one_percent": bool(basis.rank / fes.ndof < 0.01),
        "kernel_is_epsilon_free": bool(interaction_report["kernel_epsilon_m"] is None),
        "matlab_exchange_written": bool(Path(exchange_file).is_file()),
    }
    return {
        "timing_breakdown_s": timings,
        "problem": {
            "frequency_hz": FREQUENCY_HZ,
            "conductivity_S_per_m": SIGMA_AL,
            "disk_radius_m": 0.065,
            "disk_thickness_m": 0.003,
            "disk_bottom_m": 0.0108,
            "parent_space": "HCurl",
            "parent_order": 6,
            "mesh_maxh_m": float(maxh_m),
            "outer_quadrature_order": int(outer_quad),
        },
        "coil_builder": coil_report,
        "reduced_model": {
            "mesh_elements": int(mesh.ne),
            "parent_ndof": int(fes.ndof),
            "evrs_rank": int(basis.rank),
            "runtime_reduction_ratio": float(basis.rank / fes.ndof),
            "current_sample_count": int(basis.current_basis.n_samples),
            "interaction": interaction_report,
            "cln_handoff": cln_report,
            "eddy_bubble": basis.eddy_bubbling.diagnostics(),
        },
        "observables": {
            "target_upward_force_N": TARGET_PHYSICAL_FORCE_N,
            "reference_current_force_N": reference_force.tolist(),
            "upward_force_relative_error": force_relative_error,
            "transverse_force_ratio": transverse_ratio,
            "joule_loss_W": joule_loss_W,
            "magnetic_energy_J": magnetic_energy_J,
            "current_sweep": force_rows,
            "force_over_current_squared_relative_spread": quadratic_spread,
        },
        "checks": checks,
    }


def _make_success_artifact(details, output_file, exchange_file, duration_s):
    checks = {
        "ran_to_completion": True,
        "result_files_exist": bool(
            Path(exchange_file).is_file() and Path(output_file).parent.is_dir()
        ),
        "validation_passed": bool(all(details["checks"].values())),
        **details["checks"],
    }
    force_error = details["observables"]["upward_force_relative_error"]
    field_error = details["coil_builder"]["field_cross_check"]["flux_density_relative_l2"]
    return {
        "radia_version": getattr(radia, "__version__", "unknown"),
        "schema": "cae-ai-lab.solver-run.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "case": "TEAM 28 CoilBuilder coupled to HCurl eddy-bubble CLN",
        "solver": "radia-ngsolve",
        "source_artifact": _repo_path(__file__),
        "pass": bool(all(checks.values())),
        "run": {
            "command": ("python validation_test/maglev/" "team28_coilbuilder_eddy_bubble.py"),
            "workdir": ".",
            "exit_code": 0,
            "duration_s": float(duration_s),
        },
        "result_files": [
            _repo_path(output_file),
            _repo_path(exchange_file),
        ],
        "checks": checks,
        "tolerances": {
            "coil_field_max_rel": 1.0e-3,
            "force_max_rel": 1.0e-2,
            "transverse_force_ratio_max": 5.0e-3,
            "projection_residual_max": 1.0e-10,
        },
        "errors": {
            "max_rel": float(max(force_error, field_error)),
            "max_abs": float(
                abs(
                    details["observables"]["reference_current_force_N"][2] - TARGET_PHYSICAL_FORCE_N
                )
            ),
        },
        "tool_versions": {
            "python": platform.python_version(),
            "radia": getattr(radia, "__version__", "unknown"),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
            "numpy": np.__version__,
        },
        "timing_breakdown_s": details["timing_breakdown_s"],
        "verification": {
            "method": (
                "independent winding-pack field quadrature, passive CLN gates, "
                "stored TEAM force, symmetry, and current-squared scaling"
            ),
            "command": (
                "pytest -q tests/test_team28_coilbuilder_source.py "
                "validation_test/maglev/"
                "test_team28_coilbuilder_eddy_bubble_evidence.py"
            ),
        },
        "environment": {
            "platform": platform.platform(),
        },
        "details": details,
    }


def _make_failure_artifact(error, duration_s):
    return {
        "radia_version": getattr(radia, "__version__", "unknown"),
        "schema": "cae-ai-lab.solver-run.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "case": "TEAM 28 CoilBuilder coupled to HCurl eddy-bubble CLN",
        "solver": "radia-ngsolve",
        "pass": False,
        "run": {
            "command": ("python validation_test/maglev/" "team28_coilbuilder_eddy_bubble.py"),
            "workdir": ".",
            "exit_code": 1,
            "duration_s": float(duration_s),
        },
        "checks": {
            "ran_to_completion": False,
            "validation_passed": False,
        },
        "failure": {
            "stage": "solve",
            "message": f"{type(error).__name__}: {error}",
            "next_action": "inspect the failing source, basis, or CLN validation stage",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxh", type=float, default=0.025)
    parser.add_argument("--outer-quad", type=int, default=4)
    parser.add_argument("--arc-max-segment-length", type=float, default=0.002)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--export-model", type=Path, default=DEFAULT_EXCHANGE)
    args = parser.parse_args()

    started = time.perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.export_model.parent.mkdir(parents=True, exist_ok=True)
    try:
        details = run_case(
            maxh_m=args.maxh,
            outer_quad=args.outer_quad,
            arc_max_segment_length_m=args.arc_max_segment_length,
            exchange_file=args.export_model,
        )
        artifact = _make_success_artifact(
            details,
            args.output,
            args.export_model,
            time.perf_counter() - started,
        )
    except Exception as error:
        artifact = _make_failure_artifact(error, time.perf_counter() - started)
        args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
        raise

    write_started = time.perf_counter()
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    artifact["timing_breakdown_s"]["write_results"] = time.perf_counter() - write_started
    artifact["run"]["duration_s"] = time.perf_counter() - started
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(artifact, indent=2))
    if not artifact["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
