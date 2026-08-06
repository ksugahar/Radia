"""Build the moving TEAM 28 HCurl eddy-bubble family from CoilBuilder fields.

One p=6 HCurl parent mesh and one EVRS current basis are reused at every
height. Only the translated incident CoilBuilder A/B fields, port projection,
and force operator vary. The exported family can therefore be interpolated
without a hidden state transfer and can drive an averaged mechanical plant.
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

from team28_coilbuilder_eddy_bubble import (
    REPO_ROOT,
    _evaluate_force,
    _force_operator,
    _repo_path,
    compare_coil_fields,
)
from team28_hcurl_vim_force import (
    FREQUENCY_HZ,
    REFERENCE_COIL_CURRENT_A,
    SIGMA_AL,
    _build_eddy_basis,
)

HERE = Path(__file__).resolve().parent
DEFAULT_RESULT = HERE / "team28_coilbuilder_height_family_results.json"
DEFAULT_FAMILY = HERE / "team28_coilbuilder_hcurl_eddy_cln_family.json"
REFERENCE_CURVE = (
    REPO_ROOT / "docs" / "maglev" / "demos" / "team28" / "team28_cln_sweep_results.json"
)
DEFAULT_HEIGHT_OFFSETS_MM = tuple(range(-7, 18))
DISK_WEIGHT_N = 1.055


def _load_reference_curve(height_offsets_mm):
    payload = json.loads(REFERENCE_CURVE.read_text(encoding="utf-8"))
    stored = {
        float(height): -0.5 * float(force)
        for height, force in zip(payload["dZ_mm"], payload["fz_cln_N"])
    }
    return np.asarray([stored[float(value)] for value in height_offsets_mm])


def _equilibrium_offset_m(height_m, upward_force_N, weight_N=DISK_WEIGHT_N):
    residual = np.asarray(upward_force_N) - float(weight_N)
    for index in range(len(height_m) - 1):
        left = residual[index]
        right = residual[index + 1]
        if left == 0.0:
            return float(height_m[index])
        if left * right <= 0.0:
            fraction = -left / (right - left)
            return float(height_m[index] + fraction * (height_m[index + 1] - height_m[index]))
    raise RuntimeError("the height family does not bracket force equals weight")


def run_family(
    height_offsets_mm=DEFAULT_HEIGHT_OFFSETS_MM,
    *,
    maxh_m=0.025,
    outer_quad=4,
    arc_max_segment_length_m=0.002,
    family_file=DEFAULT_FAMILY,
):
    offsets_mm = np.asarray(height_offsets_mm, dtype=float)
    if offsets_mm.ndim != 1 or len(offsets_mm) < 3:
        raise ValueError("height_offsets_mm must contain at least three values")
    if np.any(~np.isfinite(offsets_mm)) or np.any(np.diff(offsets_mm) <= 0.0):
        raise ValueError("height_offsets_mm must be finite and strictly increasing")
    offsets_m = offsets_mm * 1.0e-3

    timings = {}
    stage = time.perf_counter()
    ng, mesh, fes, basis = _build_eddy_basis(maxh_m)
    timings["mesh_and_common_eddy_basis"] = time.perf_counter() - stage

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
    timings["shared_vim_assembly"] = time.perf_counter() - stage

    stage = time.perf_counter()
    snapshots = []
    rows = []
    field_a_errors = []
    field_b_errors = []
    cross_check_indices = np.linspace(
        0,
        basis.current_basis.n_samples - 1,
        1024,
        dtype=int,
    )
    s = 2.0j * np.pi * FREQUENCY_HZ
    for offset_mm, offset_m in zip(offsets_mm, offsets_m):
        external_a, external_b, source_report = compare_coil_fields(
            basis.current_basis.points,
            coil_current_A=REFERENCE_COIL_CURRENT_A,
            arc_max_segment_length_m=arc_max_segment_length_m,
            height_offset_m=offset_m,
            cross_check_indices=cross_check_indices,
            reference_radial_order=24,
            reference_axial_order=32,
        )
        field_check = source_report["field_cross_check"]
        field_a_errors.append(field_check["vector_potential_relative_l2"])
        field_b_errors.append(field_check["flux_density_relative_l2"])
        unit_external_a = external_a / REFERENCE_COIL_CURRENT_A
        unit_external_b = external_b / REFERENCE_COIL_CURRENT_A
        rhs = vim.ExternalVectorPotentialRHS(
            basis.current_basis,
            unit_external_a,
        )
        model = vim.HCurlEddyCLNFromVIM(system, rhs)
        force_operator = _force_operator(basis.current_basis, unit_external_b)
        coefficients = model.solve_vector_potential_drive(
            s,
            REFERENCE_COIL_CURRENT_A,
        )
        force_N = _evaluate_force(
            force_operator,
            coefficients,
            REFERENCE_COIL_CURRENT_A,
        )
        snapshots.append(
            {
                "height_m": float(offset_m),
                "model": model,
                "force_operator": force_operator,
                "metadata": {
                    "height_offset_m": float(offset_m),
                    "height_offset_mm": float(offset_mm),
                    "coil_source": "radia.coil_builder.CoilBuilder",
                },
            }
        )
        rows.append(
            {
                "height_offset_mm": float(offset_mm),
                "height_offset_m": float(offset_m),
                "absolute_disk_bottom_mm": float(10.8 + offset_mm),
                "force_N": force_N.tolist(),
                "upward_force_N": float(force_N[2]),
                "vector_potential_relative_l2": float(field_a_errors[-1]),
                "flux_density_relative_l2": float(field_b_errors[-1]),
            }
        )
    timings["coilbuilder_height_sweep_and_solves"] = time.perf_counter() - stage

    stage = time.perf_counter()
    vim.ExportHCurlEddyCLNFamilyJSON(
        snapshots,
        family_file,
        metadata={
            "benchmark": "TEAM 28 electrodynamic levitation",
            "frequency_hz": FREQUENCY_HZ,
            "coil_current_reference_A": REFERENCE_COIL_CURRENT_A,
            "conductivity_S_per_m": SIGMA_AL,
            "coil_source": "radia.coil_builder.CoilBuilder",
            "coil_count": 2,
            "winding_directions": "counter-wound",
            "parent_space": "HCurl",
            "parent_order": 6,
            "height_coordinate": "offset from 10.8 mm disk-bottom reference",
            "state_basis": "common p=6 local disk EVRS basis",
            "force_convention": "positive z is upward physical time-average force",
        },
    )
    timings["family_export"] = time.perf_counter() - stage

    predicted_lift = np.asarray([row["upward_force_N"] for row in rows])
    reference_lift = _load_reference_curve(offsets_mm)
    curve_difference = predicted_lift - reference_lift
    curve_scale = max(float(np.max(np.abs(reference_lift))), 1.0e-30)
    normalized_max_abs_error = float(np.max(np.abs(curve_difference)) / curve_scale)
    significant = np.abs(reference_lift) > 0.05
    significant_max_relative_error = float(
        np.max(np.abs(curve_difference[significant]) / np.abs(reference_lift[significant]))
    )
    predicted_equilibrium_m = _equilibrium_offset_m(offsets_m, predicted_lift)
    reference_equilibrium_m = _equilibrium_offset_m(offsets_m, reference_lift)
    equilibrium_error_m = abs(predicted_equilibrium_m - reference_equilibrium_m)
    interaction_report = interaction.diagnostics()
    checks = {
        "twenty_five_strictly_ordered_height_snapshots": bool(
            len(rows) == 25 and np.all(np.diff(offsets_m) > 0.0)
        ),
        "all_snapshots_share_rank_three_state_basis": bool(
            basis.rank == 3 and all(snapshot["model"].state_order == 3 for snapshot in snapshots)
        ),
        "all_cln_snapshots_passive": bool(
            all(snapshot["model"].diagnostics()["passive"] for snapshot in snapshots)
        ),
        "all_coilbuilder_A_checks_below_0p1_percent": bool(max(field_a_errors) < 1.0e-3),
        "all_coilbuilder_B_checks_below_0p5_percent": bool(max(field_b_errors) < 5.0e-3),
        "lift_curve_global_error_below_two_percent": bool(normalized_max_abs_error < 0.02),
        "lift_curve_is_monotonically_nonincreasing": bool(
            np.all(np.diff(predicted_lift) <= 1.0e-9)
        ),
        "force_equals_weight_is_bracketed": bool(
            predicted_lift[0] > DISK_WEIGHT_N > predicted_lift[-1]
        ),
        "equilibrium_offset_matches_below_0p5_mm": bool(equilibrium_error_m < 0.5e-3),
        "interaction_projection_residual_below_1e-10": bool(
            interaction_report["projection_relative_residual"] < 1.0e-10
        ),
        "kernel_is_epsilon_free": bool(interaction_report["kernel_epsilon_m"] is None),
        "family_exchange_written": bool(Path(family_file).is_file()),
    }
    return {
        "timing_breakdown_s": timings,
        "problem": {
            "frequency_hz": FREQUENCY_HZ,
            "coil_current_reference_A": REFERENCE_COIL_CURRENT_A,
            "disk_weight_N": DISK_WEIGHT_N,
            "reference_disk_bottom_mm": 10.8,
            "mesh_maxh_m": float(maxh_m),
            "parent_space": "HCurl",
            "parent_order": 6,
            "parent_ndof": int(fes.ndof),
            "evrs_rank": int(basis.rank),
            "mesh_elements": int(mesh.ne),
        },
        "height_family": rows,
        "curve_comparison": {
            "reference": "stored TEAM 28 physical force-height regression",
            "reference_upward_force_N": reference_lift.tolist(),
            "predicted_upward_force_N": predicted_lift.tolist(),
            "normalized_max_absolute_error": normalized_max_abs_error,
            "significant_point_max_relative_error": significant_max_relative_error,
            "predicted_equilibrium_offset_mm": float(1.0e3 * predicted_equilibrium_m),
            "reference_equilibrium_offset_mm": float(1.0e3 * reference_equilibrium_m),
            "equilibrium_offset_error_mm": float(1.0e3 * equilibrium_error_m),
            "predicted_absolute_disk_bottom_mm": float(10.8 + 1.0e3 * predicted_equilibrium_m),
        },
        "field_cross_check": {
            "maximum_vector_potential_relative_l2": float(max(field_a_errors)),
            "maximum_flux_density_relative_l2": float(max(field_b_errors)),
        },
        "interaction": interaction_report,
        "checks": checks,
    }


def make_artifact(details, result_file, family_file, duration_s):
    checks = {
        "ran_to_completion": True,
        "result_files_exist": bool(Path(family_file).is_file()),
        "validation_passed": bool(all(details["checks"].values())),
        **details["checks"],
    }
    return {
        "radia_version": getattr(radia, "__version__", "unknown"),
        "schema": "cae-ai-lab.solver-run.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "case": "TEAM 28 CoilBuilder moving HCurl eddy-bubble family",
        "solver": "radia-ngsolve",
        "source_artifact": _repo_path(__file__),
        "pass": bool(all(checks.values())),
        "run": {
            "command": ("python validation_test/maglev/" "team28_coilbuilder_height_family.py"),
            "workdir": ".",
            "exit_code": 0,
            "duration_s": float(duration_s),
        },
        "result_files": [_repo_path(result_file), _repo_path(family_file)],
        "checks": checks,
        "tolerances": {
            "coil_vector_potential_max_rel": 1.0e-3,
            "coil_flux_density_max_rel": 5.0e-3,
            "lift_curve_global_max_rel": 0.02,
            "equilibrium_offset_max_abs_m": 0.5e-3,
        },
        "errors": {
            "max_rel": details["curve_comparison"]["normalized_max_absolute_error"],
            "max_abs": details["curve_comparison"]["equilibrium_offset_error_mm"] * 1.0e-3,
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
                "common-basis identity, source-field quadrature, stored force-height "
                "regression, monotonicity, passivity, and force-weight equilibrium"
            ),
            "command": (
                "pytest -q validation_test/maglev/"
                "test_team28_coilbuilder_height_family_evidence.py"
            ),
        },
        "details": details,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--height-offsets-mm",
        type=float,
        nargs="+",
        default=DEFAULT_HEIGHT_OFFSETS_MM,
    )
    parser.add_argument("--maxh", type=float, default=0.025)
    parser.add_argument("--outer-quad", type=int, default=4)
    parser.add_argument("--arc-max-segment-length", type=float, default=0.002)
    parser.add_argument("--output", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--export-family", type=Path, default=DEFAULT_FAMILY)
    args = parser.parse_args()

    started = time.perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.export_family.parent.mkdir(parents=True, exist_ok=True)
    details = run_family(
        args.height_offsets_mm,
        maxh_m=args.maxh,
        outer_quad=args.outer_quad,
        arc_max_segment_length_m=args.arc_max_segment_length,
        family_file=args.export_family,
    )
    artifact = make_artifact(
        details,
        args.output,
        args.export_family,
        time.perf_counter() - started,
    )
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
