"""Drafted housing cross-kernel and source-replay validation gates."""

from __future__ import annotations

import math
from typing import Mapping


def _relative_error(measured: float, reference: float) -> float:
    if not math.isfinite(reference) or reference <= 0.0 or not math.isfinite(measured):
        return math.inf
    return abs(measured - reference) / reference


def _finite_nonnegative(name: str, value: float) -> float:
    result = float(value)
    if result < 0.0 or not math.isfinite(result):
        raise ValueError(f"{name} must be finite and nonnegative")
    return result


def drafted_housing_cross_kernel_gate(
    summary: Mapping[str, object],
    *,
    tessellated_rtol: float = 5.0e-4,
    step_volume_rtol: float = 5.0e-5,
    step_area_rtol: float = 1.0e-5,
    external_rtol: float = 1.0e-4,
    mesh_volume_rtol: float = 5.0e-3,
) -> dict[str, object]:
    """Gate a drafted, filleted, multiply-perforated solid across CAD kernels."""

    tolerances = {
        "tessellated": _finite_nonnegative("tessellated_rtol", tessellated_rtol),
        "step_volume": _finite_nonnegative("step_volume_rtol", step_volume_rtol),
        "step_area": _finite_nonnegative("step_area_rtol", step_area_rtol),
        "external": _finite_nonnegative("external_rtol", external_rtol),
        "mesh_volume": _finite_nonnegative("mesh_volume_rtol", mesh_volume_rtol),
    }
    runs = summary.get("runs") or []
    external_rows = summary.get("external_rows") or []
    mesh = summary.get("mesh_evidence") or {}
    gmsh = summary.get("gmsh_inventory") or {}
    if not isinstance(runs, list) or len(runs) != 2:
        raise ValueError("runs must contain exactly two native replays")
    if not isinstance(external_rows, list) or len(external_rows) != 4:
        raise ValueError("external_rows must contain two heal/noheal pairs")

    native_rows = []
    for run in runs:
        native = run.get("native") or {}
        tessellated = run.get("tessellated") or {}
        step = run.get("self_roundtrip") or {}
        native_volume = float(native.get("volume_mm3", math.nan))
        native_area = float(native.get("area_mm2", math.nan))
        native_rows.append(
            {
                "run_index": int(run.get("run_index", -1)),
                "native_valid": native.get("is_valid") is True,
                "native_solid_count": int(native.get("solid_count", -1)),
                "native_volume_mm3": native_volume,
                "native_area_mm2": native_area,
                "native_topology": [
                    int(native.get(key, -1))
                    for key in ("solid_count", "shell_count", "face_count", "edge_count", "vertex_count")
                ],
                "native_bbox_extent_mm": [float(value) for value in native.get("bbox_extent_mm", [])],
                "tessellated_volume_relative_error": _relative_error(float(tessellated.get("volume_mm3", math.nan)), native_volume),
                "tessellated_area_relative_error": _relative_error(float(tessellated.get("area_mm2", math.nan)), native_area),
                "step_valid": step.get("is_valid") is True,
                "step_solid_count": int(step.get("solid_count", -1)),
                "step_volume_relative_error": _relative_error(float(step.get("volume_mm3", math.nan)), native_volume),
                "step_area_relative_error": _relative_error(float(step.get("area_mm2", math.nan)), native_area),
            }
        )

    native_by_run = {row["run_index"]: row for row in native_rows}
    external = []
    for raw in external_rows:
        run_index = int(raw.get("run_index", -1))
        native = native_by_run.get(run_index, {})
        native_volume = float(native.get("native_volume_mm3", math.nan))
        native_area = float(native.get("native_area_mm2", math.nan))
        external.append(
            {
                "run_index": run_index,
                "import_mode": str(raw.get("import_mode", "")),
                "volume_count": int(raw.get("volume_count", -1)),
                "positive_volume_count": int(raw.get("positive_volume_count", -1)),
                "surface_count": int(raw.get("surface_count", -1)),
                "curve_count": int(raw.get("curve_count", -1)),
                "vertex_count": int(raw.get("vertex_count", -1)),
                "total_volume_mm3": float(raw.get("total_volume_mm3", math.nan)),
                "total_area_mm2": float(raw.get("total_area_mm2", math.nan)),
                "native_volume_relative_error": _relative_error(float(raw.get("total_volume_mm3", math.nan)), native_volume),
                "native_area_relative_error": _relative_error(float(raw.get("total_area_mm2", math.nan)), native_area),
            }
        )

    pairs = {
        run_index: sorted(
            [row for row in external if row["run_index"] == run_index],
            key=lambda row: row["import_mode"],
        )
        for run_index in native_by_run
    }
    native_reference = native_rows[0]
    mesh_tet_count = int(mesh.get("tet_count", -1))
    gmsh_counts = gmsh.get("element_family_counts") or {}
    checks = {
        "native_replay_indices_are_exact": sorted(native_by_run) == [1, 2],
        "native_runs_are_valid_single_solids": all(row["native_valid"] and row["native_solid_count"] == 1 and row["native_volume_mm3"] > 0.0 and row["native_area_mm2"] > 0.0 for row in native_rows),
        "native_mass_topology_bbox_replay_exact": all(_relative_error(row["native_volume_mm3"], native_reference["native_volume_mm3"]) <= 1.0e-12 and _relative_error(row["native_area_mm2"], native_reference["native_area_mm2"]) <= 1.0e-12 and row["native_topology"] == native_reference["native_topology"] and row["native_bbox_extent_mm"] == native_reference["native_bbox_extent_mm"] for row in native_rows),
        "tessellation_supports_native_mass": all(row["tessellated_volume_relative_error"] <= tolerances["tessellated"] and row["tessellated_area_relative_error"] <= tolerances["tessellated"] for row in native_rows),
        "same_kernel_step_handoff_is_bounded": all(row["step_valid"] and row["step_solid_count"] == 1 and row["step_volume_relative_error"] <= tolerances["step_volume"] and row["step_area_relative_error"] <= tolerances["step_area"] for row in native_rows),
        "external_modes_cover_both_replays": all({row["import_mode"] for row in pair} == {"heal", "noheal"} and len(pair) == 2 for pair in pairs.values()),
        "external_imports_are_positive_single_solids": all(row["volume_count"] == 1 and row["positive_volume_count"] == 1 for row in external),
        "external_mass_is_bounded_from_native": all(row["native_volume_relative_error"] <= tolerances["external"] and row["native_area_relative_error"] <= tolerances["external"] for row in external),
        "external_topology_is_invariant": len({(row["surface_count"], row["curve_count"], row["vertex_count"]) for row in external}) == 1,
        "heal_noheal_observables_are_invariant": all(len(pair) == 2 and _relative_error(pair[0]["total_volume_mm3"], pair[1]["total_volume_mm3"]) <= 1.0e-12 and _relative_error(pair[0]["total_area_mm2"], pair[1]["total_area_mm2"]) <= 1.0e-12 for pair in pairs.values()),
        "tet_mesh_is_positive_first_order": mesh_tet_count > 0 and list(mesh.get("connectivity_orders") or []) == [4] and float(mesh.get("min_scaled_jacobian", math.nan)) > 0.0,
        "tet_integrated_volume_closes_external_cad": float(mesh.get("cad_volume_relative_error", math.inf)) <= tolerances["mesh_volume"],
        "gmsh_v41_contains_the_complete_tet_block": gmsh.get("status") == "ok" and gmsh.get("mesh_format") == "4.1" and gmsh.get("binary") is False and int(gmsh.get("node_count", -1)) == int(mesh.get("node_count", -2)) and int(gmsh_counts.get("tet", -1)) == mesh_tet_count and not gmsh.get("connectivity_mismatches"),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_drafted_housing_cross_kernel_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "solver_ready": not issues,
        "checks": checks,
        "issues": issues,
        "native_runs": native_rows,
        "external_rows": external,
        "mesh": {
            "tet_count": mesh_tet_count,
            "node_count": int(mesh.get("node_count", -1)),
            "cad_volume_relative_error": float(mesh.get("cad_volume_relative_error", math.inf)),
            "min_scaled_jacobian": float(mesh.get("min_scaled_jacobian", math.nan)),
        },
        "notes": [
            "Use the native B-rep as the center of a bounded multi-kernel spread; do not nominate one STEP reader as exact truth.",
            "Draft, fillet, counterbore, and repeated through-hole combinations require mass and topology checks together.",
            "A Gmsh v4.1 file is complete only when its element block reproduces the measured tet inventory, not merely when nodes are present.",
        ],
    }

def drafted_housing_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate tagged source identity and the headless solver-handoff workflow."""

    source = summary.get("source") or {}
    contract = summary.get("source_contract") or {}
    execution = summary.get("external_execution") or {}
    companions = summary.get("gmsh_companions") or {}
    timing = summary.get("timing_breakdown_s") or {}
    operations = {str(value) for value in contract.get("operations", [])}
    parameters = contract.get("parameters") or {}
    checks = {
        "tagged_upstream_source_identity_recorded": str(source.get("kind", "")).startswith("upstream_native_example_tag_") and str(source.get("version", "")) == "0.10.0" and len(str(source.get("sha256", ""))) == 64 and len(str(source.get("commit", ""))) == 40,
        "source_preserved_except_display_stub": source.get("preserved") is True and source.get("display_stubbed_only") is True,
        "drafted_perforated_source_contract_recorded": {"draft", "fillet", "counterbore", "through_holes"}.issubset(operations) and float(parameters.get("draft_angle_deg", math.nan)) > 0.0 and int(parameters.get("mounting_hole_count", -1)) == 2 and float(parameters.get("counterbore_depth_mm", math.nan)) > 0.0,
        "two_native_replays_and_two_step_artifacts_recorded": int(summary.get("native_replay_count", -1)) == 2 and len(summary.get("step_sha256") or []) == 2 and all(len(str(value)) == 64 for value in summary.get("step_sha256") or []),
        "headless_external_cad_execution": execution.get("mode") == "python_api_headless" and {"-nographics", "-batch"}.issubset(set(execution.get("headless_flags") or [])) and execution.get("gui_daemon_enabled") is False,
        "fresh_result_and_owned_process_cleanup": execution.get("result_artifact_fresh") is True and int(execution.get("owned_processes_remaining", -1)) == 0,
        "gmsh_mesh_and_launch_companions_recorded": companions.get("msh") is True and companions.get("geo") is True and companions.get("geo_opt") is True and companions.get("msh_opt") is True,
        "cross_kernel_gate_accepts_solver_handoff": summary.get("cross_kernel_gate_status") == "ok" and summary.get("solver_ready") is True,
        "four_dominant_timing_stages_recorded": len(timing) == 4 and all(float(value) >= 0.0 for value in timing.values()),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_drafted_housing_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "operations": sorted(operations),
        "notes": [
            "Run the unmodified tagged example with only its interactive display call stubbed for headless execution.",
            "Persist both STEP digests because equivalent replay metrics do not imply byte-identical exports.",
            "Require .geo/.geo.opt launch companions beside the Gmsh v4.1 mesh so the artifact opens with the intended view contract.",
        ],
    }
