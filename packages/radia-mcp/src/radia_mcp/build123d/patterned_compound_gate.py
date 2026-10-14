"""Patterned compound and source-native curved-surface validation gates."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping


_STARTUP_SUFFIXES = ("/plugins", "-commandplugindir", "-nojournal")


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def patterned_compound_translation_gate(
    summary: Mapping[str, object],
    *,
    same_kernel_volume_rtol: float = 1.0e-6,
    external_solver_volume_rtol: float = 1.0e-5,
    max_patterned_member_volume_rtol: float = 1.0e-3,
    min_bias_localization_ratio: float = 10.0,
) -> dict[str, object]:
    """Diagnose an external STEP bias without treating topology as mass closure.

    ``status=ok`` means that the translator-bias diagnosis is complete.  The
    independent ``solver_ready`` field remains false while the external volume
    errors exceed the solver handoff tolerance.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    tolerances = (
        float(same_kernel_volume_rtol),
        float(external_solver_volume_rtol),
        float(max_patterned_member_volume_rtol),
        float(min_bias_localization_ratio),
    )
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances and localization ratio must be finite and nonnegative")
    if external_solver_volume_rtol == 0.0 or min_bias_localization_ratio == 0.0:
        raise ValueError("external tolerance and localization ratio must be positive")

    authoring = _mapping(summary.get("authoring"), "authoring")
    step = _mapping(authoring.get("step"), "authoring.step")
    external = _mapping(summary.get("external"), "external")
    raw_rows = external.get("rows")
    if not isinstance(raw_rows, list) or len(raw_rows) != 2:
        raise ValueError("external.rows must contain exactly heal and noheal rows")
    if not all(isinstance(row, Mapping) for row in raw_rows):
        raise ValueError("external.rows must contain mappings")
    rows = list(raw_rows)
    expected_count = int(authoring.get("authoring_solid_count", 0))
    topology_keys = ("body_count", "volume_count", "surface_count", "curve_count", "vertex_count")
    modes = [str(row.get("mode", "")) for row in rows]
    total_errors = [float(row.get("total_volume_relative_error", math.inf)) for row in rows]
    dominant_errors = [float(row.get("tire_volume_relative_error", math.inf)) for row in rows]
    patterned_errors = [
        float(row.get("tread_total_volume_relative_error", math.inf)) for row in rows
    ]
    localization_ratios = [
        dominant / max(patterned, 1.0e-30)
        for dominant, patterned in zip(dominant_errors, patterned_errors)
    ]
    same_kernel_error = float(step.get("roundtrip_volume_relative_error", math.inf))
    external_solver_ready = all(
        total <= external_solver_volume_rtol
        and float(row.get("maximum_sorted_body_volume_relative_error", math.inf))
        <= external_solver_volume_rtol
        for total, row in zip(total_errors, rows)
    )
    checks = {
        "authoring_and_same_kernel_roundtrip_passed": authoring.get("pass") is True
        and int(step.get("roundtrip_solid_count", -1)) == expected_count
        and same_kernel_error <= same_kernel_volume_rtol,
        "same_step_digest_reached_external_cad": len(str(step.get("sha256", ""))) == 64
        and str(step.get("sha256", "")) == str(external.get("step_sha256", "")),
        "heal_and_noheal_rows_complete": modes == ["noheal", "heal"],
        "all_patterned_bodies_preserved": expected_count > 1
        and all(
            int(row.get("body_count", -1))
            == int(row.get("volume_count", -1))
            == expected_count
            for row in rows
        ),
        "heal_noheal_topology_identical": all(
            rows[0].get(key) == rows[1].get(key) for key in topology_keys
        ),
        "heal_noheal_bias_identical": abs(total_errors[0] - total_errors[1]) <= 1.0e-12,
        "external_translation_bias_is_material": all(
            error > external_solver_volume_rtol for error in total_errors
        ),
        "bias_is_localized_to_dominant_curved_body": all(
            patterned <= max_patterned_member_volume_rtol
            and ratio >= min_bias_localization_ratio
            for patterned, ratio in zip(patterned_errors, localization_ratios)
        ),
        "external_solver_handoff_is_rejected": external.get("external_solver_ready") is False
        and external_solver_ready is False
        and external.get("disposition") == "translation_bias_detected_not_solver_ready",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_patterned_compound_translation_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "diagnosis": "dominant_curved_body_translation_bias" if not issues else "incomplete_evidence",
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "expected_body_count": expected_count,
        "same_kernel_volume_relative_error": same_kernel_error,
        "external_total_volume_relative_errors": total_errors,
        "dominant_body_volume_relative_errors": dominant_errors,
        "patterned_member_total_volume_relative_errors": patterned_errors,
        "bias_localization_ratios": localization_ratios,
        "external_solver_volume_tolerance": float(external_solver_volume_rtol),
        "notes": [
            "Preserved body count and heal/noheal agreement do not prove external mass-property closure.",
            "Compare the dominant curved body separately from repeated surface-wrapped members to localize translator bias.",
            "An evidence gate may pass while solver_ready remains false; never turn a diagnosed bias into a loose acceptance tolerance.",
        ],
    }


def wrap_faces_rotational_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate immutable upstream wrap/thicken/rotation behavior and replay evidence."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    authoring = _mapping(summary.get("authoring"), "authoring")
    external = _mapping(summary.get("external"), "external")
    process = _mapping(summary.get("external_process"), "external_process")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    groups = authoring.get("prototype_groups")
    if not isinstance(groups, list) or not all(isinstance(row, Mapping) for row in groups):
        raise ValueError("authoring.prototype_groups must contain mappings")
    public_gate = patterned_compound_translation_gate(summary)
    required_operations = {"Bezier", "revolve", "wrap_faces", "thicken", "Rot"}
    operations = {str(value) for value in authoring.get("source_operations") or []}
    diagnostics = [str(value) for value in process.get("startup_diagnostics") or []]
    normalized_suffixes = {
        suffix
        for row in diagnostics
        for suffix in _STARTUP_SUFFIXES
        if row.rstrip().endswith(suffix)
    }
    startup_diagnostics_only = (
        len(diagnostics) == len(_STARTUP_SUFFIXES)
        and normalized_suffixes == set(_STARTUP_SUFFIXES)
        and all("Could not open file:" in row for row in diagnostics)
        and not list(process.get("script_error_lines") or [])
    )
    source_before = str(authoring.get("source_sha256_before", "")).lower()
    source_after = str(authoring.get("source_sha256_after", "")).lower()
    timing_values = list(timing.values())
    checks = {
        "upstream_source_identity_and_digest_recorded": authoring.get("upstream_tag") == "v0.10.0"
        and len(str(authoring.get("upstream_commit", ""))) == 40
        and Path(str(authoring.get("source_example", ""))).as_posix()
        == "examples/bicycle_tire.py"
        and len(source_before) == 64,
        "source_preserved_with_viewer_stub_only": source_before == source_after
        and authoring.get("viewer_stubbed") is True
        and authoring.get("source_kind")
        == "upstream_native_build123d_example_with_viewer_stub_only",
        "wrap_thicken_rotation_sequence_recorded": required_operations == operations,
        "six_by_180_pattern_identity": int(authoring.get("prototype_count", 0)) == 6
        and len(groups) == 6
        and int(authoring.get("rotation_count_per_prototype", 0)) == 180
        and int(authoring.get("rotation_step_degrees", 0)) == 2
        and int(authoring.get("tread_solid_count", 0)) == 1080,
        "every_rotational_copy_is_valid_and_invariant": all(
            int(row.get("copy_count", 0)) == 180
            and int(row.get("valid_copy_count", -1)) == 180
            and float(row.get("copy_to_prototype_max_relative_error", math.inf)) <= 1.0e-12
            and float(row.get("radial_center_relative_spread", math.inf)) <= 1.0e-12
            and float(row.get("x_center_spread", math.inf)) <= 1.0e-12
            for row in groups
        ),
        "headless_external_replay_recorded": external.get("execution_mode")
        == "python_api_headless_synchronous_commands"
        and {"-nographics", "-batch"}.issubset(set(external.get("headless_flags") or []))
        and external.get("gui_daemon_enabled") is False,
        "startup_only_nonzero_exit_explained": int(process.get("exit_code", -1)) in {0, 2, 3}
        and (int(process.get("exit_code", -1)) == 0 or startup_diagnostics_only),
        "fresh_artifact_and_no_owned_process_leak": process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "four_dominant_timing_stages_recorded": len(timing_values) == 4
        and all(math.isfinite(float(value)) and float(value) >= 0.0 for value in timing_values),
        "translation_diagnosis_passed_without_solver_ready_overclaim": public_gate["status"] == "ok"
        and public_gate["solver_ready"] is False,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_wrap_faces_rotational_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "prototype_count": len(groups),
        "copy_count": sum(int(row.get("copy_count", 0)) for row in groups),
        "public_gate_status": public_gate["status"],
        "solver_ready": False,
        "notes": [
            "Stub only the viewer when replaying an upstream example; geometry statements and source digest must remain unchanged.",
            "For Rot copies, verify volume, radial center, axial center, validity, prototype count, and angular coverage.",
            "A known batch startup diagnostic can explain a nonzero launcher code only with no script errors, a fresh result, and zero leaked owned processes.",
        ],
    }
