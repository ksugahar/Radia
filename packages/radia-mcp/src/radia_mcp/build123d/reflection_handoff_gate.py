"""Gates for reflection-sensitive STEP handoffs and proper-rotation recovery."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _positive(row: Mapping[str, object], name: str) -> float:
    value = float(row.get(name, float("nan")))
    if not isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _relative_error(actual: float, reference: float) -> float:
    return abs(actual - reference) / max(abs(reference), 1.0)


def build123d_reflection_rotation_handoff_gate(
    summary: Mapping[str, object],
    *,
    half_volume_rtol: float = 2.0e-6,
    rotation_compound_rtol: float = 2.0e-6,
    minimum_reflection_bias: float = 5.0e-3,
    minimum_fused_mirror_bias: float = 0.2,
    feature_delta_rtol: float = 1.0e-3,
) -> dict[str, object]:
    """Require failed reflection controls and a passing proper-rotation handoff."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    tolerances = (
        half_volume_rtol,
        rotation_compound_rtol,
        minimum_reflection_bias,
        minimum_fused_mirror_bias,
        feature_delta_rtol,
    )
    if any(not isfinite(float(value)) or float(value) <= 0.0 for value in tolerances):
        raise ValueError("all tolerances and minimum biases must be finite and positive")

    formula = _mapping(summary.get("volume_contract") or {}, "volume_contract")
    raw_rows = list(summary.get("checkpoints") or [])
    if not raw_rows or not all(isinstance(row, Mapping) for row in raw_rows):
        raise ValueError("checkpoints must contain mappings")
    rows = {str(row.get("name") or ""): row for row in raw_rows}
    required_names = (
        "half_before_fillet",
        "half_after_fillet",
        "full_after_mirror",
        "mirrored_half_alone",
        "rotated_half_alone",
        "two_body_rotation_compound",
    )
    if set(required_names) - set(rows):
        raise ValueError("all six named reflection/rotation checkpoints are required")
    if len(rows) != len(raw_rows):
        raise ValueError("checkpoint names must be unique")

    native: dict[str, float] = {}
    external: dict[str, float] = {}
    errors: dict[str, float] = {}
    for name in required_names:
        row = rows[name]
        native[name] = _positive(row, "build123d_volume_mm3")
        external[name] = _positive(row, "external_volume_mm3")
        errors[name] = _relative_error(external[name], native[name])

    unfilleted = _positive(formula, "unfilleted_volume_mm3")
    feature_delta = _positive(formula, "feature_delta_mm3")
    formula_final = _positive(formula, "formula_final_volume_mm3")
    build_final = _positive(formula, "build123d_final_volume_mm3")
    native_feature_delta = native["half_after_fillet"] - native["half_before_fillet"]
    external_feature_delta = external["half_after_fillet"] - external["half_before_fillet"]
    feature_delta_error = _relative_error(external_feature_delta, native_feature_delta)
    native_transform_spread = (
        max(
            native["half_after_fillet"],
            native["mirrored_half_alone"],
            native["rotated_half_alone"],
        )
        - min(
            native["half_after_fillet"],
            native["mirrored_half_alone"],
            native["rotated_half_alone"],
        )
    ) / native["half_after_fillet"]
    external_compound_additivity_error = _relative_error(
        external["two_body_rotation_compound"],
        external["half_after_fillet"] + external["rotated_half_alone"],
    )

    digest_ok = all(len(str(rows[name].get("step_sha256") or "")) == 64 for name in required_names)
    topology_ok = all(
        int(rows[name].get("body_count", -1)) == 1
        and int(rows[name].get("volume_count", -1)) == 1
        for name in required_names[:-1]
    ) and int(rows["two_body_rotation_compound"].get("body_count", -1)) == int(
        rows["two_body_rotation_compound"].get("volume_count", -1)
    ) == 2
    self_roundtrip_ok = all(
        _relative_error(
            _positive(rows[name], "step_roundtrip_volume_mm3"),
            native[name],
        )
        <= 1.0e-9
        for name in required_names
    )
    checks = {
        "base_plus_feature_formula_closes": _relative_error(
            unfilleted + feature_delta, formula_final
        )
        <= 1.0e-12
        and _relative_error(build_final, formula_final) <= 1.0e-10,
        "same_kernel_step_roundtrips_close": self_roundtrip_ok,
        "transform_native_volumes_are_invariant": native_transform_spread <= 1.0e-12,
        "checkpoint_digests_and_topology_are_bound": digest_ok and topology_ok,
        "positive_half_imports_are_healthy": errors["half_before_fillet"]
        <= float(half_volume_rtol)
        and errors["half_after_fillet"] <= float(half_volume_rtol),
        "fillet_delta_is_not_primary_translation_loss": feature_delta_error
        <= float(feature_delta_rtol),
        "reflection_negative_controls_expose_bias": errors["mirrored_half_alone"]
        >= float(minimum_reflection_bias)
        and errors["full_after_mirror"] >= float(minimum_fused_mirror_bias),
        "proper_rotation_half_restores_volume": errors["rotated_half_alone"]
        <= float(half_volume_rtol),
        "proper_rotation_compound_restores_volume": errors[
            "two_body_rotation_compound"
        ]
        <= float(rotation_compound_rtol),
        "proper_rotation_compound_is_additive": external_compound_additivity_error
        <= 1.0e-10,
        "solver_ready_disposition_is_explicit": summary.get("disposition")
        == "reflection_translation_bias_proper_rotation_compound_handoff_ready",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_reflection_rotation_handoff_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "external_volume_relative_errors": errors,
        "native_transform_volume_spread": native_transform_spread,
        "feature_delta_relative_error": feature_delta_error,
        "rotation_compound_additivity_error": external_compound_additivity_error,
        "recommended_handoff": "two_body_rotation_compound" if not issues else None,
        "notes": [
            "A same-kernel STEP roundtrip does not prove that a reflected BREP is portable to an independent CAD kernel.",
            "Use reflection artifacts as negative controls; for geometrically rotationally symmetric halves, a proper 180-degree rotation can preserve orientation and external mass properties.",
            "Keep the two bodies explicit at handoff and verify their external volume sum before meshing or solver-ready promotion.",
        ],
    }


def build123d_heat_exchanger_source_recovery_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Gate the upstream heat-exchanger replay and reflection recovery evidence."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    build = _mapping(summary.get("build") or {}, "build")
    checkpoints = _mapping(summary.get("checkpoint_run") or {}, "checkpoint_run")
    external = _mapping(summary.get("external_run") or {}, "external_run")
    process = _mapping(external.get("process") or {}, "external_run.process")
    handoff = _mapping(summary.get("handoff") or {}, "handoff")
    public_gate = build123d_reflection_rotation_handoff_gate(handoff)
    parameters = _mapping(build.get("parameters_mm") or {}, "build.parameters_mm")
    topology = _mapping(build.get("topology") or {}, "build.topology")
    build_checks = _mapping(build.get("checks") or {}, "build.checks")
    checkpoint_checks = _mapping(checkpoints.get("checks") or {}, "checkpoint_run.checks")
    process_exit = int(process.get("exit_code", -1))
    unexpected = list(process.get("unexpected_error_lines") or [])
    process_ok = process_exit == 0 or (
        process_exit in {1, 2, 3}
        and process.get("process_exit_policy")
        == "artifact_evidence_over_known_headless_diagnostics"
        and process.get("result_artifact_fresh") is True
        and not unexpected
        and process.get("acceptable") is True
    )
    source_sha_before = str(build.get("source_sha256_before", "")).lower()
    source_sha_after = str(build.get("source_sha256_after", "")).lower()
    build_timing = _mapping(build.get("timing") or {}, "build.timing")
    external_rows = list(external.get("rows") or [])
    checks = {
        "upstream_v0100_source_is_bound": build.get("source_kind")
        == "upstream_native_build123d_example_with_viewer_stub_only"
        and build.get("source_example") == "examples/heat_exchanger.py"
        and build.get("upstream_tag") == "v0.10.0"
        and len(str(build.get("upstream_git_blob_sha1") or "")) == 40,
        "source_sha256_is_preserved": len(source_sha_before) == 64
        and source_sha_before == source_sha_after,
        "runtime_topology_overrides_stale_comment": int(parameters.get("tube_count", 0))
        == int(parameters.get("tube_location_count", -1))
        == int(topology.get("runtime_tube_count", -2))
        == 148
        and int(topology.get("source_comment_tube_count", 0)) == 149
        and topology.get("source_comment_runtime_mismatch") is True,
        "upstream_formula_fillet_and_roundtrip_checks_passed": bool(build_checks)
        and all(value is True for value in build_checks.values()),
        "checkpoint_instrumentation_and_roundtrips_passed": checkpoints.get(
            "instrumentation"
        )
        == "runtime copy immediately before and after the upstream fillet call"
        and bool(checkpoint_checks)
        and all(value is True for value in checkpoint_checks.values()),
        "six_checkpoint_external_replay_completed": len(external_rows) == 6
        and {str(row.get("name") or "") for row in external_rows}
        == {
            "half_before_fillet",
            "half_after_fillet",
            "full_after_mirror",
            "mirrored_half_alone",
            "rotated_half_alone",
            "two_body_rotation_compound",
        },
        "headless_external_cad_process_is_classified": external.get("execution_mode")
        == "python_api_headless_synchronous_commands"
        and {"-nographics", "-batch"}.issubset(set(external.get("headless_flags") or []))
        and external.get("gui_daemon_enabled") is False
        and process_ok,
        "fresh_result_and_no_owned_process_leak": process.get("result_artifact_fresh")
        is True
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "exactly_four_build_timing_stages": len(build_timing) == 4
        and all(isfinite(float(value)) and float(value) >= 0.0 for value in build_timing.values()),
        "independent_reflection_rotation_gate_passed": public_gate["status"] == "ok",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_heat_exchanger_source_recovery_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "process_exit_code": process_exit,
        "runtime_tube_count": int(parameters.get("tube_count", 0)),
        "source_comment_tube_count": int(topology.get("source_comment_tube_count", 0)),
        "public_gate_status": public_gate["status"],
        "recommended_handoff": public_gate.get("recommended_handoff"),
        "notes": [
            "Derive array cardinality from the executed topology; a source comment may drift from current HexLocations behavior.",
            "Checkpoint immediately before and after the feature operation to localize translation loss instead of blaming fillets from the final model alone.",
            "A classified nonzero headless exit is evidence transport only; solver-ready status comes from the independent rotation-handoff gate.",
        ],
    }
