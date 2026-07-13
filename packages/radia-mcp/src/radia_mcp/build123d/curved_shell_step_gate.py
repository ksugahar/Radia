"""Diagnostic gates for curved shell/sweep STEP semantics across CAD kernels."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _mapping(value: object, field: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be a mapping")
    return value


def _rows(value: object, field: str) -> list[Mapping[str, object]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError(f"{field} must be a sequence of mappings")
    rows = list(value)
    if not all(isinstance(row, Mapping) for row in rows):
        raise ValueError(f"{field} must contain mappings")
    return rows


def _finite(value: object, field: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{field} must be a finite number")
    return result


def _topology(row: Mapping[str, object], *, external: bool = False) -> tuple[int, int, int, int]:
    prefix = "" if external else ""
    return (
        int(row.get(f"{prefix}solid_count", row.get("volume_count", -1))),
        int(row.get(f"{prefix}face_count", row.get("surface_count", -1))),
        int(row.get(f"{prefix}edge_count", row.get("curve_count", -1))),
        int(row.get(f"{prefix}vertex_count", -1)),
    )


def build123d_curved_shell_step_semantics_gate(
    summary: Mapping[str, object],
    *,
    max_same_kernel_volume_relative_error: float = 1.0e-5,
    max_same_kernel_area_relative_error: float = 1.0e-6,
    min_external_volume_relative_error: float = 1.0e-2,
    min_external_area_relative_error: float = 1.0e-3,
    min_external_to_same_kernel_ratio: float = 1000.0,
) -> dict[str, object]:
    """Classify topology-preserving mass loss without solver-ready overclaim.

    A STEP roundtrip can retain one solid and identical face/edge/vertex counts
    while changing curved surface semantics.  This gate separates a small
    same-kernel drift from a much larger independent-kernel drift and requires
    heal/noheal invariance before diagnosing a semantic portability failure.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    thresholds = [
        _finite(max_same_kernel_volume_relative_error, "max_same_kernel_volume_relative_error"),
        _finite(max_same_kernel_area_relative_error, "max_same_kernel_area_relative_error"),
        _finite(min_external_volume_relative_error, "min_external_volume_relative_error"),
        _finite(min_external_area_relative_error, "min_external_area_relative_error"),
        _finite(min_external_to_same_kernel_ratio, "min_external_to_same_kernel_ratio"),
    ]
    if any(value <= 0.0 for value in thresholds):
        raise ValueError("all thresholds must be positive")

    build = _mapping(summary.get("build"), "build")
    native = _mapping(build.get("native"), "build.native")
    reimport = _mapping(build.get("step_reimport"), "build.step_reimport")
    same_errors = _mapping(build.get("errors"), "build.errors")
    external = _mapping(summary.get("external"), "external")
    imports = _rows(external.get("imports"), "external.imports")
    if len(imports) < 2:
        raise ValueError("external.imports must contain at least noheal and heal rows")

    same_volume_error = _finite(same_errors.get("volume_relative"), "build.errors.volume_relative")
    same_area_error = _finite(
        same_errors.get("surface_area_relative"), "build.errors.surface_area_relative"
    )
    external_volume_errors = [
        _finite(row.get("native_volume_relative_error"), f"external.imports[{index}].native_volume_relative_error")
        for index, row in enumerate(imports)
    ]
    external_area_errors = [
        _finite(row.get("native_area_relative_error"), f"external.imports[{index}].native_area_relative_error")
        for index, row in enumerate(imports)
    ]
    mode_volume_spread = _finite(
        external.get("mode_volume_relative_spread"), "external.mode_volume_relative_spread"
    )
    mode_area_spread = _finite(
        external.get("mode_area_relative_spread"), "external.mode_area_relative_spread"
    )
    native_topology = _topology(native)
    reimport_topology = _topology(reimport)
    external_topologies = [_topology(row, external=True) for row in imports]
    modes = {str(row.get("mode", "")) for row in imports}
    export_kinds = {
        str(row.get("export_kind", ""))
        for row in imports
        if row.get("export_kind") is not None
    }
    export_mode_pairs = {
        (str(row.get("export_kind", "")), str(row.get("mode", "")))
        for row in imports
        if row.get("export_kind") is not None
    }
    export_variants_complete = not export_kinds or (
        export_kinds == {"compound", "solid"}
        and export_mode_pairs
        == {
            ("compound", "noheal"),
            ("compound", "heal"),
            ("solid", "noheal"),
            ("solid", "heal"),
        }
    )
    ratio_volume = min(external_volume_errors) / max(same_volume_error, 1.0e-300)
    ratio_area = min(external_area_errors) / max(same_area_error, 1.0e-300)

    checks = {
        "native_valid_single_solid": native.get("is_valid") is True and native_topology[0] == 1,
        "same_kernel_topology_preserved": native_topology == reimport_topology
        and min(native_topology) > 0,
        "same_kernel_drift_small_but_nonzero": 0.0 < same_volume_error
        <= float(max_same_kernel_volume_relative_error)
        and 0.0 < same_area_error <= float(max_same_kernel_area_relative_error),
        "external_modes_cover_noheal_and_heal": modes == {"noheal", "heal"},
        "external_export_variants_are_complete": export_variants_complete,
        "external_topology_still_matches_native": all(
            topology == native_topology for topology in external_topologies
        ),
        "external_mass_loss_is_material": all(
            error >= float(min_external_volume_relative_error) for error in external_volume_errors
        )
        and all(error >= float(min_external_area_relative_error) for error in external_area_errors),
        "external_loss_dominates_same_kernel_drift": ratio_volume
        >= float(min_external_to_same_kernel_ratio)
        and ratio_area >= float(min_external_to_same_kernel_ratio),
        "healing_does_not_change_external_result": mode_volume_spread <= 1.0e-12
        and mode_area_spread <= 1.0e-12,
        "source_and_external_diagnostics_completed": build.get("pass") is True
        and external.get("pass") is True,
        "non_solver_ready_disposition_recorded": external.get("disposition")
        == "external_kernel_mass_loss_with_topology_preserved_not_solver_ready",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_curved_shell_step_semantics_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "solver_ready": False,
        "diagnosis": "topology_preserved_external_curved_surface_mass_loss" if not issues else "incomplete_diagnosis",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "same_kernel_volume_relative_error": same_volume_error,
            "same_kernel_area_relative_error": same_area_error,
            "external_volume_relative_errors": external_volume_errors,
            "external_area_relative_errors": external_area_errors,
            "external_to_same_kernel_volume_ratio": ratio_volume,
            "external_to_same_kernel_area_ratio": ratio_area,
            "native_topology": list(native_topology),
            "external_topologies": [list(row) for row in external_topologies],
            "external_modes": sorted(modes),
            "external_export_kinds": sorted(export_kinds),
        },
        "notes": [
            "Equal solid/face/edge/vertex counts do not prove equal curved geometry or mass properties.",
            "Heal/noheal invariance rules out a simple healing toggle as the recovery path.",
            "Retain the STEP as a negative portability control; do not route it to a solver until an independent mass gate passes.",
        ],
    }


def build123d_vase_external_solid_contract_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Bind an immutable upstream vase replay to a four-path rejection gate."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    build = _mapping(summary.get("build"), "build")
    external = _mapping(summary.get("external"), "external")
    source = _mapping(build.get("source_native_example"), "build.source_native_example")
    native = _mapping(build.get("native"), "build.native")
    checks_source = _mapping(build.get("checks"), "build.checks")
    process = _mapping(external.get("process"), "external.process")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    imports = _rows(external.get("imports"), "external.imports")
    source_sha = str(source.get("source_sha256", "")).lower()
    rows_by_pair = {
        (str(row.get("export_kind", "")), str(row.get("mode", ""))): row
        for row in imports
    }
    required_pairs = {
        ("compound", "noheal"),
        ("compound", "heal"),
        ("solid", "noheal"),
        ("solid", "heal"),
    }
    native_volume = _finite(native.get("volume"), "build.native.volume")
    native_area = _finite(native.get("area"), "build.native.area")
    public_gate = build123d_curved_shell_step_semantics_gate(summary)

    checks = {
        "upstream_source_identity_bound": summary.get("source_kind")
        == "upstream-tagged-example-exact-replay"
        and source.get("repository") == "gumyr/build123d"
        and source.get("tag") == "v0.10.0"
        and source.get("path") == "examples/vase.py"
        and len(source_sha) == 64
        and all(character in "0123456789abcdef" for character in source_sha),
        "viewer_bridge_only": source.get("viewer_bridge")
        == "no-op ocp_vscode show/show_object module; source text unchanged",
        "native_replay_is_deterministic_and_valid": checks_source.get(
            "source_runs_are_deterministic_with_kernel_tolerance"
        )
        is True
        and checks_source.get("native_shape_valid") is True
        and checks_source.get("single_solid_shell") is True
        and native_volume > 0.0
        and native_area > 0.0,
        "same_kernel_step_and_brep_close": all(
            checks_source.get(name) is True
            for name in (
                "compound_step_shape_valid",
                "solid_step_shape_valid",
                "brep_shape_valid",
                "compound_step_volume_matches_native",
                "solid_step_volume_matches_native",
                "brep_volume_matches_native",
                "compound_step_area_matches_native",
                "solid_step_area_matches_native",
                "brep_area_matches_native",
                "compound_step_topology_matches_native",
                "solid_step_topology_matches_native",
                "brep_topology_matches_native",
            )
        ),
        "four_external_import_paths_complete": len(imports) == 4
        and set(rows_by_pair) == required_pairs,
        "external_entity_counts_repeat_but_measure_is_zero": all(
            int(row.get("volume_count", -1)) == 1
            and int(row.get("surface_count", -1)) == int(native.get("face_count", -2))
            and int(row.get("curve_count", -1)) == int(native.get("edge_count", -2))
            and int(row.get("vertex_count", -1)) == int(native.get("vertex_count", -2))
            and _finite(row.get("api_volume"), "external.imports.api_volume") == 0.0
            and _finite(row.get("wrapper_volume"), "external.imports.wrapper_volume")
            == 0.0
            for row in imports
        ),
        "self_roundtrip_positive_for_every_external_source": all(
            row.get("self_roundtrip_valid") is True
            and _finite(
                row.get("self_roundtrip_volume"),
                "external.imports.self_roundtrip_volume",
            )
            > 0.0
            for row in imports
        ),
        "heal_and_solid_extraction_do_not_recover_measure": all(
            _finite(row.get("native_volume_relative_error"), "native_volume_relative_error")
            == 1.0
            for row in imports
        ),
        "invalid_external_solid_is_not_meshed_or_promoted": external.get(
            "external_cad_volume_valid"
        )
        is False
        and external.get("mesh_attempted") is False
        and external.get("solver_ready") is False,
        "external_replay_is_headless_and_synchronous": external.get("execution_mode")
        == "combined_journal_headless"
        and {"-nographics", "-batch"}.issubset(
            set(external.get("headless_flags") or [])
        )
        and external.get("gui_daemon_enabled") is False,
        "external_process_is_fresh_and_classified": process.get("acceptable") is True
        and process.get("result_artifact_fresh") is True
        and not list(process.get("unexpected_error_lines") or [])
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "exactly_four_timing_stages_recorded": len(timing) == 4
        and all(
            _finite(value, f"timing_breakdown_s.{name}") >= 0.0
            for name, value in timing.items()
        ),
        "public_semantics_gate_passed": public_gate["status"] == "ok",
        "negative_portability_outcome_not_overclaimed": public_gate["solver_ready"]
        is False,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_vase_external_solid_contract_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "source": "examples/vase.py",
        "diagnosis": (
            "external_kernel_zero_volume_despite_positive_same_kernel_roundtrip"
            if not issues
            else "incomplete_diagnosis"
        ),
        "solver_ready": False,
        "notes": [
            "A named volume entity is not a positive-volume solid; measure mass before meshing.",
            "A same-kernel STEP roundtrip does not prove external-kernel portability.",
            "Heal and explicit Solid extraction are independent recovery probes, not acceptance evidence by themselves.",
        ],
    }


def build123d_tea_cup_source_contract_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Bind the upstream tea-cup source and its headless external diagnosis."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    build = _mapping(summary.get("build"), "build")
    external = _mapping(summary.get("external"), "external")
    source = _mapping(build.get("source_native_example"), "build.source_native_example")
    source_checks = _mapping(build.get("checks"), "build.checks")
    process = _mapping(external.get("process"), "external.process")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    source_sha = str(source.get("source_sha256", "")).lower()
    commit = str((build.get("versions") or {}).get("upstream_commit", "")).lower()
    features = {str(value) for value in (source.get("features") or [])}
    required_features = {
        "spline_profile_revolve",
        "open_face_shell_offset",
        "bottom_fusion",
        "edge_fillet",
        "intersection_derived_handle_contacts",
        "nonplanar_handle_sweep",
    }
    public_gate = build123d_curved_shell_step_semantics_gate(summary)
    checks = {
        "upstream_source_identity_bound": summary.get("source_kind")
        == "upstream_native_build123d_example"
        and source.get("source") == "examples/tea_cup.py"
        and len(source_sha) == 64
        and all(character in "0123456789abcdef" for character in source_sha)
        and len(commit) == 40
        and all(character in "0123456789abcdef" for character in commit),
        "viewer_only_stub_preserves_geometry": source.get("viewer_change")
        == "ocp_vscode.show stubbed; geometry executed unchanged",
        "source_feature_contract_complete": features == required_features,
        "source_immutable_and_official_assertion_reproduced": source_checks.get("source_immutable")
        is True
        and source_checks.get("official_volume_assertion_reproduced") is True,
        "same_kernel_diagnostic_completed": build.get("pass") is True
        and build.get("disposition") == "same_kernel_topology_preserved_but_mass_not_exact",
        "external_replay_is_headless_and_synchronous": external.get("execution_mode")
        == "python_api_headless_synchronous_commands"
        and {"-nographics", "-batch"}.issubset(set(external.get("headless_flags") or []))
        and external.get("persistent_gui_started") is False,
        "external_process_is_fresh_and_classified": process.get("acceptable") is True
        and process.get("result_artifact_fresh") is True
        and process.get("known_headless_diagnostics_only") is True
        and not list(process.get("unexpected_error_lines") or [])
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "exactly_four_timing_stages_recorded": len(timing) == 4
        and all(_finite(value, f"timing_breakdown_s.{name}") >= 0.0 for name, value in timing.items()),
        "public_semantics_gate_passed": public_gate["status"] == "ok",
        "negative_portability_outcome_not_overclaimed": public_gate["solver_ready"] is False,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_tea_cup_source_contract_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "source": "examples/tea_cup.py",
        "public_gate_status": public_gate["status"],
        "solver_ready": False,
        "notes": [
            "Execute upstream geometry unchanged except for a viewer stub and bind the source digest and commit.",
            "A validated negative portability result is useful MCP learning, but it is not a solver-ready handoff.",
            "Keep the original native mass properties as the authority until a repaired external representation passes independently.",
        ],
    }
