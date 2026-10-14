"""Validation gates for equivalent Builder/Algebra patterned prisms."""

from __future__ import annotations

import math
from typing import Mapping


def _relative_error(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), 1.0)


def _extent_error(a, b) -> float:
    return max((_relative_error(x, y) for x, y in zip(a, b)), default=math.inf)


def _topology(row: Mapping[str, object]) -> tuple[int, int, int, int]:
    return tuple(int(row.get(key, -1)) for key in ("solid_count", "face_count", "edge_count", "vertex_count"))


def _external_topology(row: Mapping[str, object]) -> tuple[int, int, int, int]:
    return tuple(int(row.get(key, -1)) for key in ("volume_count", "surface_count", "curve_count", "vertex_count"))


def dual_api_prismatic_pattern_gate(
    summary: Mapping[str, object],
    *,
    source_relative_tolerance: float = 1.0e-12,
    external_volume_relative_tolerance: float = 2.0e-6,
    external_pair_relative_tolerance: float = 5.0e-7,
) -> dict[str, object]:
    """Gate API parity separately from external STEP-kernel integration bias."""

    native = summary.get("native") or {}
    external = summary.get("external") or {}
    native_rows = native.get("records") or {}
    external_rows = external.get("records") or {}
    builder = native_rows.get("builder") or {}
    algebra = native_rows.get("algebra") or {}
    centered = native_rows.get("algebra_centered") or {}
    ext_builder = external_rows.get("builder") or {}
    ext_algebra = external_rows.get("algebra") or {}
    ext_centered = external_rows.get("algebra_centered") or {}
    expected_volume = float(native.get("expected_volume", math.nan))
    source_tol = float(source_relative_tolerance)
    external_tol = float(external_volume_relative_tolerance)
    pair_tol = float(external_pair_relative_tolerance)

    native_volume_error = _relative_error(builder.get("volume", math.nan), algebra.get("volume", math.nan))
    native_area_error = _relative_error(builder.get("area", math.nan), algebra.get("area", math.nan))
    native_extent_error = _extent_error(builder.get("bbox_extent", []), algebra.get("bbox_extent", []))
    external_errors = {
        name: _relative_error(row.get("volume", math.nan), native_rows.get(name, {}).get("volume", math.nan))
        for name, row in external_rows.items()
        if name in native_rows
    }
    external_pair_error = _relative_error(ext_builder.get("volume", math.nan), ext_algebra.get("volume", math.nan))
    checks = {
        "official_expected_volume_reproduced": all(
            _relative_error(row.get("volume", math.nan), expected_volume) <= source_tol
            for row in (builder, algebra)
        ),
        "native_dual_api_volume_match": native_volume_error <= source_tol,
        "native_dual_api_area_match": native_area_error <= source_tol,
        "native_dual_api_bbox_extent_match": native_extent_error <= source_tol,
        "native_dual_api_topology_match": _topology(builder) == _topology(algebra) and min(_topology(builder)) > 0,
        "native_shapes_valid": all(row.get("is_valid") is True for row in (builder, algebra, centered)),
        "native_step_roundtrips_valid": all((row.get("self_roundtrip") or {}).get("is_valid") is True for row in (builder, algebra, centered)),
        "centering_preserves_native_mass_topology": _relative_error(algebra.get("volume", math.nan), centered.get("volume", math.nan)) <= source_tol and _topology(algebra) == _topology(centered),
        "external_all_single_solid": all(int(row.get("volume_count", 0)) == 1 for row in (ext_builder, ext_algebra, ext_centered)),
        "external_topology_matches_across_exports": _external_topology(ext_builder) == _external_topology(ext_algebra) == _external_topology(ext_centered) and min(_external_topology(ext_builder)) > 0,
        "external_bbox_extents_match_native": all(_extent_error(row.get("bbox_extent", []), native_rows[name].get("bbox_extent", [])) <= source_tol for name, row in external_rows.items() if name in native_rows),
        "external_volumes_within_kernel_tolerance": len(external_errors) == 3 and all(error <= external_tol for error in external_errors.values()),
        "external_dual_api_pair_within_tolerance": external_pair_error <= pair_tol,
        "centering_does_not_hide_step_history_bias": _relative_error(ext_algebra.get("volume", math.nan), ext_centered.get("volume", math.nan)) <= source_tol,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_dual_api_prismatic_pattern_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "native_dual_api_volume_relative_error": native_volume_error,
            "native_dual_api_area_relative_error": native_area_error,
            "native_dual_api_bbox_extent_max_relative_error": native_extent_error,
            "external_volume_relative_errors": external_errors,
            "external_dual_api_pair_relative_error": external_pair_error,
            "external_centering_relative_change": _relative_error(ext_algebra.get("volume", math.nan), ext_centered.get("volume", math.nan)),
        },
        "tolerances": {
            "native_source_relative": source_tol,
            "external_volume_relative": external_tol,
            "external_pair_relative": pair_tol,
        },
        "notes": [
            "Use machine-precision parity for equivalent build123d APIs before STEP export.",
            "Use an explicit ppm-scale external-kernel volume tolerance while keeping body/topology and bbox contracts exact.",
            "A centered translation is a useful control, but it does not erase STEP construction-history differences.",
        ],
    }


def dual_api_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate immutable upstream dual-API execution and headless external replay."""

    sources = summary.get("sources") or {}
    artifacts = summary.get("artifacts") or {}
    diagnostics = [str(value) for value in summary.get("startup_diagnostics", [])]
    script_errors = [str(value) for value in summary.get("script_error_lines", [])]
    expected_suffixes = {"/plugins", "-commandplugindir", "-nojournal"}
    observed_suffixes = {
        suffix
        for row in diagnostics
        for suffix in expected_suffixes
        if row.rstrip().endswith(suffix)
    }
    process_exit_code = int(summary.get("external_process_exit_code", -1))
    process_ok = process_exit_code == 0 or (
        process_exit_code in {2, 3}
        and len(diagnostics) == 3
        and all("Could not open file:" in row for row in diagnostics)
        and observed_suffixes == expected_suffixes
        and not script_errors
        and summary.get("result_artifact_fresh") is True
    )
    checks = {
        "upstream_tag_commit_recorded": len(str(summary.get("upstream_commit", ""))) == 40,
        "builder_source_digest_recorded": len(str((sources.get("builder") or {}).get("sha256", ""))) == 64,
        "algebra_source_digest_recorded": len(str((sources.get("algebra") or {}).get("sha256", ""))) == 64,
        "source_files_preserved": all((sources.get(name) or {}).get("preserved") is True for name in ("builder", "algebra")),
        "official_assertion_reproduced": summary.get("official_assertion_reproduced") is True,
        "viewer_suppressed_for_headless_replay": summary.get("viewer_suppressed") is True,
        "source_replayed_without_rewrite": summary.get("source_replay_mode") == "runpy_with_viewer_stub",
        "three_step_digests_recorded": len(artifacts) == 3 and all(len(str(row.get("sha256", ""))) == 64 for row in artifacts.values()),
        "derived_centering_marked": (artifacts.get("algebra_centered") or {}).get("derived_control") is True,
        "external_headless": summary.get("external_execution_mode") == "python_api_headless" and {"-nographics", "-batch"}.issubset(set(summary.get("headless_flags") or [])),
        "persistent_gui_disabled": summary.get("gui_daemon_enabled") is False,
        "external_process_acceptable": process_ok,
        "owned_processes_closed": int(summary.get("owned_processes_remaining", 0)) == 0,
        "public_gate_passed": summary.get("public_gate_status") == "ok",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_dual_api_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "process": {
            "exit_code": process_exit_code,
            "startup_only_allowlisted": observed_suffixes == expected_suffixes and len(diagnostics) == 3,
            "script_errors": script_errors,
            "result_artifact_fresh": summary.get("result_artifact_fresh") is True,
        },
    }
