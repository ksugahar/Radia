"""Gates for a source-native hex/pyramid/tet export replay."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _mapping(summary: Mapping[str, object], name: str) -> Mapping[str, object]:
    value = summary.get(name)
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _counts(value: Mapping[str, object]) -> dict[str, int]:
    return {str(key).lower(): int(count) for key, count in value.items()}


def cubit_pyramid_mixed_export_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate explicit pyramid preservation across Gmsh and Nastran exports."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    counts = _counts(_mapping(summary, "element_counts"))
    connectivity = _mapping(summary, "connectivity_sizes")
    quality = _mapping(summary, "quality")
    gmsh = _mapping(summary, "gmsh_header")
    independent = _mapping(summary, "gmsh_independent")
    bdf = {str(key).upper(): int(value) for key, value in _mapping(summary, "nastran_cards").items()}
    companions = _mapping(summary, "gmsh_companions")

    gmsh_counts = _counts(_mapping(independent, "element_counts"))
    gmsh_volumes = {
        str(key).lower(): float(value)
        for key, value in _mapping(independent, "volume_by_type").items()
    }
    hex_quality = _mapping(quality, "hex")
    tet_quality = _mapping(quality, "tet")
    pyramid_quality = _mapping(quality, "pyramid")
    hex_scaled = _mapping(hex_quality, "scaled_jacobian")
    tet_scaled = _mapping(tet_quality, "scaled_jacobian")
    pyramid_scaled = _mapping(pyramid_quality, "scaled_jacobian")
    volume_error = float(summary.get("gmsh_cad_volume_relative_error", math.inf))

    checks = {
        "hex_pyramid_tet_transition_is_explicit": counts.get("hex", 0) > 0
        and counts.get("pyramid", 0) > 0
        and counts.get("tet", 0) > 0
        and counts.get("wedge", 0) == 0,
        "first_order_connectivity_is_exact": {
            key: list(connectivity.get(key, []))
            for key in ("hex", "pyramid", "tet", "wedge")
        }
        == {"hex": [8], "pyramid": [5], "tet": [4], "wedge": []},
        "supported_quality_is_positive": bool(hex_scaled.get("available"))
        and int(hex_scaled.get("count", 0)) == counts.get("hex", 0)
        and float(hex_scaled.get("min", math.nan)) > 0.0
        and bool(tet_scaled.get("available"))
        and int(tet_scaled.get("count", 0)) == counts.get("tet", 0)
        and float(tet_scaled.get("min", math.nan)) > 0.0,
        "pyramid_quality_falls_back_to_geometric_integration": pyramid_scaled.get(
            "available"
        )
        is False
        and int(pyramid_scaled.get("count", -1)) == 0
        and gmsh_volumes.get("pyramid", 0.0) > 0.0,
        "gmsh_is_complete_ascii_v41": gmsh.get("status") == "ok"
        and gmsh.get("mesh_format") == "4.1"
        and gmsh.get("binary") is False,
        "gmsh_counts_match_cubit": all(
            gmsh_counts.get(kind, -1) == counts.get(kind, 0)
            for kind in ("hex", "pyramid", "tet")
        )
        and gmsh_counts.get("other_3d", 0) == 0,
        "all_mixed_gmsh_volumes_are_positive": all(
            gmsh_volumes.get(kind, 0.0) > 0.0 for kind in ("hex", "pyramid", "tet")
        ),
        "gmsh_volume_matches_cad": math.isfinite(volume_error) and volume_error <= 1.0e-2,
        "nastran_preserves_explicit_cpyram": bdf.get("CHEXA", 0)
        == counts.get("hex", 0)
        and bdf.get("CPYRAM", 0) == counts.get("pyramid", 0)
        and bdf.get("CTETRA", 0) == counts.get("tet", 0),
        "gmsh_launch_companions_recorded": all(
            companions.get(name) is True
            for name in ("geo", "geo_opt", "msh_opt")
        ),
    }
    return {
        "policy": "cubit_pyramid_mixed_export_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "element_counts": counts,
            "gmsh_volume_by_type": gmsh_volumes,
            "gmsh_cad_volume_relative_error": volume_error,
            "hex_scaled_jacobian_minimum": float(hex_scaled.get("min", math.nan)),
            "tet_scaled_jacobian_minimum": float(tet_scaled.get("min", math.nan)),
            "nastran_card_counts": {
                key: bdf.get(key, 0) for key in ("CHEXA", "CPYRAM", "CTETRA")
            },
        },
        "lesson": (
            "Keep pyramids explicit when validating a conformal hex-to-tet transition. "
            "Match Cubit counts to Gmsh type 5/7/4 blocks, integrate each family, and "
            "require default Nastran export to retain CPYRAM rather than silently "
            "dropping or collapsing transition cells."
        ),
    }


def cubit_pyramid_source_plugin_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Gate source migration and the executable-owned plugin startup path."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    process = _mapping(summary, "process")
    transforms = _mapping(summary, "compatibility_transforms")
    timing = _mapping(summary, "timing_breakdown_s")
    public = _mapping(summary, "public_gate")
    flags = {str(value).lower() for value in summary.get("headless_flags", [])}
    known_errors = process.get("error_lines", [])
    unexpected = process.get("unexpected_error_lines", [])
    if isinstance(known_errors, (str, bytes)) or not isinstance(known_errors, Sequence):
        raise ValueError("process.error_lines must be a sequence")
    if isinstance(unexpected, (str, bytes)) or not isinstance(unexpected, Sequence):
        raise ValueError("process.unexpected_error_lines must be a sequence")

    checks = {
        "source_native_python_digest_recorded": summary.get("source_native_example")
        is True
        and summary.get("source_extension") == ".py"
        and len(str(summary.get("source_sha256", ""))) == 64,
        "combined_journal_headless_runtime": summary.get("execution_mode")
        == "combined_journal_headless"
        and {"-nographics", "-batch"}.issubset(flags)
        and summary.get("gui_daemon_enabled") is False,
        "plugin_registration_path_was_demonstrated": summary.get(
            "ordinary_python_plugin_command_rejected"
        )
        is True
        and summary.get("combined_journal_plugin_command_succeeded") is True,
        "runtime_and_export_migrations_recorded": "2024.3 Python init"
        in str(transforms.get("cubit_runtime", ""))
        and "combined journal startup" in str(transforms.get("cubit_runtime", ""))
        and "export_3D_Nastran" in str(transforms.get("legacy_nastran_function", ""))
        and "export jmag_nastran" in str(transforms.get("legacy_nastran_function", "")),
        "current_cubit_version_recorded": str(summary.get("version", "")).startswith(
            "2025.12"
        ),
        "known_nonzero_headless_exit_classified": int(process.get("exit_code", -1))
        == 3
        and process.get("known_headless_diagnostics_only") is True
        and len(known_errors) > 0
        and len(unexpected) == 0
        and process.get("acceptable") is True,
        "fresh_artifact_and_no_process_leak": process.get("result_artifact_fresh")
        is True
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "source_commands_replayed": summary.get("source_commands_present") is True,
        "public_mixed_export_gate_passed": public.get("policy")
        == "cubit_pyramid_mixed_export_gate_v1"
        and public.get("status") == "ok",
        "public_negative_control_rejected": summary.get("public_negative_status")
        == "needs_attention",
        "four_stage_timing_recorded": len(timing) == 4
        and all(float(value) >= 0.0 for value in timing.values()),
    }
    return {
        "policy": "cubit_pyramid_source_plugin_replay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "launcher_classification": (
            "known_headless_diagnostics" if checks["known_nonzero_headless_exit_classified"] else "execution_error"
        ),
        "lesson": (
            "Plugin export commands are registered by the Coreform executable startup "
            "path, not by an ordinary Python cubit.init import. Replay legacy Python "
            "examples through a combined #!cubit/#!python journal, keep the run headless, "
            "and accept a nonzero exit only with fresh artifacts and allowlisted diagnostics."
        ),
    }
