"""Gates for a symmetric swept hex/pyramid/tet Cubit mesh."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def cubit_symmetric_swept_mixed_mesh_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate symmetry, transition topology, quality, and independent export volume."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    geometry = _mapping(summary.get("geometry"), "geometry")
    counts = _mapping(summary.get("element_counts"), "element_counts")
    connectivity = _mapping(summary.get("connectivity_sizes"), "connectivity_sizes")
    quality = _mapping(summary.get("quality"), "quality")
    gmsh = _mapping(summary.get("gmsh"), "gmsh")
    gmsh_counts = _mapping(gmsh.get("element_counts"), "gmsh.element_counts")
    gmsh_volumes = _mapping(gmsh.get("volume_by_type"), "gmsh.volume_by_type")

    hex_count = int(counts.get("hex", 0))
    tet_count = int(counts.get("tet", 0))
    pyramid_count = int(counts.get("pyramid", 0))
    wedge_count = int(counts.get("wedge", 0))
    symmetry_error = float(geometry.get("left_right_aggregate_relative_error", math.inf))
    volume_error = float(gmsh.get("cad_volume_relative_error", math.inf))
    hex_quality = float(_mapping(quality.get("hex_scaled_jacobian"), "quality.hex_scaled_jacobian").get("min", math.nan))
    tet_quality = float(_mapping(quality.get("tet_scaled_jacobian"), "quality.tet_scaled_jacobian").get("min", math.nan))
    checks = {
        "five_volume_swept_partition": int(geometry.get("volume_count", 0)) == 5
        and float(geometry.get("cad_total_volume", 0.0)) > 0.0,
        "mirrored_aggregate_cad_volume_closes": math.isfinite(symmetry_error)
        and symmetry_error <= 1.0e-12,
        "hex_pyramid_tet_transition_is_explicit": hex_count > 0
        and pyramid_count > 0
        and tet_count > 0
        and wedge_count == 0,
        "connectivity_matches_first_order_families": connectivity.get("hex") == [8]
        and connectivity.get("pyramid") == [5]
        and connectivity.get("tet") == [4],
        "hex_and_tet_quality_are_positive": math.isfinite(hex_quality)
        and hex_quality > 0.0
        and math.isfinite(tet_quality)
        and tet_quality > 0.0,
        "gmsh_ascii_v41_complete": gmsh.get("mesh_format") == "4.1"
        and gmsh.get("binary") is False,
        "gmsh_topology_matches_cubit": int(gmsh_counts.get("hex", -1)) == hex_count
        and int(gmsh_counts.get("pyramid", -1)) == pyramid_count
        and int(gmsh_counts.get("tet", -1)) == tet_count
        and int(gmsh_counts.get("other_3d", -1)) == 0,
        "every_exported_volume_family_has_positive_volume": all(
            float(gmsh_volumes.get(kind, 0.0)) > 0.0 for kind in ("hex", "pyramid", "tet")
        ),
        "gmsh_integrated_volume_matches_cad": math.isfinite(volume_error)
        and volume_error <= 1.0e-3,
    }
    return {
        "policy": "cubit_symmetric_swept_mixed_mesh_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "hex_count": hex_count,
            "pyramid_count": pyramid_count,
            "tet_count": tet_count,
            "left_right_aggregate_relative_error": symmetry_error,
            "minimum_hex_scaled_jacobian": hex_quality,
            "minimum_tet_scaled_jacobian": tet_quality,
            "gmsh_cad_volume_relative_error": volume_error,
        },
        "lesson": (
            "For a symmetric swept profile, compare aggregate half-model CAD volumes before "
            "judging the mesh. Mesh sweepable regions with hexes, then allow the complex region "
            "to use explicit pyramids between quad faces and tet interiors. Validate all three "
            "families again from the exported Gmsh connectivity and integrated volume."
        ),
    }


def cubit_symmetric_swept_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate source journal, headless diagnostics, timing, and public closure."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    process = _mapping(summary.get("process"), "process")
    public = _mapping(summary.get("public_gate"), "public_gate")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    source_sha = str(summary.get("source_sha256", "")).lower()
    checks = {
        "source_native_journal_and_promotion_recorded": summary.get("source_native_journal") is True
        and summary.get("promotion") == "mirrored_quad_sections_to_swept_hex_pyramid_tet",
        "source_digest_is_sha256": len(source_sha) == 64
        and all(character in "0123456789abcdef" for character in source_sha),
        "headless_python_journal_route_recorded": summary.get("execution_mode")
        == "headless_combined_journal_then_python_inventory"
        and {"-nographics", "-batch"}.issubset(set(summary.get("headless_flags") or []))
        and summary.get("gui_daemon_enabled") is False,
        "classified_exit_has_fresh_complete_artifact": process.get("acceptable") is True
        and process.get("result_artifact_fresh") is True
        and not list(process.get("unexpected_error_lines") or [])
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "known_diagnostics_do_not_hide_mesh_failure": process.get(
            "known_headless_diagnostics_only"
        )
        is True
        and int(process.get("exit_code", -1)) in {0, 3, 4},
        "exactly_four_nonnegative_timing_stages": len(timing) == 4
        and all(math.isfinite(float(value)) and float(value) >= 0.0 for value in timing.values()),
        "public_mixed_mesh_gate_passed": public.get("policy")
        == "cubit_symmetric_swept_mixed_mesh_gate_v1"
        and public.get("status") == "ok",
        "public_negative_control_rejected": summary.get("public_negative_status")
        == "needs_attention",
    }
    return {
        "policy": "cubit_symmetric_swept_source_replay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "notes": [
            "A license diagnostic alone is not completion evidence; require a fresh numerical artifact.",
            "Treat the pyramid family as explicit topology even when a direct Cubit pyramid quality array is unavailable.",
        ],
    }
