"""Gates for a conformal Cubit hex-pyramid-tet transition cell."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path


_VOLUME_FAMILIES = ("hex", "pyramid", "tet", "wedge")


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
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a finite number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{field} must be a finite number")
    return number


def _counts(value: object, field: str) -> dict[str, int]:
    row = _mapping(value, field)
    result = {family: int(row.get(family, 0)) for family in _VOLUME_FAMILIES}
    if any(count < 0 for count in result.values()):
        raise ValueError(f"{field} counts must be nonnegative")
    return result


def cubit_conformal_hex_pyramid_tet_interface_gate(
    summary: Mapping[str, object],
    *,
    mapped_volume_id: int = 1,
    transition_volume_id: int = 2,
    min_scaled_jacobian: float = 0.1,
    max_volume_relative_error: float = 1.0e-9,
) -> dict[str, object]:
    """Validate topology, adjacency, quality, and independent volume closure.

    The gate intentionally does not require hexes to outnumber tetrahedra.  A
    minimal conformal transition can contain one mapped hex, one pyramid, and
    several tetrahedra; ownership of the shared quad is the stronger contract.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    threshold = _finite(min_scaled_jacobian, "min_scaled_jacobian")
    tolerance = _finite(max_volume_relative_error, "max_volume_relative_error")
    if threshold <= 0.0 or tolerance < 0.0:
        raise ValueError("quality threshold must be positive and volume tolerance nonnegative")
    mapped_id = str(int(mapped_volume_id))
    transition_id = str(int(transition_volume_id))
    if mapped_id == transition_id:
        raise ValueError("mapped and transition volume IDs must differ")

    totals = _counts(summary.get("element_counts"), "element_counts")
    per_volume_raw = _mapping(summary.get("per_volume_element_counts"), "per_volume_element_counts")
    per_volume = {
        str(volume_id): _counts(row, f"per_volume_element_counts[{volume_id!r}]")
        for volume_id, row in per_volume_raw.items()
    }
    if mapped_id not in per_volume or transition_id not in per_volume:
        raise ValueError("mapped and transition volume inventories are required")

    quality_raw = _mapping(summary.get("quality"), "quality")
    quality_counts: dict[str, int] = {}
    quality_minima: dict[str, float] = {}
    for family in ("hex", "pyramid", "tet"):
        family_row = _mapping(quality_raw.get(family), f"quality.{family}")
        metric = _mapping(family_row.get("scaled_jacobian"), f"quality.{family}.scaled_jacobian")
        quality_counts[family] = int(metric.get("count", -1))
        quality_minima[family] = _finite(metric.get("min"), f"quality.{family}.scaled_jacobian.min")

    interfaces = _rows(summary.get("interface_surfaces"), "interface_surfaces")
    ownership = _rows(summary.get("interface_face_ownership"), "interface_face_ownership")
    interface_face_ids = {
        int(face_id)
        for row in interfaces
        for face_id in (row.get("face_ids") or [])
    }
    ownership_face_ids = {int(row.get("face_id", -1)) for row in ownership}
    all_interface_connectivity = [
        list(connectivity)
        for row in interfaces
        for connectivity in (row.get("face_connectivity") or [])
    ]
    interface_face_incidence_counts = [
        int(row.get("face_incidence_count", 2)) for row in interfaces
    ]

    geometry = _mapping(summary.get("geometry"), "geometry")
    cad_total = _finite(geometry.get("cad_total_volume_m3"), "geometry.cad_total_volume_m3")
    analytic_total = _finite(
        geometry.get("analytic_total_volume_m3"), "geometry.analytic_total_volume_m3"
    )
    gmsh_inventory = _mapping(summary.get("gmsh_inventory"), "gmsh_inventory")
    gmsh_counts = _counts(gmsh_inventory.get("volume_family_counts"), "gmsh_inventory.volume_family_counts")
    gmsh_volume = _mapping(summary.get("gmsh_volume_inventory"), "gmsh_volume_inventory")
    reconstructed_counts = _counts(gmsh_volume.get("family_counts"), "gmsh_volume_inventory.family_counts")
    reconstructed_total = _finite(
        gmsh_volume.get("total_volume_m3"), "gmsh_volume_inventory.total_volume_m3"
    )
    family_volumes = _mapping(gmsh_volume.get("family_volumes_m3"), "gmsh_volume_inventory.family_volumes_m3")
    reconstructed_family_sum = sum(
        _finite(family_volumes.get(family), f"gmsh_volume_inventory.family_volumes_m3.{family}")
        for family in ("hex", "pyramid", "tet")
    )

    denominator = max(abs(cad_total), 1.0e-300)
    analytic_error = abs(cad_total - analytic_total) / max(abs(analytic_total), 1.0e-300)
    reconstructed_error = abs(reconstructed_total - cad_total) / denominator
    reconstructed_sum_error = abs(reconstructed_family_sum - reconstructed_total) / max(
        abs(reconstructed_total), 1.0e-300
    )
    per_volume_sum = {
        family: sum(row[family] for row in per_volume.values()) for family in _VOLUME_FAMILIES
    }
    mapped = per_volume[mapped_id]
    transition = per_volume[transition_id]
    export = _mapping(summary.get("gmsh_export"), "gmsh_export")

    checks = {
        "two_distinct_partition_volumes_recorded": set(per_volume) == {mapped_id, transition_id},
        "mapped_volume_is_hex_only": mapped["hex"] > 0
        and all(mapped[family] == 0 for family in ("pyramid", "tet", "wedge")),
        "transition_volume_is_tet_plus_pyramid": transition["tet"] > 0
        and transition["pyramid"] > 0
        and transition["hex"] == 0
        and transition["wedge"] == 0,
        "per_volume_inventory_matches_totals": per_volume_sum == totals,
        "quality_count_matches_topology": all(
            quality_counts[family] == totals[family] for family in ("hex", "pyramid", "tet")
        ),
        "all_volume_families_above_quality_threshold": all(
            quality_minima[family] >= threshold for family in ("hex", "pyramid", "tet")
        ),
        "single_shared_interface_between_partitions": len(interfaces) == 1
        and {int(value) for value in (interfaces[0].get("adjacent_volumes") or [])}
        == {int(mapped_id), int(transition_id)},
        "interface_faces_are_quads": bool(interface_face_ids)
        and len(all_interface_connectivity) == len(interface_face_ids)
        and all(len(connectivity) == 4 for connectivity in all_interface_connectivity),
        "each_interface_quad_has_one_hex_and_one_pyramid_owner": bool(ownership)
        and ownership_face_ids == interface_face_ids
        and all(
            int(row.get("node_count", 0)) == 4
            and len(list(row.get("hex_owners") or [])) == 1
            and len(list(row.get("pyramid_owners") or [])) == 1
            for row in ownership
        ),
        "interface_quads_are_two_sided_manifold": bool(interfaces)
        and all(count == 2 for count in interface_face_incidence_counts)
        and all(
            not list(row.get("tet_owners") or [])
            and not list(row.get("wedge_owners") or [])
            and not list(row.get("other_owners") or [])
            for row in ownership
        ),
        "every_pyramid_serves_the_transition": int(summary.get("matched_pyramid_count", -1))
        == totals["pyramid"],
        "gmsh_export_is_fresh_nonempty_digest": int(export.get("bytes", 0)) > 0
        and len(str(export.get("sha256", ""))) == 64,
        "gmsh_is_ascii_v41_with_valid_connectivity": gmsh_inventory.get("status") == "ok"
        and str(gmsh_inventory.get("mesh_format", "")) == "4.1"
        and gmsh_inventory.get("binary") is False
        and not list(gmsh_inventory.get("connectivity_mismatches") or []),
        "gmsh_volume_families_match_live_inventory": gmsh_counts == totals
        and reconstructed_counts == totals,
        "cad_volume_matches_analytic": analytic_error <= tolerance,
        "independent_gmsh_volume_sum_matches_cad": reconstructed_error <= tolerance
        and reconstructed_sum_error <= tolerance
        and geometry.get("element_volume_source")
        == "independent_gmsh_v41_coordinate_reconstruction",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_conformal_hex_pyramid_tet_interface_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "element_counts": totals,
        "quality_minima": quality_minima,
        "interface_face_count": len(interface_face_ids),
        "matched_pyramid_count": int(summary.get("matched_pyramid_count", -1)),
        "cad_volume_relative_error": analytic_error,
        "gmsh_reconstructed_volume_relative_error": reconstructed_error,
        "notes": [
            "Do not require hex dominance: a valid minimal transition may contain more tetrahedra than hexes.",
            "A conformal transition is proven by quad ownership on both sides, not by family counts alone.",
            "Reconstruct volume from Gmsh 4.1 coordinates so CAD closure is independent of Cubit's quality API.",
        ],
    }


def cubit_mixed_transition_source_gate(
    summary: Mapping[str, object],
    *,
    mapped_volume_id: int = 1,
    transition_volume_id: int = 2,
) -> dict[str, object]:
    """Gate the source journal and a classified synchronous headless replay."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    commands_raw = summary.get("source_commands")
    if isinstance(commands_raw, (str, bytes)) or not isinstance(commands_raw, Sequence):
        raise ValueError("source_commands must be a sequence")
    commands = [str(command).strip().lower() for command in commands_raw if str(command).strip()]
    command_text = "\n".join(commands)

    def index(fragment: str) -> int:
        return next((offset for offset, command in enumerate(commands) if fragment in command), -1)

    source_sha = str(summary.get("source_sha256", "")).lower()
    quality_probe = _mapping(summary.get("quality_probe"), "quality_probe")
    process = _mapping(summary.get("process"), "process")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    export_artifacts_value = summary.get("export_artifacts")
    export_artifacts_present = isinstance(export_artifacts_value, Mapping)
    export_artifacts_ok = True
    if export_artifacts_present:
        required_raw = export_artifacts_value.get("required")
        if isinstance(required_raw, (str, bytes)) or not isinstance(
            required_raw, Sequence
        ):
            raise ValueError("export_artifacts.required must be a sequence")
        required = [str(name) for name in required_raw]
        artifact_rows = _rows(
            export_artifacts_value.get("artifacts"), "export_artifacts.artifacts"
        )
        artifact_names = [str(row.get("name", "")) for row in artifact_rows]
        required_names = set(required)
        artifact_by_name = dict(zip(artifact_names, artifact_rows, strict=True))
        export_artifacts_ok = (
            bool(required_names)
            and all(required)
            and all(artifact_names)
            and len(required_names) == len(required)
            and len(artifact_by_name) == len(artifact_rows)
            and set(artifact_by_name) == required_names
        )
        export_artifacts_ok = export_artifacts_ok and all(
            row.get("fresh") is True
            and int(row.get("bytes", 0)) > 0
            and len(str(row.get("sha256", ""))) == 64
            for row in artifact_by_name.values()
        )

    replay_identity_value = summary.get("replay_identity")
    replay_identity_present = isinstance(replay_identity_value, Mapping)
    replay_identity_ok = True
    if replay_identity_present:
        pinned_journal = str(replay_identity_value.get("pinned_journal_sha256", ""))
        pinned_model = str(replay_identity_value.get("pinned_source_model_sha256", ""))
        replayed_journal = str(
            replay_identity_value.get("replayed_journal_sha256", "")
        )
        replayed_model = str(
            replay_identity_value.get("replayed_source_model_sha256", "")
        )
        replay_identity_ok = (
            all(
                len(value) == 64
                and all(character.lower() in "0123456789abcdef" for character in value)
                for value in (
                    pinned_journal,
                    pinned_model,
                    replayed_journal,
                    replayed_model,
                )
            )
            and pinned_journal == replayed_journal
            and pinned_model == replayed_model
            and pinned_journal == str(summary.get("source_sha256", ""))
        )
    public_gate = cubit_conformal_hex_pyramid_tet_interface_gate(
        summary,
        mapped_volume_id=mapped_volume_id,
        transition_volume_id=transition_volume_id,
    )
    block_families = ("hex", "pyramid", "tet", "tri", "face")
    ordered = [
        index("brick x 2 y 1 z 1"),
        index("webcut volume 1"),
        index("volume 1 scheme map"),
        index("mesh volume 1"),
        index("volume 2 scheme tetmesh"),
        index("mesh volume 2"),
    ]
    exit_code = int(process.get("exit_code", -1))
    unexpected = list(process.get("unexpected_error_lines") or [])
    known_diagnostics_only = process.get("known_headless_diagnostics_only") is True
    exit_explained = exit_code == 0 or (
        exit_code > 0
        and known_diagnostics_only
        and not unexpected
        and process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0
        and public_gate["status"] == "ok"
    )
    checks = {
        "source_native_journal_and_sha256_recorded": str(summary.get("source_kind", "")).startswith(
            "source_native_"
        )
        and Path(str(summary.get("source_journal", ""))).suffix.lower() == ".jou"
        and len(source_sha) == 64
        and all(character in "0123456789abcdef" for character in source_sha),
        "source_builds_and_webcuts_brick": ordered[0] >= 0
        and ordered[1] > ordered[0]
        and "xplane" in commands[ordered[1]],
        "source_meshes_map_before_tetmesh": all(offset >= 0 for offset in ordered)
        and ordered == sorted(ordered),
        "source_registers_all_mixed_and_surface_families": all(
            f"add {family} all" in command_text for family in block_families
        ),
        "source_applies_meter_scale": any(
            "volume all scale 0.001" in command for command in commands
        ),
        "headless_batch_without_gui_daemon": summary.get("execution_mode")
        == "headless_combined_journal_then_python_inventory"
        and {"-nographics", "-batch"}.issubset(set(summary.get("headless_flags") or []))
        and summary.get("gui_daemon_enabled") is False,
        "cubit_version_recorded": bool(str(summary.get("version", "")).strip()),
        "unsupported_aggregate_quality_probe_is_diagnosed": quality_probe.get("command_supported")
        is False
        and "unknown metric name volume" in str(quality_probe.get("diagnostic", "")).lower()
        and quality_probe.get("failure_interpretation") == "unsupported_api_not_zero_quality",
        "per_element_scaled_jacobian_fallback_recorded": quality_probe.get("fallback")
        == "per_element_scaled_jacobian_by_family"
        and {str(value) for value in (quality_probe.get("families") or [])}
        == {"hex", "pyramid", "tet"},
        "nonzero_exit_is_semantically_classified": exit_explained,
        "fresh_artifact_and_no_owned_process_leak": process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "all_required_export_artifacts_are_fresh": export_artifacts_ok,
        "journal_and_source_model_identity_match_replay": replay_identity_ok,
        "exactly_four_timing_stages_recorded": len(timing) == 4
        and all(_finite(value, f"timing_breakdown_s.{name}") >= 0.0 for name, value in timing.items()),
        "independent_interface_gate_passed": public_gate["status"] == "ok",
    }
    issues = [name for name, ok in checks.items() if not ok]
    warnings = []
    if not export_artifacts_present:
        warnings.append("per_artifact_export_freshness_not_recorded")
    if not replay_identity_present:
        warnings.append("journal_model_replay_identity_not_recorded")
    return {
        "policy": "cubit_mixed_transition_source_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "warnings": warnings,
        "source_journal": Path(str(summary.get("source_journal", ""))).name,
        "process_exit_code": exit_code,
        "public_gate_status": public_gate["status"],
        "notes": [
            "Replay source commands synchronously and headlessly; never infer completion from a queued GUI playback.",
            "An unsupported aggregate quality query is an API diagnostic, not evidence of zero elements or zero quality.",
            "Permit a nonzero batch exit only with classified diagnostics, a fresh passing artifact, and no leaked process.",
        ],
    }
