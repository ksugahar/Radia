"""Validation gates for region-owned hex/tet/pyramid Cubit meshes."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite
from pathlib import Path


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def cubit_region_owned_mixed_mesh_gate(
    summary: Mapping[str, object],
    *,
    min_scaled_jacobian: float = 0.1,
    max_volume_relative_error: float = 5.0e-8,
) -> dict[str, object]:
    """Gate a hex-owned region embedded in tet air with pyramid transitions."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    quality_limit = float(min_scaled_jacobian)
    volume_limit = float(max_volume_relative_error)
    if not isfinite(quality_limit) or quality_limit <= 0.0:
        raise ValueError("min_scaled_jacobian must be finite and > 0")
    if not isfinite(volume_limit) or volume_limit < 0.0:
        raise ValueError("max_volume_relative_error must be finite and nonnegative")

    counts_raw = _mapping(summary.get("element_counts") or {}, "element_counts")
    per_volume_raw = _mapping(
        summary.get("per_volume_element_counts") or {},
        "per_volume_element_counts",
    )
    quality_raw = _mapping(summary.get("quality") or {}, "quality")
    geometry = _mapping(summary.get("geometry") or {}, "geometry")
    export = _mapping(summary.get("gmsh_export") or {}, "gmsh_export")
    header = _mapping(export.get("header") or {}, "gmsh_export.header")
    inventory = _mapping(summary.get("gmsh_inventory") or {}, "gmsh_inventory")
    inventory_counts_raw = _mapping(
        inventory.get("volume_family_counts") or {},
        "gmsh_inventory.volume_family_counts",
    )
    inventory_counts = {
        str(key).lower(): int(value) for key, value in inventory_counts_raw.items()
    }

    counts = {str(key).lower(): int(value) for key, value in counts_raw.items()}
    if any(value < 0 for value in counts.values()):
        raise ValueError("element counts must be nonnegative")
    conductor_volumes = [int(value) for value in summary.get("conductor_volumes") or []]
    air_volume = int(summary.get("air_volume", 0))
    if not conductor_volumes or air_volume <= 0 or air_volume in conductor_volumes:
        raise ValueError("conductor_volumes and a distinct positive air_volume are required")

    per_volume: dict[int, dict[str, int]] = {}
    for volume in conductor_volumes + [air_volume]:
        row = _mapping(
            per_volume_raw.get(str(volume), per_volume_raw.get(volume, {})),
            f"per_volume_element_counts[{volume}]",
        )
        per_volume[volume] = {
            kind: int(row.get(kind, 0)) for kind in ("hex", "tet", "pyramid", "wedge")
        }
        if any(value < 0 for value in per_volume[volume].values()):
            raise ValueError("per-volume element counts must be nonnegative")

    quality_minima: dict[str, float] = {}
    quality_counts: dict[str, int] = {}
    for kind in ("hex", "tet", "pyramid"):
        family = _mapping(quality_raw.get(kind) or {}, f"quality.{kind}")
        if "scaled_jacobian" not in family:
            raise ValueError(f"quality.{kind}.scaled_jacobian must be a mapping")
        metric = _mapping(
            family.get("scaled_jacobian") or {},
            f"quality.{kind}.scaled_jacobian",
        )
        if "min" not in metric or "count" not in metric:
            raise ValueError(
                f"quality.{kind}.scaled_jacobian requires min and count"
            )
        quality_minima[kind] = float(metric.get("min", float("nan")))
        quality_counts[kind] = int(metric.get("count", -1))

    interface_rows = list(summary.get("conductor_air_interfaces") or [])
    if not all(isinstance(row, Mapping) for row in interface_rows):
        raise ValueError("conductor_air_interfaces must contain mappings")
    interface_by_conductor: dict[int, int] = {volume: 0 for volume in conductor_volumes}
    interfaces_conformal = True
    for row in interface_rows:
        adjacent = {int(value) for value in row.get("adjacent_volumes") or []}
        owners = adjacent.intersection(conductor_volumes)
        if air_volume not in adjacent or len(adjacent) != 2 or len(owners) != 1:
            interfaces_conformal = False
            continue
        owner = next(iter(owners))
        interface_by_conductor[owner] += 1
        face_count = int(row.get("face_count", 0))
        interfaces_conformal = interfaces_conformal and (
            face_count > 0
            and face_count == int(row.get("quad_count", 0))
            and int(row.get("tri_count", 0)) == 0
            and float(row.get("area", 0.0)) > 0.0
        )

    volume_error = float(geometry.get("volume_relative_error", float("inf")))
    checks = {
        "conductor_regions_are_hex_only": all(
            per_volume[volume]["hex"] > 0
            and all(per_volume[volume][kind] == 0 for kind in ("tet", "pyramid", "wedge"))
            for volume in conductor_volumes
        ),
        "air_region_is_tet_with_pyramid_transition": per_volume[air_volume]["tet"] > 0
        and per_volume[air_volume]["pyramid"] > 0
        and per_volume[air_volume]["hex"] == 0
        and per_volume[air_volume]["wedge"] == 0,
        "global_topology_matches_region_sum": all(
            counts.get(kind, 0)
            == sum(row[kind] for row in per_volume.values())
            for kind in ("hex", "tet", "pyramid", "wedge")
        ),
        "quality_count_matches_topology": all(
            quality_counts[kind] == counts.get(kind, 0)
            for kind in ("hex", "tet", "pyramid")
        ),
        "all_volume_families_above_quality_threshold": all(
            isfinite(quality_minima[kind]) and quality_minima[kind] >= quality_limit
            for kind in ("hex", "tet", "pyramid")
        ),
        "one_conformal_quad_interface_per_conductor_region": interfaces_conformal
        and len(interface_rows) == len(conductor_volumes)
        and all(count == 1 for count in interface_by_conductor.values()),
        "cad_partition_volume_is_conserved": isfinite(volume_error)
        and volume_error <= volume_limit,
        "gmsh_ascii_v41_handoff_complete": int(export.get("bytes", 0)) > 0
        and len(str(export.get("sha256", ""))) == 64
        and str(header.get("version", "")) == "4.1"
        and int(header.get("file_type", -1)) == 0
        and all(
            header.get(name) is True
            for name in ("has_entities_section", "has_nodes_section", "has_elements_section")
        ),
        "parsed_gmsh_topology_matches_cubit": inventory.get("status") == "ok"
        and not list(inventory.get("connectivity_mismatches") or [])
        and all(
            inventory_counts.get(kind, 0) == counts.get(kind, 0)
            for kind in ("hex", "tet", "pyramid")
        )
        and inventory_counts.get("wedge", 0) == counts.get("wedge", 0),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_region_owned_mixed_mesh_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "element_counts": counts,
        "quality_minima": quality_minima,
        "interface_count": len(interface_rows),
        "interface_count_by_conductor": interface_by_conductor,
        "volume_relative_error": volume_error,
        "parsed_gmsh_volume_family_counts": inventory_counts,
        "notes": [
            "Judge mixed topology by region ownership, not by requiring the global hex count to exceed the air tet count.",
            "Mesh sweepable conductors first; tetmesh may then insert pyramids against their quad interfaces.",
            "A nonempty export is insufficient without per-region topology, quality, interface, and CAD-volume closure.",
            "Parse the Gmsh 4.1 element blocks independently and require their volume-family counts to match Cubit's in-memory inventory.",
        ],
    }


def cubit_helical_conductor_source_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate the source-journal replay and classified HPC-tet fallback."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    process = _mapping(summary.get("process") or {}, "process")
    timing = _mapping(summary.get("timing") or {}, "timing")
    parameters = _mapping(summary.get("geometry_parameters") or {}, "geometry_parameters")
    public_gate = cubit_region_owned_mixed_mesh_gate(summary)
    categories = {str(value) for value in process.get("error_categories") or []}
    unexpected = list(process.get("unexpected_error_lines") or [])
    exit_code = int(process.get("exit_code", -1))
    classified_nonzero = exit_code == 0 or (
        exit_code == 4
        and {
            "hpc_tet_attempt_failed",
            "standard_tet_pyramid_fallback_completed",
            "session_error_summary",
        }.issubset(categories)
        and process.get("hpc_fallback_completed") is True
        and process.get("process_exit_policy")
        == "classified_hpc_fallback_plus_fresh_artifact"
        and not unexpected
        and process.get("result_artifact_fresh") is True
        and public_gate["status"] == "ok"
    )
    timing_values = list(timing.values())
    source_sha = str(summary.get("source_sha256", "")).lower()
    checks = {
        "source_native_helical_journal_identified": str(
            summary.get("source_kind", "")
        ).startswith("source_native_local_helical_conductor_journal")
        and Path(str(summary.get("source_journal", ""))).suffix.lower() == ".jou",
        "source_digest_is_sha256": len(source_sha) == 64
        and all(character in "0123456789abcdef" for character in source_sha),
        "four_turn_geometry_parameters_recorded": int(parameters.get("turns", 0)) == 4
        and all(float(parameters.get(name, 0.0)) > 0.0 for name in ("a_mm", "b_mm", "h_mm")),
        "headless_combined_replay_recorded": summary.get("execution_mode")
        == "headless_combined_journal_then_python_inventory"
        and {"-nographics", "-batch"}.issubset(set(summary.get("headless_flags") or []))
        and summary.get("gui_daemon_enabled") is False,
        "nonzero_exit_has_classified_completed_fallback": classified_nonzero,
        "fresh_artifact_and_no_owned_process_leak": process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "exactly_four_source_timing_stages": len(timing_values) == 4
        and all(isfinite(float(value)) and float(value) >= 0.0 for value in timing_values),
        "independent_region_owned_mesh_gate_passed": public_gate["status"] == "ok",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_helical_conductor_source_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "process_exit_code": exit_code,
        "error_categories": sorted(categories),
        "public_gate_status": public_gate["status"],
        "notes": [
            "Do not allowlist Cubit exit code 4 by number alone.",
            "Require the failed HPC attempt, completed standard tet/pyramid fallback, fresh artifact, no unexpected errors, and an independent mesh gate.",
            "The source lesson is a replay contract; the generic region-owned mesh gate remains vendor-neutral engineering knowledge.",
        ],
    }
