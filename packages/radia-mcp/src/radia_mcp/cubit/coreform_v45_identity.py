"""Identity gates for periodic high-order and headless mesh replays."""

from __future__ import annotations

from collections.abc import Mapping, Sequence


def _sha(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _same(identity: Mapping[str, object], *names: str) -> bool:
    return all(identity.get(f"result_{name}") == identity.get(name) for name in names)


def periodic_hex_v45_ok(summary: Mapping[str, object]) -> bool:
    identity = summary.get("hex_sweep_periodic_pairing_curvature_jacobian_block_export_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    source = identity.get("paired_source_node_ids")
    target = identity.get("paired_target_node_ids")
    result_source = identity.get("result_paired_source_node_ids")
    result_target = identity.get("result_paired_target_node_ids")
    try:
        curvature = int(identity["curvature_order"])
        jacobian = float(identity["result_minimum_scaled_jacobian"])
        threshold = float(identity["minimum_allowed_scaled_jacobian"])
    except (KeyError, TypeError, ValueError):
        return False
    generations = ("generation", "pairing_generation", "curvature_generation", "jacobian_generation", "block_generation", "export_generation")
    return (
        all(bool(str(identity.get(name) or "")) for name in generations)
        and isinstance(source, Sequence) and not isinstance(source, (str, bytes))
        and isinstance(target, Sequence) and not isinstance(target, (str, bytes))
        and list(source) == list(result_source or [])
        and list(target) == list(result_target or [])
        and len(source) == len(target) >= 2
        and curvature >= 2
        and jacobian >= threshold > 0.0
        and _same(identity, "source_surface_id", "target_surface_id", "curvature_order", "minimum_scaled_jacobian", "block_membership", "mesh_owner")
        and bool(str(identity.get("mesh_owner") or "").startswith("headless:"))
        and _sha(identity.get("mesh_export_sha256"))
        and identity.get("accepted_mesh_export_sha256") == identity.get("mesh_export_sha256")
    )


def journal_v45_ok(summary: Mapping[str, object]) -> bool:
    identity = summary.get("headless_journal_units_command_status_geometry_mesh_database_owner_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    commands = identity.get("commands")
    replay = identity.get("replay_commands")
    statuses = identity.get("command_status")
    replay_statuses = identity.get("replay_command_status")
    return (
        bool(str(identity.get("generation") or ""))
        and isinstance(commands, Sequence) and not isinstance(commands, (str, bytes))
        and list(commands) == list(replay or [])
        and isinstance(statuses, Sequence) and list(statuses) == list(replay_statuses or [])
        and bool(statuses) and all(str(status) == "success" for status in statuses)
        and identity.get("session_units") == identity.get("replay_session_units") == "mm"
        and identity.get("geometry_generation_id") == identity.get("replay_geometry_generation_id")
        and identity.get("mesh_generation_id") == identity.get("replay_mesh_generation_id")
        and identity.get("database_owner") == identity.get("replay_database_owner")
        and bool(str(identity.get("database_owner") or "").startswith("headless:"))
        and _sha(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def quality_v45_ok(summary: Mapping[str, object]) -> bool:
    identity = summary.get("quality_metric_reference_element_dimension_threshold_block_export_owner_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        dimension = int(identity["dimension"])
        replay_dimension = int(identity["replay_dimension"])
        threshold = float(identity["quality_threshold"])
        replay_threshold = float(identity["replay_quality_threshold"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(str(identity.get("generation") or ""))
        and identity.get("reference_element") == identity.get("replay_reference_element") in {"hex8", "hex27"}
        and dimension == replay_dimension == 3
        and threshold == replay_threshold > 0.0
        and identity.get("block_name") == identity.get("replay_block_name")
        and identity.get("export_generation_id") == identity.get("replay_export_generation_id")
        and identity.get("database_owner") == identity.get("replay_database_owner")
        and bool(str(identity.get("database_owner") or "").startswith("headless:"))
        and _sha(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def coreform_v45_ok(summary: Mapping[str, object]) -> bool:
    return periodic_hex_v45_ok(summary) and journal_v45_ok(summary) and quality_v45_ok(summary)
