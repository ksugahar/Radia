"""Hex topology, mesh-group, and source-replay identity checks for v50."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_HEX = "hex_node_order_face_adjacency_orientation_jacobian_mesh_owner_identity"
_GROUP = "block_sideset_entity_dimension_overlap_duplicate_membership_owner_identity"
_JOURNAL = "journal_include_relative_path_aprepro_scope_expansion_digest_owner_identity"
_RESTORE = "save_restore_session_entity_id_block_sideset_mesh_revision_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _positive_unique_ints(value: object, count: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (count is None or len(value) == count)
        and bool(value)
        and all(isinstance(item, int) and not isinstance(item, bool) and item > 0 for item in value)
        and len(set(value)) == len(value)
    )


def _hex_ok(row: Mapping[str, object]) -> bool:
    node_order = row.get("hex_node_order")
    adjacency = row.get("face_adjacency")
    orientations = row.get("face_orientation_signs")
    jacobian = row.get("minimum_scaled_jacobian")
    return (
        _generations(
            row,
            "node_order_generation",
            "face_adjacency_generation",
            "orientation_generation",
            "jacobian_generation",
            "owner_generation",
            "result_generation",
        )
        and _positive_unique_ints(node_order, 8)
        and row.get("result_hex_node_order") == node_order
        and isinstance(adjacency, Mapping)
        and bool(adjacency)
        and all(
            isinstance(face, str)
            and face.startswith("face:")
            and isinstance(neighbor, str)
            and neighbor.startswith(("boundary:", "hex:"))
            for face, neighbor in adjacency.items()
        )
        and row.get("result_face_adjacency") == adjacency
        and isinstance(orientations, Sequence)
        and not isinstance(orientations, (str, bytes))
        and len(orientations) == 6
        and all(value in {-1, 1} for value in orientations)
        and row.get("result_face_orientation_signs") == orientations
        and isinstance(jacobian, (int, float))
        and math.isfinite(float(jacobian))
        and float(jacobian) > 0.0
        and row.get("result_minimum_scaled_jacobian") == jacobian
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _result(row)
    )


def _membership_map(value: object, keys: set[str]) -> bool:
    return (
        isinstance(value, Mapping)
        and set(value) == keys
        and all(_positive_unique_ints(members) for members in value.values())
    )


def _group_ok(row: Mapping[str, object]) -> bool:
    dimensions = row.get("entity_dimensions")
    memberships = row.get("group_memberships")
    overlaps = row.get("allowed_overlaps")
    duplicates = row.get("duplicate_memberships")
    if not isinstance(dimensions, Mapping) or not dimensions:
        return False
    keys = set(dimensions)
    expected_dimensions = all(
        isinstance(name, str)
        and ((name.startswith("block:") and dimension == 3) or (name.startswith("sideset:") and dimension == 2))
        for name, dimension in dimensions.items()
    )
    if not expected_dimensions or not _membership_map(memberships, keys):
        return False
    overlap_set = set(overlaps) if isinstance(overlaps, list) else set()
    observed_overlap: set[int] = set()
    groups = list(memberships.items())
    for index, (left_name, left_members) in enumerate(groups):
        for right_name, right_members in groups[index + 1 :]:
            if left_name.split(":", 1)[0] == right_name.split(":", 1)[0]:
                observed_overlap.update(set(left_members) & set(right_members))
    return (
        _generations(
            row,
            "dimension_generation",
            "membership_generation",
            "overlap_generation",
            "duplicate_generation",
            "owner_generation",
            "result_generation",
        )
        and row.get("result_entity_dimensions") == dimensions
        and row.get("result_group_memberships") == memberships
        and isinstance(overlaps, list)
        and all(isinstance(value, int) and value > 0 for value in overlaps)
        and len(overlaps) == len(overlap_set)
        and observed_overlap <= overlap_set
        and row.get("result_allowed_overlaps") == overlaps
        and duplicates == []
        and row.get("result_duplicate_memberships") == duplicates
        and str(row.get("group_owner") or "").startswith("headless:")
        and row.get("result_group_owner") == row.get("group_owner")
        and _result(row)
    )


def _safe_relative_paths(value: object) -> bool:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        return False
    normalized: list[str] = []
    for item in value:
        if not isinstance(item, str) or not item:
            return False
        path = item.replace("\\", "/")
        parts = path.split("/")
        if path.startswith("/") or re.match(r"^[A-Za-z]:", path) or any(part in {"", ".", ".."} for part in parts):
            return False
        normalized.append(path)
    return len(normalized) == len(set(normalized))


def _string_mapping(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(isinstance(key, str) and key and isinstance(item, str) and item for key, item in value.items())
    )


def _string_sequence(value: object) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and bool(value)
        and all(isinstance(item, str) and item for item in value)
    )


def _journal_ok(row: Mapping[str, object]) -> bool:
    paths = row.get("relative_include_paths")
    scope = row.get("aprepro_scope")
    commands = row.get("expanded_commands")
    return (
        _generations(
            row,
            "include_generation",
            "scope_generation",
            "expansion_generation",
            "digest_generation",
            "owner_generation",
            "result_generation",
        )
        and _safe_relative_paths(paths)
        and row.get("result_relative_include_paths") == paths
        and _string_mapping(scope)
        and row.get("result_aprepro_scope") == scope
        and _string_sequence(commands)
        and row.get("result_expanded_commands") == commands
        and _digest(row.get("expanded_journal_sha256"))
        and row.get("result_expanded_journal_sha256") == row.get("expanded_journal_sha256")
        and str(row.get("journal_owner") or "").startswith("headless:")
        and row.get("result_journal_owner") == row.get("journal_owner")
        and _result(row)
    )


def _group_entity_map(value: object, prefix: str, valid_ids: set[int]) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(
            isinstance(name, str)
            and name.startswith(prefix)
            and _positive_unique_ints(members)
            and set(members) <= valid_ids
            for name, members in value.items()
        )
    )


def _restore_ok(row: Mapping[str, object]) -> bool:
    entity_ids = row.get("entity_ids")
    if not (
        isinstance(entity_ids, Mapping)
        and bool(entity_ids)
        and all(isinstance(name, str) and name and isinstance(identifier, int) and identifier > 0 for name, identifier in entity_ids.items())
        and len(set(entity_ids.values())) == len(entity_ids)
    ):
        return False
    valid_ids = set(entity_ids.values())
    blocks = row.get("blocks")
    sidesets = row.get("sidesets")
    return (
        _generations(
            row,
            "session_generation",
            "entity_generation",
            "block_generation",
            "sideset_generation",
            "mesh_generation",
            "owner_generation",
            "result_generation",
        )
        and str(row.get("session_id") or "").startswith("session:")
        and row.get("result_session_id") == row.get("session_id")
        and row.get("result_entity_ids") == entity_ids
        and _group_entity_map(blocks, "block:", valid_ids)
        and row.get("result_blocks") == blocks
        and _group_entity_map(sidesets, "sideset:", valid_ids)
        and row.get("result_sidesets") == sidesets
        and str(row.get("mesh_revision") or "").startswith("mesh-revision:")
        and row.get("result_mesh_revision") == row.get("mesh_revision")
        and str(row.get("model_owner") or "").startswith("headless:")
        and row.get("result_model_owner") == row.get("model_owner")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    hex_identity = payload.get(_HEX)
    group_identity = payload.get(_GROUP)
    if hex_identity is not None:
        checks["v50_hex_node_order_adjacency_orientation_jacobian_owner"] = (
            isinstance(hex_identity, Mapping) and _hex_ok(hex_identity)
        )
    if group_identity is not None:
        checks["v50_block_sideset_dimension_overlap_duplicate_owner"] = (
            isinstance(group_identity, Mapping) and _group_ok(group_identity)
        )
    return _report("cubit_v50_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    journal = payload.get(_JOURNAL)
    restore = payload.get(_RESTORE)
    if journal is not None:
        checks["v50_journal_include_aprepro_expansion_digest_owner"] = (
            isinstance(journal, Mapping) and _journal_ok(journal)
        )
    if restore is not None:
        checks["v50_save_restore_session_entities_groups_mesh_owner"] = (
            isinstance(restore, Mapping) and _restore_ok(restore)
        )
    return _report("cubit_v50_source_identity_v1", checks) if checks else {}
