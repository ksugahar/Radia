"""Cross-artifact mesh ownership and headless replay checks for v47."""

from __future__ import annotations

from collections.abc import Mapping


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _result_ok(row: Mapping[str, object]) -> bool:
    return (
        bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _interface_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    faces = row.get("interface_face_ids")
    owners = row.get("interface_owner_pairs")
    return (
        bool(generation)
        and row.get("interface_generation") == generation
        and row.get("block_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(faces, list)
        and bool(faces)
        and all(isinstance(face, str) and bool(face.strip()) for face in faces)
        and len(set(faces)) == len(faces)
        and row.get("result_interface_face_ids") == faces
        and isinstance(owners, list)
        and len(owners) == len(faces)
        and all(isinstance(pair, list) and len(pair) == 2 and all(isinstance(owner, str) and bool(owner.strip()) for owner in pair) and pair[0] != pair[1] for pair in owners)
        and row.get("result_interface_owner_pairs") == owners
        and row.get("duplicate_interface_face_count") == row.get("result_duplicate_interface_face_count") == 0
        and row.get("unowned_interface_face_count") == row.get("result_unowned_interface_face_count") == 0
        and _result_ok(row)
    )


def _remap_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    remap = row.get("entity_remap")
    required = {"blocks", "nodesets", "sidesets"}
    if not isinstance(remap, Mapping) or set(remap) != required:
        return False
    targets: list[str] = []
    for family in required:
        mapping = remap.get(family)
        if not isinstance(mapping, Mapping) or not mapping:
            return False
        targets.extend(str(value) for value in mapping.values())
    return (
        bool(generation)
        and row.get("merge_generation") == generation
        and row.get("imprint_generation") == generation
        and row.get("set_remap_generation") == generation
        and row.get("result_generation") == generation
        and len(set(targets)) == len(targets)
        and row.get("result_entity_remap") == remap
        and row.get("orphan_set_count") == row.get("result_orphan_set_count") == 0
        and row.get("duplicate_target_count") == row.get("result_duplicate_target_count") == 0
        and _result_ok(row)
    )


def _journal_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    commands = row.get("commands")
    dependencies = row.get("dependency_order")
    entities = row.get("entity_generations")
    valid_dependencies = (
        isinstance(commands, list)
        and isinstance(dependencies, list)
        and len(commands) == len(dependencies)
        and all(isinstance(deps, list) and all(isinstance(index, int) and 0 <= index < step for index in deps) for step, deps in enumerate(dependencies))
    )
    return (
        bool(generation)
        and row.get("command_generation") == generation
        and row.get("entity_generation") == generation
        and row.get("result_generation") == generation
        and valid_dependencies
        and row.get("result_commands") == commands
        and row.get("result_dependency_order") == dependencies
        and isinstance(entities, Mapping)
        and bool(entities)
        and all(value == generation for value in entities.values())
        and row.get("result_entity_generations") == entities
        and row.get("stale_entity_reference_count") == row.get("result_stale_entity_reference_count") == 0
        and _result_ok(row)
    )


def _export_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    parts = row.get("part_order")
    groups = row.get("physical_group_map")
    return (
        bool(generation)
        and row.get("manifest_generation") == generation
        and row.get("physical_group_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(parts, list)
        and bool(parts)
        and len(set(parts)) == len(parts)
        and isinstance(groups, Mapping)
        and set(groups) == set(parts)
        and all(isinstance(value, int) and value > 0 for value in groups.values())
        and len(set(groups.values())) == len(groups)
        and row.get("result_part_order") == parts
        and row.get("result_physical_group_map") == groups
        and row.get("duplicate_physical_group_count") == row.get("result_duplicate_physical_group_count") == 0
        and _result_ok(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    interface = payload.get("hex_tet_transition_interface_face_owner_conservation_identity")
    remap = payload.get("block_nodeset_sideset_remap_after_merge_generation_identity")
    if interface is not None:
        checks["v47_transition_interface_face_ownership"] = isinstance(interface, Mapping) and _interface_ok(interface)
    if remap is not None:
        checks["v47_block_nodeset_sideset_remap"] = isinstance(remap, Mapping) and _remap_ok(remap)
    return _report("cubit_v47_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    journal = payload.get("headless_journal_dependency_order_entity_generation_identity")
    export = payload.get("export_manifest_part_order_physical_group_mapping_identity")
    if journal is not None:
        checks["v47_headless_journal_dependency_generation"] = isinstance(journal, Mapping) and _journal_ok(journal)
    if export is not None:
        checks["v47_export_part_physical_group_mapping"] = isinstance(export, Mapping) and _export_ok(export)
    return _report("cubit_v47_source_identity_v1", checks) if checks else {}
