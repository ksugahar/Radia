from __future__ import annotations

from radia_mcp.cubit.cross_artifact_mesh_lineage_v47 import validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v47_public_hex_tet_transition_interface_face_owner_conservation_mismatch",
    "v47_public_block_nodeset_sideset_remap_after_merge_generation_mismatch",
    "v47_source_tool_headless_journal_dependency_order_entity_generation_mismatch",
    "v47_source_tool_export_manifest_part_order_physical_group_mapping_mismatch",
}


def _records() -> dict[str, object]:
    interface_generation = "transition-interface-v47-901"
    remap_generation = "set-remap-v47-901"
    journal_generation = "journal-dependency-v47-901"
    export_generation = "export-manifest-v47-901"
    faces = ["face:101", "face:102", "face:103"]
    owners = [["hex:block1", "tet:block2"] for _ in faces]
    remap = {"blocks": {"1": "11", "2": "12"}, "nodesets": {"20": "120"}, "sidesets": {"30": "130"}}
    commands = ["reset", "create brick 1", "imprint all", "mesh volume all"]
    dependencies = [[], [0], [1], [2]]
    entities = {"volume:1": journal_generation, "surface:1": journal_generation}
    parts = ["rotor", "stator", "airgap"]
    groups = {"rotor": 101, "stator": 102, "airgap": 103}
    return {
        "hex_tet_transition_interface_face_owner_conservation_identity": {
            "generation": interface_generation, "interface_generation": interface_generation, "block_generation": interface_generation, "result_generation": interface_generation,
            "interface_face_ids": faces, "result_interface_face_ids": faces, "interface_owner_pairs": owners, "result_interface_owner_pairs": owners,
            "duplicate_interface_face_count": 0, "result_duplicate_interface_face_count": 0, "unowned_interface_face_count": 0, "result_unowned_interface_face_count": 0,
            "owner": "headless:transition-v47", "accepted_owner": "headless:transition-v47", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        "block_nodeset_sideset_remap_after_merge_generation_identity": {
            "generation": remap_generation, "merge_generation": remap_generation, "imprint_generation": remap_generation, "set_remap_generation": remap_generation, "result_generation": remap_generation,
            "entity_remap": remap, "result_entity_remap": remap, "orphan_set_count": 0, "result_orphan_set_count": 0, "duplicate_target_count": 0, "result_duplicate_target_count": 0,
            "owner": "headless:remap-v47", "accepted_owner": "headless:remap-v47", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        "headless_journal_dependency_order_entity_generation_identity": {
            "generation": journal_generation, "command_generation": journal_generation, "entity_generation": journal_generation, "result_generation": journal_generation,
            "commands": commands, "result_commands": commands, "dependency_order": dependencies, "result_dependency_order": dependencies, "entity_generations": entities, "result_entity_generations": entities,
            "stale_entity_reference_count": 0, "result_stale_entity_reference_count": 0, "owner": "headless:journal-v47", "accepted_owner": "headless:journal-v47", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "export_manifest_part_order_physical_group_mapping_identity": {
            "generation": export_generation, "manifest_generation": export_generation, "physical_group_generation": export_generation, "result_generation": export_generation,
            "part_order": parts, "result_part_order": parts, "physical_group_map": groups, "result_physical_group_map": groups,
            "duplicate_physical_group_count": 0, "result_duplicate_physical_group_count": 0, "owner": "headless:export-v47", "accepted_owner": "headless:export-v47", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }


def test_v47_positive_replays_are_accepted() -> None:
    records = _records()
    assert validate_public_identity(records)["status"] == "ok"
    assert validate_source_identity(records)["status"] == "ok"


def test_v47_interface_and_set_mutations_are_rejected() -> None:
    records = _records()
    records["hex_tet_transition_interface_face_owner_conservation_identity"]["result_interface_face_ids"] = ["face:101", "face:101", "face:103"]
    records["block_nodeset_sideset_remap_after_merge_generation_identity"]["set_remap_generation"] = "old"
    assert validate_public_identity(records)["status"] == "needs_attention"


def test_v47_journal_and_export_mutations_are_rejected() -> None:
    records = _records()
    records["headless_journal_dependency_order_entity_generation_identity"]["result_stale_entity_reference_count"] = 1
    records["export_manifest_part_order_physical_group_mapping_identity"]["result_part_order"] = ["airgap", "rotor", "stator"]
    assert validate_source_identity(records)["status"] == "needs_attention"
