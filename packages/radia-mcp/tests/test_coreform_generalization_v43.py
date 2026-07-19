from __future__ import annotations

from test_coreform_generalization_v37 import _public_result, _source_result, summary
from test_coreform_generalization_v39 import _generations


_HYBRID = (
    "hybrid_tet_hex_pyramid_transition_interface_orientation_quality_block_export_generation_identity"
)
_WEBCUT = (
    "webcut_multivolume_sharedface_block_connectivity_hexcount_euler_quality_export_generation_identity"
)
_SIDESET = (
    "sideset_normal_propagation_merge_mesh_export_entity_owner_result_generation_identity"
)
_CAD = (
    "cad_import_healing_tolerance_units_bodycount_topology_database_owner_result_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v43_public_hybrid_tet_hex_pyramid_transition_interface_orientation_quality_export_mismatch",
    "v43_public_webcut_multivolume_sharedface_block_connectivity_hexcount_euler_export_mismatch",
    "v43_source_sideset_normal_propagation_merge_mesh_export_entity_owner_mismatch",
    "v43_source_cad_import_healing_tolerance_units_bodycount_topology_database_owner_mismatch",
)


def _with_v43_coreform_identity(row: dict) -> dict:
    generation = "hybrid-transition-726"
    row[_HYBRID] = {
        "hybrid_generation": generation,
        **_generations(
            generation,
            "transition_generation",
            "interface_generation",
            "orientation_generation",
            "quality_generation",
            "block_generation",
            "export_generation",
            "result_generation",
        ),
        "hex_element_count": 64,
        "result_hex_element_count": 64,
        "tet_element_count": 128,
        "result_tet_element_count": 128,
        "pyramid_element_count": 24,
        "result_pyramid_element_count": 24,
        "hex_pyramid_interface_face_count": 24,
        "result_hex_pyramid_interface_face_count": 24,
        "pyramid_tet_interface_face_count": 96,
        "result_pyramid_tet_interface_face_count": 96,
        "interface_orientation_dot_products": [-1.0, -1.0, -1.0],
        "result_interface_orientation_dot_products": [-1.0, -1.0, -1.0],
        "minimum_hex_scaled_jacobian": 0.42,
        "result_minimum_hex_scaled_jacobian": 0.42,
        "minimum_pyramid_scaled_jacobian": 0.31,
        "result_minimum_pyramid_scaled_jacobian": 0.31,
        "minimum_tet_scaled_jacobian": 0.28,
        "result_minimum_tet_scaled_jacobian": 0.28,
        "minimum_allowed_scaled_jacobian": 0.20,
        "result_minimum_allowed_scaled_jacobian": 0.20,
        "block_membership": {"block:hex": [1], "block:transition": [2], "block:tet": [3]},
        "result_block_membership": {"block:hex": [1], "block:transition": [2], "block:tet": [3]},
        "mesh_owner": "headless:hybrid-726",
        "result_mesh_owner": "headless:hybrid-726",
        "hybrid_export_sha256": "1" * 64,
        "accepted_hybrid_export_sha256": "1" * 64,
    }

    generation = "webcut-blocks-726"
    row[_WEBCUT] = {
        "webcut_generation": generation,
        **_generations(
            generation,
            "topology_generation",
            "sharedface_generation",
            "connectivity_generation",
            "mesh_generation",
            "quality_generation",
            "block_generation",
            "export_generation",
            "result_generation",
        ),
        "volume_ids": [1, 2, 3],
        "result_volume_ids": [1, 2, 3],
        "shared_face_pairs": [[1, 2, 11], [2, 3, 12]],
        "result_shared_face_pairs": [[1, 2, 11], [2, 3, 12]],
        "block_connectivity": {"block:10": [1, 2], "block:20": [3]},
        "result_block_connectivity": {"block:10": [1, 2], "block:20": [3]},
        "hex_element_count": 384,
        "result_hex_element_count": 384,
        "topology_euler_characteristic": 1,
        "result_topology_euler_characteristic": 1,
        "minimum_scaled_jacobian": 0.37,
        "result_minimum_scaled_jacobian": 0.37,
        "minimum_allowed_scaled_jacobian": 0.20,
        "result_minimum_allowed_scaled_jacobian": 0.20,
        "mesh_owner": "headless:webcut-726",
        "result_mesh_owner": "headless:webcut-726",
        "webcut_export_sha256": "2" * 64,
        "accepted_webcut_export_sha256": "2" * 64,
    }

    generation = "sideset-propagation-726"
    row[_SIDESET] = {
        "sideset_generation": generation,
        **_generations(
            generation,
            "normal_generation",
            "merge_generation",
            "mesh_generation",
            "entity_generation",
            "block_generation",
            "database_generation",
            "export_generation",
            "result_generation",
        ),
        "sideset_membership": {"sideset:20": [11, 12]},
        "replay_sideset_membership": {"sideset:20": [11, 12]},
        "outward_normals": {"surface:11": [1.0, 0.0, 0.0], "surface:12": [0.0, 1.0, 0.0]},
        "replay_outward_normals": {"surface:11": [1.0, 0.0, 0.0], "surface:12": [0.0, 1.0, 0.0]},
        "merge_entity_map": {"surface:101": "surface:11", "surface:102": "surface:12"},
        "replay_merge_entity_map": {"surface:101": "surface:11", "surface:102": "surface:12"},
        "mesh_generation_id": 726,
        "replay_mesh_generation_id": 726,
        "block_membership": {"block:10": [1, 2]},
        "replay_block_membership": {"block:10": [1, 2]},
        "database_owner": "headless:sideset-726",
        "replay_database_owner": "headless:sideset-726",
        "mesh_export_sha256": "3" * 64,
        "replay_mesh_export_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }

    generation = "cad-healing-726"
    row[_CAD] = {
        "cad_generation": generation,
        **_generations(
            generation,
            "unit_generation",
            "healing_generation",
            "topology_generation",
            "entity_generation",
            "database_generation",
            "owner_generation",
            "result_generation",
        ),
        "source_unit": "mm",
        "replay_source_unit": "mm",
        "model_unit": "m",
        "replay_model_unit": "m",
        "unit_scale_to_m": 0.001,
        "replay_unit_scale_to_m": 0.001,
        "healing_tolerance_m": 1.0e-6,
        "replay_healing_tolerance_m": 1.0e-6,
        "body_count": 3,
        "replay_body_count": 3,
        "watertight_body_ids": [1, 2, 3],
        "replay_watertight_body_ids": [1, 2, 3],
        "topology_counts": {"volume": 3, "surface": 18, "curve": 36, "vertex": 24},
        "replay_topology_counts": {"volume": 3, "surface": 18, "curve": 36, "vertex": 24},
        "entity_generation_id": 726,
        "replay_entity_generation_id": 726,
        "database_owner": "headless:cad-726",
        "replay_database_owner": "headless:cad-726",
        "source_cad_sha256": "5" * 64,
        "replay_source_cad_sha256": "5" * 64,
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    return row


def test_v43_positive_public_and_source_contracts() -> None:
    row = _with_v43_coreform_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v43_rejects_hybrid_transition_mismatch() -> None:
    row = _with_v43_coreform_identity(summary())
    row[_HYBRID]["result_pyramid_element_count"] = 0
    row[_HYBRID]["result_mesh_owner"] = "gui:old"
    assert _public_result(row)["status"] == "needs_attention"


def test_v43_rejects_webcut_connectivity_mismatch() -> None:
    row = _with_v43_coreform_identity(summary())
    row[_WEBCUT]["result_shared_face_pairs"] = []
    row[_WEBCUT]["result_topology_euler_characteristic"] = 3
    assert _public_result(row)["status"] == "needs_attention"


def test_v43_rejects_sideset_normal_replay_mismatch() -> None:
    row = _with_v43_coreform_identity(summary())
    row[_SIDESET]["replay_outward_normals"] = {"surface:11": [-1.0, 0.0, 0.0]}
    row[_SIDESET]["replay_database_owner"] = "gui:old"
    assert _source_result(row)["status"] == "needs_attention"


def test_v43_rejects_cad_healing_units_and_topology_mismatch() -> None:
    row = _with_v43_coreform_identity(summary())
    row[_CAD]["replay_source_unit"] = "in"
    row[_CAD]["replay_unit_scale_to_m"] = 0.0254
    row[_CAD]["replay_body_count"] = 1
    row[_CAD]["replay_watertight_body_ids"] = []
    assert _source_result(row)["status"] == "needs_attention"


def test_v43_rejects_self_consistent_bad_pyramid_orientation() -> None:
    row = _with_v43_coreform_identity(summary())
    row[_HYBRID]["interface_orientation_dot_products"] = [-1.0, 1.0, -1.0]
    row[_HYBRID]["result_interface_orientation_dot_products"] = [-1.0, 1.0, -1.0]
    assert _public_result(row)["status"] == "needs_attention"


def test_v43_rejects_self_consistent_nonconformal_webcut() -> None:
    row = _with_v43_coreform_identity(summary())
    row[_WEBCUT]["shared_face_pairs"] = [[1, 2, 11]]
    row[_WEBCUT]["result_shared_face_pairs"] = [[1, 2, 11]]
    assert _public_result(row)["status"] == "needs_attention"
