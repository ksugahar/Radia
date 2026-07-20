from copy import deepcopy

from radia_mcp.cubit.volume_frame_identity_v56 import BATCH, CURVED, EXODUS, VOLUME, validate_public_identity, validate_source_identity


CASE_IDS = {
    "v56_public_hexmesh_blockvolume_sum_cadvolume_orientation_owner_mismatch",
    "v56_public_curvedhex_faceconformity_nodeorder_jacobian_geometry_owner_mismatch",
    "v56_source_tool_batchjournal_errorstatus_rollback_outputdatabase_owner_mismatch",
    "v56_source_tool_exodus_coordinateframe_map_qarecord_timestep_owner_mismatch",
}


def _records() -> dict[str, object]:
    generation = "coreform-v56-test"
    generations = lambda fields: {field: generation for field in fields}
    volumes = {"block:1": [0.25, 0.25], "block:2": [0.2, 0.3]}
    sums = {"block:1": 0.5, "block:2": 0.5}
    faces = [{"face": "face:12", "side_a_nodes": [1, 2, 3, 4, 9, 10, 11, 12, 21], "side_b_nodes": [4, 3, 2, 1, 11, 10, 9, 12, 21]}]
    transform = [[1.0, 0.0, 0.0, 0.1], [0.0, 1.0, 0.0, 0.2], [0.0, 0.0, 1.0, 0.3], [0.0, 0.0, 0.0, 1.0]]
    qa = [{"program": "headless-export", "version": "v56", "date": "2026-07-20", "time": "09:00:00"}]
    steps = [{"index": 0, "time_s": 0.0}, {"index": 1, "time_s": 0.5}, {"index": 2, "time_s": 1.0}]
    return {
        VOLUME: {"generation": generation, **generations(("element_generation", "block_generation", "cad_generation", "orientation_generation", "owner_generation", "result_generation")), "signed_element_volume_m3": volumes, "result_signed_element_volume_m3": volumes, "block_volume_sum_m3": sums, "result_block_volume_sum_m3": sums, "cad_volume_m3": 1.0, "result_cad_volume_m3": 1.0, "orientation": "positive", "result_orientation": "positive", "mesh_owner": "headless:volume-v56", "result_mesh_owner": "headless:volume-v56", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
        CURVED: {"generation": generation, **generations(("face_generation", "order_generation", "jacobian_generation", "geometry_generation", "owner_generation", "result_generation")), "shared_face_nodes": faces, "result_shared_face_nodes": faces, "element_order": "HEX27", "result_element_order": "HEX27", "high_order_jacobian_samples": [0.42, 0.51], "result_high_order_jacobian_samples": [0.42, 0.51], "geometry_revision_sha256": "3" * 64, "result_geometry_revision_sha256": "3" * 64, "mesh_owner": "headless:curved-v56", "result_mesh_owner": "headless:curved-v56", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64},
        BATCH: {"generation": generation, **generations(("status_generation", "rollback_generation", "database_generation", "owner_generation", "result_generation")), "error_status": "none", "replayed_error_status": "none", "rollback_checkpoint": "checkpoint:pre-v56", "replayed_rollback_checkpoint": "checkpoint:pre-v56", "rollback_applied": False, "replayed_rollback_applied": False, "output_database_revision": "database:post-v56", "replayed_output_database_revision": "database:post-v56", "session_owner": "headless:batch-v56", "replayed_session_owner": "headless:batch-v56", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64},
        EXODUS: {"generation": generation, **generations(("frame_generation", "qa_generation", "timestep_generation", "owner_generation", "result_generation")), "coordinate_transform_4x4": transform, "replayed_coordinate_transform_4x4": transform, "qa_records": qa, "replayed_qa_records": qa, "time_steps": steps, "replayed_time_steps": steps, "file_owner": "headless:exodus-v56", "replayed_file_owner": "headless:exodus-v56", "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64},
    }


def test_v56_positive_identities_are_accepted() -> None:
    records = _records()
    assert validate_public_identity(records)["status"] == "ok"
    assert validate_source_identity(records)["status"] == "ok"


def test_v56_frozen_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    records[VOLUME]["result_cad_volume_m3"] = 2.0
    records[CURVED]["result_element_order"] = "HEX8"
    records[BATCH]["replayed_error_status"] = "failed"
    records[EXODUS]["replayed_time_steps"] = []
    assert validate_public_identity(records)["status"] == "needs_attention"
    assert validate_source_identity(records)["status"] == "needs_attention"


def test_v56_self_consistent_geometry_and_transaction_contradictions_are_rejected() -> None:
    records = deepcopy(_records())
    records[VOLUME]["orientation"] = records[VOLUME]["result_orientation"] = "negative"
    records[CURVED]["high_order_jacobian_samples"] = records[CURVED]["result_high_order_jacobian_samples"] = [-0.1]
    records[BATCH]["rollback_applied"] = records[BATCH]["replayed_rollback_applied"] = True
    records[EXODUS]["time_steps"] = records[EXODUS]["replayed_time_steps"] = [{"index": 1, "time_s": 0.0}]
    assert validate_public_identity(records)["status"] == "needs_attention"
    assert validate_source_identity(records)["status"] == "needs_attention"


def test_v56_malformed_values_reject_without_raising() -> None:
    records = deepcopy(_records())
    records[VOLUME]["signed_element_volume_m3"] = {"block:1": [[0.5]]}
    records[EXODUS]["coordinate_transform_4x4"] = [[1.0]]
    assert validate_public_identity(records)["status"] == "needs_attention"
    assert validate_source_identity(records)["status"] == "needs_attention"
