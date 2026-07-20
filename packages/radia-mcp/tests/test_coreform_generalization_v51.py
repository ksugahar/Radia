from copy import deepcopy

from radia_mcp.cubit.quality_parallel_identity_v51 import validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v51_public_hex_quality_metric_reference_jacobian_sample_points_order_mesh_owner_mismatch",
    "v51_public_periodic_highorder_node_edge_face_parametric_rotation_owner_mismatch",
    "v51_source_tool_parallel_sculpt_partition_ghost_overlap_seed_merge_revision_owner_mismatch",
    "v51_source_tool_exodus_int64_idmap_names_qa_time_global_order_owner_mismatch",
}


def _records() -> dict[str, object]:
    quality = "hex-quality-v51-2001"
    periodic = "periodic-v51-2001"
    sculpt = "sculpt-v51-2001"
    exodus = "exodus-v51-2001"
    samples = [[-1.0, -1.0, -1.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
    node_pairs = [[101, 201], [102, 202], [103, 203]]
    edge = {"edge:11": [0.0, 0.5, 1.0]}
    face = {"face:21": [[0.0, 0.0], [0.5, 0.5], [1.0, 1.0]]}
    rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    partitions = {"rank:0": [1, 2], "rank:1": [3, 4]}
    overlap = [[2, 3]]
    ids = {"node:1": 4294967301, "element:1": 4294967401}
    names = {"block:1": "steel", "sideset:1": "wall"}
    qa = [["CAE-AI Lab", "export", "2026-07-20", "v51"]]
    times = [0.0, 0.5, 1.0]
    globals_ = ["energy", "volume"]
    return {
        "hex_quality_metric_reference_jacobian_samples_order_mesh_owner_identity": {"generation": quality, **{name: quality for name in ("metric_generation", "jacobian_generation", "sample_generation", "order_generation", "owner_generation", "result_generation")}, "quality_metric": "scaled_jacobian", "result_quality_metric": "scaled_jacobian", "reference_jacobian_sha256": "1" * 64, "result_reference_jacobian_sha256": "1" * 64, "sample_points": samples, "result_sample_points": samples, "element_order": 1, "result_element_order": 1, "mesh_owner": "headless:quality-v51", "result_mesh_owner": "headless:quality-v51", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
        "periodic_highorder_node_edge_face_parametric_rotation_owner_identity": {"generation": periodic, **{name: periodic for name in ("node_generation", "edge_generation", "face_generation", "rotation_generation", "owner_generation", "result_generation")}, "node_pairs": node_pairs, "result_node_pairs": node_pairs, "edge_parametric_coordinates": edge, "result_edge_parametric_coordinates": edge, "face_parametric_coordinates": face, "result_face_parametric_coordinates": face, "rotation_map": rotation, "result_rotation_map": rotation, "periodic_owner": "headless:periodic-v51", "result_periodic_owner": "headless:periodic-v51", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64},
        "parallel_sculpt_partition_ghost_seed_merge_revision_owner_identity": {"generation": sculpt, **{name: sculpt for name in ("partition_generation", "ghost_generation", "seed_generation", "merge_generation", "revision_generation", "owner_generation", "result_generation")}, "partitions": partitions, "result_partitions": partitions, "ghost_overlap": overlap, "result_ghost_overlap": overlap, "random_seed": 51001, "result_random_seed": 51001, "merge_order": ["rank:0", "rank:1"], "result_merge_order": ["rank:0", "rank:1"], "mesh_revision": "mesh-revision:v51-r8", "result_mesh_revision": "mesh-revision:v51-r8", "job_owner": "headless:sculpt-v51", "result_job_owner": "headless:sculpt-v51", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64},
        "exodus_int64_idmap_names_qa_time_global_order_owner_identity": {"generation": exodus, **{name: exodus for name in ("id_generation", "name_generation", "qa_generation", "time_generation", "global_generation", "owner_generation", "result_generation")}, "int64_ids": ids, "result_int64_ids": ids, "entity_names": names, "result_entity_names": names, "qa_records": qa, "result_qa_records": qa, "time_values": times, "result_time_values": times, "global_variable_order": globals_, "result_global_variable_order": globals_, "database_owner": "headless:exodus-v51", "result_database_owner": "headless:exodus-v51", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64},
    }


def test_v51_positive_public_and_source_replays_are_accepted() -> None:
    assert validate_public_identity(_records())["status"] == "ok"
    assert validate_source_identity(_records())["status"] == "ok"


def test_v51_public_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    records["hex_quality_metric_reference_jacobian_samples_order_mesh_owner_identity"]["result_quality_metric"] = "condition_number"
    records["periodic_highorder_node_edge_face_parametric_rotation_owner_identity"]["result_node_pairs"] = [[101, 202], [102, 201]]
    assert validate_public_identity(records)["status"] == "needs_attention"


def test_v51_source_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    records["parallel_sculpt_partition_ghost_seed_merge_revision_owner_identity"]["result_random_seed"] = 51002
    records["exodus_int64_idmap_names_qa_time_global_order_owner_identity"]["result_global_variable_order"] = ["volume", "energy"]
    assert validate_source_identity(records)["status"] == "needs_attention"


def test_v51_invalid_canonical_records_are_rejected() -> None:
    records = deepcopy(_records())
    records["periodic_highorder_node_edge_face_parametric_rotation_owner_identity"]["rotation_map"] = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]
    records["parallel_sculpt_partition_ghost_seed_merge_revision_owner_identity"]["partitions"] = {"rank:0": [1, 2], "rank:1": [2, 3]}
    records["exodus_int64_idmap_names_qa_time_global_order_owner_identity"]["int64_ids"] = {"node:1": 101}
    assert validate_public_identity(records)["status"] == "needs_attention"
    assert validate_source_identity(records)["status"] == "needs_attention"
