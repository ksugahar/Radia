from __future__ import annotations

from copy import deepcopy

from radia_mcp.cubit.semantic_mesh_identity_v48 import validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v48_public_periodic_node_pair_orientation_translation_tolerance_partition_owner_mismatch",
    "v48_public_high_order_curved_edge_face_control_point_cad_signature_mismatch",
    "v48_source_tool_mesh_scheme_fallback_provenance_sweep_source_target_map_owner_mismatch",
    "v48_source_tool_exodus_timestep_global_variable_qa_record_mesh_revision_mismatch",
}


def _result(owner: str, digest: str) -> dict[str, object]:
    return {
        "owner": owner,
        "accepted_owner": owner,
        "result_sha256": digest * 64,
        "accepted_result_sha256": digest * 64,
    }


def _records() -> dict[str, object]:
    periodic_generation = "periodic-interface-v48-901"
    curved_generation = "curved-control-points-v48-901"
    scheme_generation = "mesh-scheme-fallback-v48-901"
    exodus_generation = "exodus-transient-v48-901"
    pairs = [[101, 201], [102, 202], [103, 203], [104, 204]]
    fallback = {"requested": "sweep", "applied": "tetmesh", "reason": "source_target_topology_incompatible"}
    sweep_map = {"volume:1": ["surface:11", "surface:12"]}
    times = [0.0, 0.1, 0.2]
    names = ["energy", "work"]
    rows = [[1.0, 0.0], [0.95, 0.05], [0.90, 0.10]]
    qa = [["cae-ai-lab", "cubit-export", "2026-07-20", "00:00:00"]]
    return {
        "periodic_node_pair_orientation_translation_tolerance_partition_owner_identity": {
            "generation": periodic_generation,
            "pair_generation": periodic_generation,
            "orientation_generation": periodic_generation,
            "translation_generation": periodic_generation,
            "partition_generation": periodic_generation,
            "result_generation": periodic_generation,
            "node_pairs": pairs,
            "result_node_pairs": pairs,
            "face_orientation": [1, -1],
            "result_face_orientation": [1, -1],
            "translation": [0.025, 0.0, 0.0],
            "result_translation": [0.025, 0.0, 0.0],
            "pair_tolerance": 1.0e-9,
            "maximum_pair_error": 2.0e-12,
            "result_maximum_pair_error": 2.0e-12,
            "partition_owner": "headless:periodic-partition-v48-901",
            "result_partition_owner": "headless:periodic-partition-v48-901",
            **_result("headless:periodic-interface-v48-901", "5"),
        },
        "high_order_curved_edge_face_control_point_cad_signature_identity": {
            "generation": curved_generation,
            "edge_generation": curved_generation,
            "face_generation": curved_generation,
            "cad_generation": curved_generation,
            "mesh_generation": curved_generation,
            "result_generation": curved_generation,
            "element_order": 3,
            "result_element_order": 3,
            "edge_control_point_order": {"edge:11": [1101, 1102], "edge:12": [1201, 1202]},
            "result_edge_control_point_order": {"edge:11": [1101, 1102], "edge:12": [1201, 1202]},
            "face_control_point_order": {"face:21": [2101, 2102, 2103, 2104]},
            "result_face_control_point_order": {"face:21": [2101, 2102, 2103, 2104]},
            "cad_geometry_signature": "6" * 64,
            "result_cad_geometry_signature": "6" * 64,
            "maximum_projection_error": 4.0e-10,
            "projection_tolerance": 1.0e-8,
            "result_maximum_projection_error": 4.0e-10,
            **_result("headless:curved-control-points-v48-901", "7"),
        },
        "mesh_scheme_fallback_provenance_sweep_source_target_map_owner_identity": {
            "generation": scheme_generation,
            "scheme_generation": scheme_generation,
            "fallback_generation": scheme_generation,
            "sweep_map_generation": scheme_generation,
            "volume_generation": scheme_generation,
            "result_generation": scheme_generation,
            "scheme_fallback": fallback,
            "result_scheme_fallback": fallback,
            "sweep_source_target_map": sweep_map,
            "result_sweep_source_target_map": sweep_map,
            "fallback_count": 1,
            "result_fallback_count": 1,
            **_result("headless:mesh-scheme-v48-901", "8"),
        },
        "exodus_timestep_global_variable_qa_record_mesh_revision_identity": {
            "generation": exodus_generation,
            "timestep_generation": exodus_generation,
            "global_variable_generation": exodus_generation,
            "qa_generation": exodus_generation,
            "mesh_generation": exodus_generation,
            "result_generation": exodus_generation,
            "timesteps": times,
            "result_timesteps": times,
            "global_variable_names": names,
            "result_global_variable_names": names,
            "global_variable_rows": rows,
            "result_global_variable_rows": rows,
            "qa_records": qa,
            "result_qa_records": qa,
            "mesh_revision": "mesh-revision-v48-901",
            "result_mesh_revision": "mesh-revision-v48-901",
            **_result("headless:exodus-transient-v48-901", "9"),
        },
    }


def test_v48_positive_public_and_source_replays_are_accepted() -> None:
    records = _records()
    assert validate_public_identity(records)["status"] == "ok"
    assert validate_source_identity(records)["status"] == "ok"


def test_v48_periodic_and_curved_semantic_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    periodic = records["periodic_node_pair_orientation_translation_tolerance_partition_owner_identity"]
    periodic["result_node_pairs"] = [[101, 201], [102, 202], [103, 204], [104, 203]]
    periodic["result_face_orientation"] = [1, 1]
    periodic["result_translation"] = [0.024, 0.0, 0.0]
    periodic["result_maximum_pair_error"] = 2.0e-6
    periodic["result_partition_owner"] = "headless:periodic-partition-v48-old"
    curved = records["high_order_curved_edge_face_control_point_cad_signature_identity"]
    curved["result_edge_control_point_order"] = {"edge:11": [1102, 1101], "edge:12": [1201, 1202]}
    curved["result_face_control_point_order"] = {"face:21": [2101, 2103, 2102, 2104]}
    curved["result_cad_geometry_signature"] = "a" * 64
    result = validate_public_identity(records)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) == {
        "v48_periodic_node_pair_transform_partition_owner",
        "v48_curved_control_point_cad_signature",
    }


def test_v48_scheme_and_exodus_semantic_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    scheme = records["mesh_scheme_fallback_provenance_sweep_source_target_map_owner_identity"]
    scheme["fallback_generation"] = "mesh-scheme-fallback-v48-old"
    scheme["result_sweep_source_target_map"] = {"volume:1": ["surface:12", "surface:11"]}
    scheme["accepted_owner"] = "headless:mesh-scheme-v48-old"
    exodus = records["exodus_timestep_global_variable_qa_record_mesh_revision_identity"]
    exodus["result_timesteps"] = [0.0, 0.2, 0.1]
    exodus["result_global_variable_rows"] = [[1.0, 0.0], [0.90, 0.10], [0.95, 0.05]]
    exodus["result_qa_records"] = [["cae-ai-lab", "old-export", "2026-07-19", "23:59:59"]]
    exodus["result_mesh_revision"] = "mesh-revision-v48-old"
    result = validate_source_identity(records)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) == {
        "v48_mesh_scheme_fallback_sweep_map_owner",
        "v48_exodus_timestep_global_qa_mesh_revision",
    }
