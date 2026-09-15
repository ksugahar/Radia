from __future__ import annotations

from test_coreform_generalization_v37 import _public_result, _source_result, summary
from test_coreform_generalization_v39 import _generations
from test_coreform_generalization_v41 import _with_v41_coreform_identity


_PERIODIC = "periodic_hex_surface_pair_transform_nodemap_orientation_interval_quality_block_export_generation_identity"
_BOUNDARY = "boundarylayer_hex_thickness_growth_layers_transition_normal_quality_block_export_generation_identity"
_ROLLBACK = "python_batch_exception_rollback_entity_generation_database_session_commandlog_result_generation_identity"
_RECIPE = "meshrecipe_dependency_dag_parameter_execution_scheme_set_export_owner_result_generation_identity"
_PROMOTED_CASE_IDS = (
    "v42_public_periodic_hex_surface_pair_transform_node_map_orientation_quality_export_mismatch",
    "v42_public_boundarylayer_hex_thickness_growth_layers_transition_quality_block_export_mismatch",
    "v42_source_python_batch_exception_rollback_entity_generation_database_session_owner_mismatch",
    "v42_source_meshrecipe_dependency_dag_parameter_hash_execution_order_export_owner_mismatch",
)


def _with_v42_coreform_identity(row: dict) -> dict:
    row = _with_v41_coreform_identity(row)
    generation = "periodic-hex-725"
    transform = [[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    node_map = [[101, 201], [102, 202], [103, 203], [104, 204]]
    row[_PERIODIC] = {
        "periodic_generation": generation,
        **_generations(generation, "surface_generation", "transform_generation", "nodemap_generation", "orientation_generation", "interval_generation", "quality_generation", "block_generation", "export_generation", "result_generation"),
        "source_surface_id": 11, "result_source_surface_id": 11,
        "target_surface_id": 21, "result_target_surface_id": 21,
        "rigid_transform": transform, "result_rigid_transform": transform,
        "node_map": node_map, "result_node_map": node_map,
        "source_face_normal": [1.0, 0.0, 0.0], "result_source_face_normal": [1.0, 0.0, 0.0],
        "target_face_normal": [-1.0, 0.0, 0.0], "result_target_face_normal": [-1.0, 0.0, 0.0],
        "source_interval_counts": [8, 6], "result_source_interval_counts": [8, 6],
        "target_interval_counts": [8, 6], "result_target_interval_counts": [8, 6],
        "minimum_scaled_jacobian": 0.38, "result_minimum_scaled_jacobian": 0.38,
        "minimum_allowed_scaled_jacobian": 0.20, "result_minimum_allowed_scaled_jacobian": 0.20,
        "block_owner": "block:periodic-725", "result_block_owner": "block:periodic-725",
        "periodic_export_sha256": "1" * 64, "accepted_periodic_export_sha256": "1" * 64,
    }

    generation = "boundary-layer-hex-725"
    first, growth, count = 1.0e-3, 1.2, 4
    total = first * (growth**count - 1.0) / (growth - 1.0)
    row[_BOUNDARY] = {
        "boundarylayer_generation": generation,
        **_generations(generation, "thickness_generation", "growth_generation", "layer_generation", "transition_generation", "normal_generation", "quality_generation", "block_generation", "export_generation", "result_generation"),
        "first_layer_thickness_m": first, "result_first_layer_thickness_m": first,
        "growth_ratio": growth, "result_growth_ratio": growth,
        "layer_count": count, "result_layer_count": count,
        "total_layer_thickness_m": total, "result_total_layer_thickness_m": total,
        "transition_topology": ["hex"] * count, "result_transition_topology": ["hex"] * count,
        "wall_normal": [0.0, 1.0, 0.0], "result_wall_normal": [0.0, 1.0, 0.0],
        "layer_direction": [0.0, 1.0, 0.0], "result_layer_direction": [0.0, 1.0, 0.0],
        "minimum_scaled_jacobian": 0.34, "result_minimum_scaled_jacobian": 0.34,
        "minimum_allowed_scaled_jacobian": 0.20, "result_minimum_allowed_scaled_jacobian": 0.20,
        "block_owner": "block:boundary-layer-725", "result_block_owner": "block:boundary-layer-725",
        "boundarylayer_export_sha256": "2" * 64, "accepted_boundarylayer_export_sha256": "2" * 64,
    }

    generation = "python-rollback-725"
    commands = ["reset", "create brick x 1 y 1 z 1", "mesh volume 1"]
    entities = {"volume:1": 724}
    row[_ROLLBACK] = {
        "batch_generation": generation,
        **_generations(generation, "exception_generation", "rollback_generation", "entity_generation", "model_generation", "database_generation", "session_generation", "commandlog_generation", "result_generation"),
        "commands": commands, "replay_commands": commands,
        "exception_boundary_index": 2, "replay_exception_boundary_index": 2,
        "exception_type": "RuntimeError", "replay_exception_type": "RuntimeError",
        "transaction_committed": False, "replay_transaction_committed": False,
        "entity_generations_before": entities, "replay_entity_generations_before": entities,
        "entity_generations_after": entities, "replay_entity_generations_after": entities,
        "active_model": "rollback-model-725", "replay_active_model": "rollback-model-725",
        "database_owner": "headless:rollback-725", "replay_database_owner": "headless:rollback-725",
        "session_owner": "batch:rollback-725", "replay_session_owner": "batch:rollback-725",
        "command_log_sha256": "3" * 64, "replay_command_log_sha256": "3" * 64,
        "database_sha256": "4" * 64, "replay_database_sha256": "4" * 64,
        "batch_result_sha256": "5" * 64, "accepted_batch_result_sha256": "5" * 64,
    }

    generation = "mesh-recipe-725"
    dag = {"geometry": [], "intervals": ["geometry"], "mesh": ["intervals"], "export": ["mesh"]}
    order = ["geometry", "intervals", "mesh", "export"]
    schemes = {"volume:1": "map", "surface:11": "pave"}
    blocks, sidesets = {"block:10": ["volume:1"]}, {"sideset:20": ["surface:11"]}
    row[_RECIPE] = {
        "recipe_generation": generation,
        **_generations(generation, "dag_generation", "parameter_generation", "execution_generation", "scheme_generation", "set_generation", "export_generation", "owner_generation", "result_generation"),
        "dependency_dag": dag, "replay_dependency_dag": dag,
        "parameter_sha256": "6" * 64, "replay_parameter_sha256": "6" * 64,
        "execution_order": order, "replay_execution_order": order,
        "scheme_assignments": schemes, "replay_scheme_assignments": schemes,
        "block_owners": blocks, "replay_block_owners": blocks,
        "sideset_owners": sidesets, "replay_sideset_owners": sidesets,
        "export_generation_id": 725, "replay_export_generation_id": 725,
        "recipe_owner": "headless:recipe-725", "replay_recipe_owner": "headless:recipe-725",
        "recipe_export_sha256": "7" * 64, "accepted_recipe_export_sha256": "7" * 64,
    }
    return row


def test_v42_positive_public_and_source_contracts() -> None:
    row = _with_v42_coreform_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 4


def test_v42_public_periodic_hex_mismatch() -> None:
    row = _with_v42_coreform_identity(summary())
    row[_PERIODIC].update({"nodemap_generation": "periodic-hex-724", "result_target_surface_id": 99, "result_rigid_transform": [[1.0]], "result_node_map": [[101, 999]], "result_target_face_normal": [1.0, 0.0, 0.0], "result_target_interval_counts": [7, 6], "result_minimum_scaled_jacobian": -0.1, "result_block_owner": "block:old", "accepted_periodic_export_sha256": "a" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["periodic_hex_surfaces_use_current_pair_transform_nodemap_orientation_intervals_quality_block_and_export"]


def test_v42_public_boundary_layer_hex_mismatch() -> None:
    row = _with_v42_coreform_identity(summary())
    row[_BOUNDARY].update({"transition_generation": "boundary-layer-hex-724", "result_first_layer_thickness_m": -1.0, "result_growth_ratio": 0.5, "result_layer_count": 0, "result_total_layer_thickness_m": -1.0, "result_transition_topology": ["tet"], "result_wall_normal": [0.0, -1.0, 0.0], "result_minimum_scaled_jacobian": -0.2, "result_block_owner": "block:old", "accepted_boundarylayer_export_sha256": "b" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["boundary_layer_hexes_use_current_thickness_growth_layers_transition_normal_quality_block_and_export"]


def test_v42_source_python_batch_rollback_mismatch() -> None:
    row = _with_v42_coreform_identity(summary())
    row[_ROLLBACK].update({"rollback_generation": "python-rollback-724", "replay_commands": ["reset"], "replay_exception_boundary_index": 0, "replay_exception_type": "ValueError", "replay_transaction_committed": True, "replay_entity_generations_after": {"volume:1": 725}, "replay_active_model": "old", "replay_database_owner": "gui:old", "replay_session_owner": "gui:old", "replay_command_log_sha256": "c" * 64, "accepted_batch_result_sha256": "d" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["python_batches_rollback_current_exception_entities_model_database_session_log_and_result"]


def test_v42_source_mesh_recipe_mismatch() -> None:
    row = _with_v42_coreform_identity(summary())
    row[_RECIPE].update({"dag_generation": "mesh-recipe-724", "replay_dependency_dag": {"export": []}, "replay_parameter_sha256": "e" * 64, "replay_execution_order": ["export", "mesh"], "replay_scheme_assignments": {}, "replay_block_owners": {}, "replay_sideset_owners": {}, "replay_export_generation_id": 724, "replay_recipe_owner": "gui:old", "accepted_recipe_export_sha256": "f" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["mesh_recipes_use_current_dependency_dag_parameters_order_schemes_sets_export_owner_and_result"]


def test_v42_rejects_self_consistent_wrong_periodic_face_orientation() -> None:
    row = _with_v42_coreform_identity(summary())
    row[_PERIODIC]["target_face_normal"] = row[_PERIODIC]["result_target_face_normal"] = [1.0, 0.0, 0.0]
    assert _public_result(row)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_shrinking_boundary_layers() -> None:
    row = _with_v42_coreform_identity(summary())
    identity = row[_BOUNDARY]
    growth = 0.5
    total = identity["first_layer_thickness_m"] * (growth**identity["layer_count"] - 1.0) / (growth - 1.0)
    identity["growth_ratio"] = identity["result_growth_ratio"] = growth
    identity["total_layer_thickness_m"] = identity["result_total_layer_thickness_m"] = total
    assert _public_result(row)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_committed_failed_batch() -> None:
    row = _with_v42_coreform_identity(summary())
    row[_ROLLBACK]["transaction_committed"] = row[_ROLLBACK]["replay_transaction_committed"] = True
    assert _source_result(row)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_non_topological_recipe_order() -> None:
    row = _with_v42_coreform_identity(summary())
    order = ["export", "mesh", "intervals", "geometry"]
    row[_RECIPE]["execution_order"] = row[_RECIPE]["replay_execution_order"] = order
    assert _source_result(row)["status"] == "needs_attention"
