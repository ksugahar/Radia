from __future__ import annotations

from test_coreform_generalization_v37 import _public_result, _source_result, summary
from test_coreform_generalization_v39 import _generations
from test_coreform_generalization_v40 import _with_v40_coreform_identity


_MULTISWEEP = "multisweep_source_target_chain_section_correspondence_twist_interval_quality_block_export_generation_identity"
_IMPRINT = "imprint_merge_tolerance_coincident_topology_volume_block_sideset_quality_export_generation_identity"
_JOURNAL = "journal_undo_redo_transaction_entityid_selection_database_owner_generation_identity"
_EXODUS = "exodus_transient_timestep_nodal_element_variable_truth_sideset_mesh_database_result_generation_identity"
_PROMOTED_CASE_IDS = (
    "v41_public_multisweep_source_target_chain_section_correspondence_twist_interval_hexquality_export_mismatch",
    "v41_public_imprint_merge_tolerance_coincident_topology_block_sideset_hex_export_mismatch",
    "v41_source_journal_undo_redo_transaction_entityid_selection_database_owner_mismatch",
    "v41_source_exodus_transient_timestep_nodal_element_variable_truth_table_owner_mismatch",
)


def _with_v41_coreform_identity(row: dict) -> dict:
    row = _with_v40_coreform_identity(row)
    generation = "multisweep-hex-724"
    row[_MULTISWEEP] = {
        "multisweep_generation": generation,
        **_generations(generation, "source_generation", "target_generation", "section_generation", "twist_generation", "interval_generation", "quality_generation", "block_generation", "export_generation", "result_generation"),
        "source_surface_ids": [11, 12], "result_source_surface_ids": [11, 12],
        "target_surface_ids": [31, 32], "result_target_surface_ids": [31, 32],
        "sweep_volume_chain": [1, 2], "result_sweep_volume_chain": [1, 2],
        "section_node_correspondence": [[101, 201], [102, 202], [103, 203], [104, 204]],
        "result_section_node_correspondence": [[101, 201], [102, 202], [103, 203], [104, 204]],
        "section_twist_deg": [0.0, 15.0], "result_section_twist_deg": [0.0, 15.0],
        "chain_interval_counts": [8, 8], "result_chain_interval_counts": [8, 8],
        "hex_element_count": 512, "result_hex_element_count": 512,
        "minimum_scaled_jacobian": 0.35, "result_minimum_scaled_jacobian": 0.35,
        "minimum_allowed_scaled_jacobian": 0.20, "result_minimum_allowed_scaled_jacobian": 0.20,
        "block_owner": "block:multisweep-41", "result_block_owner": "block:multisweep-41",
        "multisweep_export_sha256": "1" * 64, "accepted_multisweep_export_sha256": "1" * 64,
    }
    generation = "imprint-merge-724"
    row[_IMPRINT] = {
        "imprint_merge_generation": generation,
        **_generations(generation, "tolerance_generation", "coincident_generation", "topology_generation", "volume_generation", "set_generation", "quality_generation", "export_generation", "result_generation"),
        "merge_tolerance_m": 1.0e-7, "result_merge_tolerance_m": 1.0e-7,
        "coincident_vertex_pairs": [[1, 101], [2, 102]], "result_coincident_vertex_pairs": [[1, 101], [2, 102]],
        "coincident_surface_pairs": [[11, 21]], "result_coincident_surface_pairs": [[11, 21]],
        "merged_entity_map": {"vertex:101": "vertex:1", "surface:21": "surface:11"},
        "result_merged_entity_map": {"vertex:101": "vertex:1", "surface:21": "surface:11"},
        "volume_owners": {"volume:1": "rotor", "volume:2": "stator"},
        "result_volume_owners": {"volume:1": "rotor", "volume:2": "stator"},
        "block_membership": {"block:10": ["volume:1"], "block:20": ["volume:2"]},
        "result_block_membership": {"block:10": ["volume:1"], "block:20": ["volume:2"]},
        "sideset_membership": {"sideset:30": ["surface:11"]}, "result_sideset_membership": {"sideset:30": ["surface:11"]},
        "hex_element_count": 640, "result_hex_element_count": 640,
        "minimum_scaled_jacobian": 0.31, "result_minimum_scaled_jacobian": 0.31,
        "imprint_export_sha256": "2" * 64, "accepted_imprint_export_sha256": "2" * 64,
    }
    generation = "journal-replay-724"
    row[_JOURNAL] = {
        "journal_generation": generation,
        **_generations(generation, "transaction_generation", "undo_generation", "redo_generation", "entity_generation", "selection_generation", "database_generation", "owner_generation", "result_generation"),
        "transactions": ["create volume 1", "mesh volume 1", "undo", "redo"],
        "replay_transactions": ["create volume 1", "mesh volume 1", "undo", "redo"],
        "undo_boundary_index": 2, "replay_undo_boundary_index": 2,
        "redo_boundary_index": 3, "replay_redo_boundary_index": 3,
        "entity_id_remap": {"volume:1": "volume:1", "surface:11": "surface:11"},
        "replay_entity_id_remap": {"volume:1": "volume:1", "surface:11": "surface:11"},
        "active_selection": ["volume:1"], "replay_active_selection": ["volume:1"],
        "geometry_generation_id": 724, "replay_geometry_generation_id": 724,
        "database_owner": "headless:journal-724", "replay_database_owner": "headless:journal-724",
        "session_owner": "batch:journal-724", "replay_session_owner": "batch:journal-724",
        "database_sha256": "3" * 64, "replay_database_sha256": "3" * 64,
        "journal_result_sha256": "4" * 64, "accepted_journal_result_sha256": "4" * 64,
    }
    generation = "exodus-transient-724"
    row[_EXODUS] = {
        "exodus_generation": generation,
        **_generations(generation, "time_generation", "nodal_generation", "element_generation", "truth_generation", "sideset_generation", "mesh_generation", "database_generation", "result_generation"),
        "time_steps_s": [0.0, 0.1, 0.2], "replay_time_steps_s": [0.0, 0.1, 0.2],
        "nodal_variable_names": ["temperature"], "replay_nodal_variable_names": ["temperature"],
        "nodal_values": [[300.0, 301.0], [302.0, 303.0], [304.0, 305.0]],
        "replay_nodal_values": [[300.0, 301.0], [302.0, 303.0], [304.0, 305.0]],
        "element_variable_names": ["energy"], "replay_element_variable_names": ["energy"],
        "element_values": [[1.0], [1.1], [1.2]], "replay_element_values": [[1.0], [1.1], [1.2]],
        "block_truth_table": {"block:10": [True]}, "replay_block_truth_table": {"block:10": [True]},
        "sideset_values": {"sideset:20": [0.0, 0.1, 0.2]}, "replay_sideset_values": {"sideset:20": [0.0, 0.1, 0.2]},
        "mesh_generation_id": 724, "replay_mesh_generation_id": 724,
        "database_owner": "headless:exodus-724", "replay_database_owner": "headless:exodus-724",
        "database_sha256": "5" * 64, "replay_database_sha256": "5" * 64,
        "exodus_result_sha256": "6" * 64, "accepted_exodus_result_sha256": "6" * 64,
    }
    return row


def test_v41_positive_public_and_source_contracts() -> None:
    row = _with_v41_coreform_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 4


def test_v41_public_multisweep_mismatch() -> None:
    row = _with_v41_coreform_identity(summary())
    row[_MULTISWEEP].update({"section_generation": "multisweep-hex-723", "result_target_surface_ids": [99], "result_sweep_volume_chain": [2, 1], "result_section_node_correspondence": [[101, 204]], "result_chain_interval_counts": [7, 8], "result_minimum_scaled_jacobian": -0.1, "accepted_multisweep_export_sha256": "a" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["multisweep_hexes_use_current_sources_targets_chain_sections_twist_intervals_quality_block_and_export"]


def test_v41_public_imprint_merge_mismatch() -> None:
    row = _with_v41_coreform_identity(summary())
    row[_IMPRINT].update({"topology_generation": "imprint-merge-723", "result_merge_tolerance_m": 1.0e-2, "result_coincident_surface_pairs": [], "result_merged_entity_map": {}, "result_block_membership": {}, "result_sideset_membership": {}, "result_hex_element_count": 0, "accepted_imprint_export_sha256": "b" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["imprint_merge_hexes_use_current_tolerance_coincidence_topology_sets_quality_and_export"]


def test_v41_source_journal_mismatch() -> None:
    row = _with_v41_coreform_identity(summary())
    row[_JOURNAL].update({"undo_generation": "journal-replay-723", "replay_transactions": ["create volume 1", "redo"], "replay_entity_id_remap": {"volume:1": "volume:9"}, "replay_active_selection": ["volume:9"], "replay_database_owner": "gui:old", "accepted_journal_result_sha256": "d" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["journal_undo_redo_replays_use_current_transactions_entities_selection_database_owner_and_result"]


def test_v41_source_exodus_mismatch() -> None:
    row = _with_v41_coreform_identity(summary())
    row[_EXODUS].update({"time_generation": "exodus-transient-723", "replay_time_steps_s": [0.2, 0.0], "replay_nodal_variable_names": ["old"], "replay_nodal_values": [[999.0]], "replay_block_truth_table": {"block:99": [False]}, "replay_database_owner": "gui:old", "accepted_exodus_result_sha256": "f" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["exodus_transients_use_current_times_variables_truth_sets_mesh_database_and_result"]


def test_v41_rejects_self_consistent_incompatible_multisweep_intervals() -> None:
    row = _with_v41_coreform_identity(summary())
    row[_MULTISWEEP]["chain_interval_counts"] = row[_MULTISWEEP]["result_chain_interval_counts"] = [8, 7]
    assert _public_result(row)["status"] == "needs_attention"


def test_v41_rejects_self_consistent_oversized_merge_tolerance() -> None:
    row = _with_v41_coreform_identity(summary())
    row[_IMPRINT]["merge_tolerance_m"] = row[_IMPRINT]["result_merge_tolerance_m"] = 1.0e-2
    assert _public_result(row)["status"] == "needs_attention"


def test_v41_rejects_self_consistent_redo_before_undo() -> None:
    row = _with_v41_coreform_identity(summary())
    transactions = ["create volume 1", "redo", "mesh volume 1", "undo"]
    row[_JOURNAL]["transactions"] = row[_JOURNAL]["replay_transactions"] = transactions
    row[_JOURNAL]["undo_boundary_index"] = row[_JOURNAL]["replay_undo_boundary_index"] = 3
    row[_JOURNAL]["redo_boundary_index"] = row[_JOURNAL]["replay_redo_boundary_index"] = 1
    assert _source_result(row)["status"] == "needs_attention"


def test_v41_rejects_self_consistent_nonmonotone_exodus_time() -> None:
    row = _with_v41_coreform_identity(summary())
    row[_EXODUS]["time_steps_s"] = row[_EXODUS]["replay_time_steps_s"] = [0.0, 0.2, 0.1]
    assert _source_result(row)["status"] == "needs_attention"
