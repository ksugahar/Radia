from __future__ import annotations

import json

from test_cubit_conformal_mixed_transition_gate import (
    _with_v28_sweep_refinement_sideset_python_identity,
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v29_public_hex_map_curve_interval_parity_vertex_valence_face_pairing_jacobian_block_mismatch",
    "v29_public_mesh_morph_boundary_displacement_constraint_frame_jacobian_sideset_export_mismatch",
    "v29_source_exodus_truth_table_element_variable_block_time_step_order_checksum_mismatch",
    "v29_source_cad_import_body_sheet_lump_transform_unit_heal_transaction_session_mismatch",
)


def _with_v29_map_morph_exodus_cad_identity(row: dict) -> dict:
    row = _with_v28_sweep_refinement_sideset_python_identity(row)
    generation = "hex-map-161"
    row[
        "hex_map_curve_interval_parity_vertex_valence_face_pairing_logical_jacobian_block_generation_identity"
    ] = {
        "map_generation": generation,
        "interval_map_generation": generation,
        "valence_map_generation": generation,
        "pairing_map_generation": generation,
        "logical_map_generation": generation,
        "jacobian_map_generation": generation,
        "block_map_generation": generation,
        "result_map_generation": generation,
        "boundary_curve_ids": [1, 2, 3, 4],
        "result_boundary_curve_ids": [1, 2, 3, 4],
        "curve_intervals": [8, 8, 8, 8],
        "result_curve_intervals": [8, 8, 8, 8],
        "corner_vertex_ids": [11, 12, 13, 14, 15, 16, 17, 18],
        "result_corner_vertex_ids": [11, 12, 13, 14, 15, 16, 17, 18],
        "corner_valences": [3, 3, 3, 3, 3, 3, 3, 3],
        "result_corner_valences": [3, 3, 3, 3, 3, 3, 3, 3],
        "source_target_face_pairs": [[101, 201]],
        "result_source_target_face_pairs": [[101, 201]],
        "logical_corner_coordinates": [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        "result_logical_corner_coordinates": [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        "scaled_jacobians": [0.72, 0.68, 0.61, 0.57],
        "result_scaled_jacobians": [0.72, 0.68, 0.61, 0.57],
        "block_id": 30,
        "result_block_id": 30,
        "map_mesh_sha256": "1" * 64,
        "result_map_mesh_sha256": "1" * 64,
    }

    generation = "mesh-morph-161"
    row[
        "mesh_morph_boundary_displacement_constraint_frame_smoothing_jacobian_sideset_export_generation_identity"
    ] = {
        "morph_generation": generation,
        "boundary_morph_generation": generation,
        "constraint_morph_generation": generation,
        "frame_morph_generation": generation,
        "smoothing_morph_generation": generation,
        "jacobian_morph_generation": generation,
        "sideset_morph_generation": generation,
        "export_morph_generation": generation,
        "result_morph_generation": generation,
        "boundary_node_ids": [101, 102, 103],
        "result_boundary_node_ids": [101, 102, 103],
        "boundary_displacements": [
            [0.1, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.1, 0.0, 0.0],
        ],
        "result_boundary_displacements": [
            [0.1, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.1, 0.0, 0.0],
        ],
        "fixed_node_ids": [1, 2, 3],
        "result_fixed_node_ids": [1, 2, 3],
        "coordinate_frame": "global_cartesian",
        "result_coordinate_frame": "global_cartesian",
        "interior_smoothing": "winslow_volume",
        "result_interior_smoothing": "winslow_volume",
        "minimum_scaled_jacobian": 0.41,
        "result_minimum_scaled_jacobian": 0.41,
        "sideset_ids": [100, 200],
        "result_sideset_ids": [100, 200],
        "morphed_mesh_sha256": "2" * 64,
        "result_morphed_mesh_sha256": "2" * 64,
        "export_sha256": "3" * 64,
        "accepted_export_sha256": "3" * 64,
    }

    generation = "exodus-truth-161"
    row[
        "exodus_truth_table_element_variable_block_timestep_layout_file_generation_identity"
    ] = {
        "truth_generation": generation,
        "variable_truth_generation": generation,
        "block_truth_generation": generation,
        "table_truth_generation": generation,
        "timestep_truth_generation": generation,
        "layout_truth_generation": generation,
        "file_truth_generation": generation,
        "result_truth_generation": generation,
        "element_variable_names": ["stress_xx", "energy"],
        "decoded_element_variable_names": ["stress_xx", "energy"],
        "block_ids": [10, 20],
        "decoded_block_ids": [10, 20],
        "truth_table": [[1, 1], [1, 0]],
        "decoded_truth_table": [[1, 1], [1, 0]],
        "time_step_indices": [1, 2],
        "decoded_time_step_indices": [1, 2],
        "value_layout": "time_block_variable_element",
        "decoded_value_layout": "time_block_variable_element",
        "exodus_sha256": "4" * 64,
        "decoded_exodus_sha256": "4" * 64,
        "truth_table_sha256": "5" * 64,
        "decoded_truth_table_sha256": "5" * 64,
    }

    generation = "cad-import-161"
    row[
        "cad_import_body_sheet_lump_transform_unit_heal_topology_session_result_generation_identity"
    ] = {
        "cad_generation": generation,
        "classification_cad_generation": generation,
        "transform_cad_generation": generation,
        "unit_cad_generation": generation,
        "heal_cad_generation": generation,
        "topology_cad_generation": generation,
        "session_cad_generation": generation,
        "result_cad_generation": generation,
        "body_ids": [1],
        "result_body_ids": [1],
        "sheet_ids": [2],
        "result_sheet_ids": [2],
        "lump_ids": [11, 12],
        "result_lump_ids": [11, 12],
        "placement_matrix": [
            [1, 0, 0, 10],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ],
        "result_placement_matrix": [
            [1, 0, 0, 10],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ],
        "length_unit": "mm",
        "result_length_unit": "mm",
        "heal_transaction_id": "heal-161",
        "result_heal_transaction_id": "heal-161",
        "topology_revision": "topology-161",
        "result_topology_revision": "topology-161",
        "cubit_session_id": "headless-session-161",
        "result_cubit_session_id": "headless-session-161",
        "cad_sha256": "6" * 64,
        "loaded_cad_sha256": "6" * 64,
        "import_result_sha256": "7" * 64,
        "accepted_import_result_sha256": "7" * 64,
    }
    return row


def test_v29_positive_map_morph_exodus_and_cad_identity() -> None:
    row = _with_v29_map_morph_exodus_cad_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v29_public_hex_map_curve_interval_parity_vertex_valence_face_pairing_jacobian_block_mismatch() -> None:
    row = _with_v29_map_morph_exodus_cad_identity(summary())
    row[
        "hex_map_curve_interval_parity_vertex_valence_face_pairing_logical_jacobian_block_generation_identity"
    ].update(
        {
            "interval_map_generation": "hex-map-160",
            "pairing_map_generation": "hex-map-159",
            "result_map_generation": "hex-map-158",
            "result_boundary_curve_ids": [1, 2, 3, 5],
            "result_curve_intervals": [8, 7, 8, 6],
            "result_corner_vertex_ids": [11, 12, 13, 14, 15, 16, 17],
            "result_corner_valences": [3, 4, 3, 3, 3, 3, 2],
            "result_source_target_face_pairs": [[201, 101]],
            "result_logical_corner_coordinates": [[0, 0, 0]],
            "result_scaled_jacobians": [0.72, -0.1],
            "result_block_id": 40,
            "result_map_mesh_sha256": "8" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mapped_hexes_use_current_intervals_valence_face_pairing_logical_coordinates_jacobians_and_blocks"
    ]


def test_v29_public_mesh_morph_boundary_displacement_constraint_frame_jacobian_sideset_export_mismatch() -> None:
    row = _with_v29_map_morph_exodus_cad_identity(summary())
    row[
        "mesh_morph_boundary_displacement_constraint_frame_smoothing_jacobian_sideset_export_generation_identity"
    ].update(
        {
            "boundary_morph_generation": "mesh-morph-160",
            "frame_morph_generation": "mesh-morph-159",
            "result_morph_generation": "mesh-morph-158",
            "result_boundary_node_ids": [101, 103],
            "result_boundary_displacements": [[-0.1, 0.0, 0.0]],
            "result_fixed_node_ids": [2, 3, 4],
            "result_coordinate_frame": "local_rotated",
            "result_interior_smoothing": "none",
            "result_minimum_scaled_jacobian": -0.2,
            "result_sideset_ids": [100, 300],
            "result_morphed_mesh_sha256": "9" * 64,
            "accepted_export_sha256": "a" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mesh_morphs_use_current_boundaries_constraints_frame_smoothing_jacobians_sidesets_and_export"
    ]


def test_v29_source_exodus_truth_table_element_variable_block_time_step_order_checksum_mismatch() -> None:
    row = _with_v29_map_morph_exodus_cad_identity(summary())
    row[
        "exodus_truth_table_element_variable_block_timestep_layout_file_generation_identity"
    ].update(
        {
            "variable_truth_generation": "exodus-truth-160",
            "timestep_truth_generation": "exodus-truth-159",
            "result_truth_generation": "exodus-truth-158",
            "decoded_element_variable_names": ["energy", "stress_xx"],
            "decoded_block_ids": [20, 10],
            "decoded_truth_table": [[1, 0], [0, 1]],
            "decoded_time_step_indices": [2, 1],
            "decoded_value_layout": "variable_time_element_block",
            "decoded_exodus_sha256": "b" * 64,
            "decoded_truth_table_sha256": "c" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "exodus_element_variables_use_current_truth_table_blocks_timesteps_layout_and_file"
    ]


def test_v29_source_cad_import_body_sheet_lump_transform_unit_heal_transaction_session_mismatch() -> None:
    row = _with_v29_map_morph_exodus_cad_identity(summary())
    row[
        "cad_import_body_sheet_lump_transform_unit_heal_topology_session_result_generation_identity"
    ].update(
        {
            "classification_cad_generation": "cad-import-160",
            "heal_cad_generation": "cad-import-159",
            "result_cad_generation": "cad-import-158",
            "result_body_ids": [2],
            "result_sheet_ids": [1],
            "result_lump_ids": [12, 11],
            "result_placement_matrix": [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ],
            "result_length_unit": "m",
            "result_heal_transaction_id": "heal-old",
            "result_topology_revision": "topology-160",
            "result_cubit_session_id": "gui-session-old",
            "loaded_cad_sha256": "d" * 64,
            "accepted_import_result_sha256": "e" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cad_imports_use_current_classification_transform_units_heal_topology_session_and_result"
    ]
