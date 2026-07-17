from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v18 import _public_v18, _source_v18


def _public_v19():
    reference, measured = _public_v18()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("a" if index == 0 else "b") * 64
            row["boolean_history_subshape_label_fillet_order_identity"] = {
                "boolean_generation": "boolean-53",
                "fillet_generation": "fillet-53",
                "history_boolean_generation": "boolean-53",
                "label_fillet_generation": "fillet-53",
                "edge_order_fillet_generation": "fillet-53",
                "subshape_labels": ["inlet", "outlet"],
                "resolved_subshape_labels": ["inlet", "outlet"],
                "history_face_ids": [51, 52],
                "resolved_face_ids": [51, 52],
                "fillet_edge_ids": [71, 72, 73],
                "resolved_fillet_edge_ids": [71, 72, 73],
                "subshape_history_sha256": digest,
                "resolved_subshape_history_sha256": digest,
            }
            row["assembly_mate_frame_unit_location_generation_identity"] = {
                "assembly_generation": "assembly-53",
                "mate_frame_assembly_generation": "assembly-53",
                "parent_location_assembly_generation": "assembly-53",
                "unit_assembly_generation": "assembly-53",
                "length_unit": "mm",
                "mate_frame_length_unit": "mm",
                "parent_location_length_unit": "mm",
                "mate_names": ["shaft_axis", "housing_axis"],
                "resolved_mate_names": ["shaft_axis", "housing_axis"],
                "local_frame_sha256": ["1" * 64, "2" * 64],
                "resolved_local_frame_sha256": ["1" * 64, "2" * 64],
                "parent_location_sha256": ["3" * 64, "4" * 64],
                "resolved_parent_location_sha256": ["3" * 64, "4" * 64],
                "mate_resolution_sha256": digest,
                "resolved_mate_resolution_sha256": digest,
            }
    return reference, measured


def _source_v19():
    row = _source_v18()
    identity = row["replay_identity"]
    identity["step_import_tolerance_unit_healing_generation_identity"] = {
        "step_import_generation": "step-import-53",
        "healing_generation": "healing-53",
        "tolerance_import_generation": "step-import-53",
        "tolerance_healing_generation": "healing-53",
        "healed_edge_map_import_generation": "step-import-53",
        "healed_edge_map_healing_generation": "healing-53",
        "source_length_unit": "mm",
        "tolerance_length_unit": "mm",
        "tolerance_value": 1.0e-5,
        "tolerance_si_m": 1.0e-8,
        "healed_edge_ids": [101, 102, 103],
        "imported_healed_edge_ids": [101, 102, 103],
        "healed_edge_map_sha256": "c" * 64,
        "imported_healed_edge_map_sha256": "c" * 64,
    }
    identity["tessellation_vertex_index_normal_transform_generation_identity"] = {
        "shape_generation": "shape-53",
        "location_transform_generation": "location-53",
        "vertex_shape_generation": "shape-53",
        "index_shape_generation": "shape-53",
        "normal_shape_generation": "shape-53",
        "vertex_location_transform_generation": "location-53",
        "index_location_transform_generation": "location-53",
        "normal_location_transform_generation": "location-53",
        "vertex_count": 4,
        "triangle_indices": [[0, 1, 2], [0, 2, 3]],
        "normal_count": 4,
        "transformed_triangle_indices": [[0, 1, 2], [0, 2, 3]],
        "tessellation_sha256": "d" * 64,
        "rendered_tessellation_sha256": "d" * 64,
    }
    return row


def test_v19_positive_public_and_source_cad_identity():
    reference, measured = _public_v19()
    public = _public_result(reference, measured)
    assert public["status"] == "ok"
    assert public["checks"][
        "boolean_subshape_labels_follow_current_fillet_edge_order"
    ]
    assert public["checks"][
        "assembly_mates_share_current_frame_unit_and_parent_location"
    ]
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v19()))
    )
    assert source["status"] == "ok"
    assert source["checks"][
        "step_import_tolerances_and_healed_edges_share_current_generation"
    ]
    assert source["checks"][
        "tessellation_vertices_indices_and_normals_use_final_transform"
    ]


def test_v19_public_boolean_history_subshape_label_fillet_reorder_mismatch():
    reference, measured = _public_v19()
    measured["external_cad"][0][
        "boolean_history_subshape_label_fillet_order_identity"
    ].update(
        {
            "label_fillet_generation": "fillet-52",
            "edge_order_fillet_generation": "fillet-52",
            "resolved_subshape_labels": ["outlet", "inlet"],
            "resolved_face_ids": [52, 51],
            "resolved_fillet_edge_ids": [73, 71, 72],
            "resolved_subshape_history_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "boolean_subshape_labels_follow_current_fillet_edge_order"
    ] is False


def test_v19_public_assembly_mate_frame_unit_location_generation_mismatch():
    reference, measured = _public_v19()
    measured["external_cad"][0][
        "assembly_mate_frame_unit_location_generation_identity"
    ].update(
        {
            "mate_frame_assembly_generation": "assembly-52",
            "unit_assembly_generation": "assembly-52",
            "mate_frame_length_unit": "m",
            "resolved_local_frame_sha256": ["2" * 64, "1" * 64],
            "resolved_parent_location_sha256": ["4" * 64, "3" * 64],
            "resolved_mate_resolution_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "assembly_mates_share_current_frame_unit_and_parent_location"
    ] is False


def test_v19_source_step_import_tolerance_unit_healing_generation_mismatch():
    row = _source_v19()
    row["replay_identity"][
        "step_import_tolerance_unit_healing_generation_identity"
    ].update(
        {
            "tolerance_import_generation": "step-import-52",
            "healed_edge_map_healing_generation": "healing-52",
            "tolerance_length_unit": "m",
            "tolerance_si_m": 1.0e-5,
            "imported_healed_edge_ids": [103, 102, 101],
            "imported_healed_edge_map_sha256": "e" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "step_import_tolerances_and_healed_edges_share_current_generation"
    ] is False


def test_v19_source_tessellation_vertex_index_normal_transform_generation_mismatch():
    row = _source_v19()
    row["replay_identity"][
        "tessellation_vertex_index_normal_transform_generation_identity"
    ].update(
        {
            "vertex_location_transform_generation": "location-52",
            "normal_location_transform_generation": "location-52",
            "transformed_triangle_indices": [[0, 2, 1], [0, 3, 2]],
            "rendered_tessellation_sha256": "e" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "tessellation_vertices_indices_and_normals_use_final_transform"
    ] is False
