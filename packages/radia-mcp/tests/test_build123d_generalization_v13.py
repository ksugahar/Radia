from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v12 import _public_v12, _source_v12


def _public_v13():
    reference, measured = _public_v12()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["tessellation_tolerance_unit_identity"] = {
                "tessellation_generation": "tessellation-generation-45", "surface_area_generation": "tessellation-generation-45",
                "linear_deflection_value": 0.05, "linear_deflection_unit": "mm", "linear_deflection_scale_to_m": 1.0e-3,
                "area_evaluation_deflection_unit": "mm", "area_evaluation_deflection_scale_to_m": 1.0e-3,
            }
            row["compound_label_topology_identity"] = {
                "boolean_generation": "boolean-generation-45", "label_table_boolean_generation": "boolean-generation-45", "selector_boolean_generation": "boolean-generation-45",
                "label": "mounting_face", "topology_index": 7, "label_topology_index": 7,
                "final_shape_sha256": digest, "selected_subshape_parent_sha256": digest,
            }
    return reference, measured


def _source_v13():
    row = _source_v12()
    identity = row["replay_identity"]
    identity["step_geometry_unit_scale_identity"] = {
        "step_import_generation": "step-import-45", "geometry_coordinate_generation": "step-import-45", "metadata_generation": "step-import-45",
        "geometry_length_unit": "m", "metadata_length_unit": "m", "geometry_scale_to_m": 1.0, "metadata_scale_to_m": 1.0,
    }
    identity["selector_cache_shape_identity"] = {
        "active_shape_generation": "shape-generation-45", "selector_cache_shape_generation": "shape-generation-45", "selected_face_shape_generation": "shape-generation-45",
        "selector_query_sha256": "3" * 64, "cached_selector_query_sha256": "3" * 64, "selected_face_ids": [7, 8], "live_face_ids": [7, 8],
    }
    return row


def test_v13_public_tessellation_tolerance_length_unit_mismatch():
    reference, measured = _public_v13()
    measured["external_cad"][0]["tessellation_tolerance_unit_identity"].update({"area_evaluation_deflection_unit": "m", "area_evaluation_deflection_scale_to_m": 1.0})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["tessellated_area_uses_one_length_unit_tolerance"] is False


def test_v13_public_compound_label_topology_index_previous_boolean():
    reference, measured = _public_v13()
    measured["external_cad"][0]["compound_label_topology_identity"].update({"label_table_boolean_generation": "boolean-generation-44", "selected_subshape_parent_sha256": "4" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["compound_labels_resolve_on_final_boolean_topology"] is False


def test_v13_source_step_geometry_unit_scale_metadata_mismatch():
    row = _source_v13()
    row["replay_identity"]["step_geometry_unit_scale_identity"].update({"metadata_generation": "step-import-44", "metadata_length_unit": "mm", "metadata_scale_to_m": 1.0e-3})
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["step_geometry_and_metadata_share_one_length_unit"] is False


def test_v13_source_selector_cache_previous_shape_generation():
    row = _source_v13()
    row["replay_identity"]["selector_cache_shape_identity"].update({"selector_cache_shape_generation": "shape-generation-44", "cached_selector_query_sha256": "5" * 64})
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["selector_cache_belongs_to_active_shape_generation"] is False
