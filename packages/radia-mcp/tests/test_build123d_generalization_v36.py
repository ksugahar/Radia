from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v35 import _public_v35, _source_v35


_PROMOTED_CASE_IDS = (
    "v36_public_fillet_chamfer_edge_selection_radius_tolerance_topology_volume_owner_mismatch",
    "v36_public_loft_section_order_orientation_continuity_volume_centroid_owner_mismatch",
    "v36_source_step_assembly_hierarchy_transform_name_unit_occurrence_digest_mismatch",
    "v36_source_boolean_tolerance_healing_sliver_nonmanifold_history_owner_mismatch",
)


def _public_v36():
    reference, measured = _public_v35()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "fillet-chamfer-contract-231"
            row[
                "fillet_chamfer_edge_selection_radius_distance_tolerance_euler_volume_area_owner_brep_result_generation_identity"
            ] = {
                "feature_generation": generation,
                **{key: generation for key in ("selection_generation", "size_generation", "tolerance_generation", "topology_generation", "volume_generation", "area_generation", "owner_generation", "brep_generation", "result_generation")},
                "feature_kind": "fillet",
                "result_feature_kind": "fillet",
                "selected_edge_ids": [11, 12, 13, 14],
                "result_selected_edge_ids": [11, 12, 13, 14],
                "radius_m": 0.01,
                "result_radius_m": 0.01,
                "distance_m": 0.0,
                "result_distance_m": 0.0,
                "modeling_tolerance_m": 1.0e-6,
                "result_modeling_tolerance_m": 1.0e-6,
                "topology_before_v_e_f": [8, 12, 6],
                "result_topology_before_v_e_f": [8, 12, 6],
                "topology_after_v_e_f": [16, 24, 10],
                "result_topology_after_v_e_f": [16, 24, 10],
                "volume_before_m3": 1.0,
                "result_volume_before_m3": 1.0,
                "volume_after_m3": 0.99,
                "result_volume_after_m3": 0.99,
                "surface_area_before_m2": 6.0,
                "result_surface_area_before_m2": 6.0,
                "surface_area_after_m2": 5.95,
                "result_surface_area_after_m2": 5.95,
                "shape_owner": "part:fillet31",
                "result_shape_owner": "part:fillet31",
                "feature_brep_sha256": suffix * 64,
                "accepted_feature_brep_sha256": suffix * 64,
            }
            generation = "loft-contract-231"
            row[
                "loft_section_order_orientation_guide_continuity_volume_centroid_owner_brep_result_generation_identity"
            ] = {
                "loft_generation": generation,
                **{key: generation for key in ("section_generation", "orientation_generation", "guide_generation", "continuity_generation", "volume_generation", "centroid_generation", "owner_generation", "brep_generation", "result_generation")},
                "section_ids": ["wire:z0", "wire:z1", "wire:z2"],
                "result_section_ids": ["wire:z0", "wire:z1", "wire:z2"],
                "section_parameters": [0.0, 0.5, 1.0],
                "result_section_parameters": [0.0, 0.5, 1.0],
                "section_orientation_signs": [1, 1, 1],
                "result_section_orientation_signs": [1, 1, 1],
                "guide_correspondence": [[0, 0], [1, 1], [2, 2], [3, 3]],
                "result_guide_correspondence": [[0, 0], [1, 1], [2, 2], [3, 3]],
                "continuity": "G1",
                "result_continuity": "G1",
                "loft_volume_m3": 0.75,
                "result_loft_volume_m3": 0.75,
                "loft_centroid_m": [0.0, 0.0, 0.5],
                "result_loft_centroid_m": [0.0, 0.0, 0.5],
                "loft_owner": "part:loft31",
                "result_loft_owner": "part:loft31",
                "loft_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_loft_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v36():
    row = _source_v35()
    identity = row["replay_identity"]
    generation = "step-assembly-contract-231"
    identity[
        "step_assembly_hierarchy_occurrence_transform_repeated_part_name_color_unit_owner_file_result_generation_identity"
    ] = {
        "assembly_generation": generation,
        **{key: generation for key in ("hierarchy_generation", "occurrence_generation", "transform_generation", "part_generation", "name_generation", "color_generation", "unit_generation", "owner_generation", "file_generation", "result_generation")},
        "occurrence_paths": ["root/base", "root/arm:0", "root/arm:1"],
        "decoded_occurrence_paths": ["root/base", "root/arm:0", "root/arm:1"],
        "part_ids": ["base", "arm", "arm"],
        "decoded_part_ids": ["base", "arm", "arm"],
        "transform_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "decoded_transform_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "product_names": ["base", "arm", "arm"],
        "decoded_product_names": ["base", "arm", "arm"],
        "colors_rgb": [[0.2, 0.3, 0.4], [0.8, 0.1, 0.1], [0.8, 0.1, 0.1]],
        "decoded_colors_rgb": [[0.2, 0.3, 0.4], [0.8, 0.1, 0.1], [0.8, 0.1, 0.1]],
        "length_unit": "m",
        "decoded_length_unit": "m",
        "assembly_owner": "assembly:root31",
        "decoded_assembly_owner": "assembly:root31",
        "step_file_sha256": "8" * 64,
        "decoded_step_file_sha256": "8" * 64,
    }
    generation = "boolean-history-contract-231"
    identity[
        "boolean_tolerance_healing_sliver_nonmanifold_operation_history_input_output_owner_brep_result_generation_identity"
    ] = {
        "boolean_generation": generation,
        **{key: generation for key in ("tolerance_generation", "healing_generation", "sliver_generation", "manifold_generation", "history_generation", "input_generation", "output_generation", "owner_generation", "brep_generation", "result_generation")},
        "operation": "fuse",
        "decoded_operation": "fuse",
        "fuzzy_tolerance_m": 1.0e-6,
        "decoded_fuzzy_tolerance_m": 1.0e-6,
        "healing_actions": ["same_domain", "sew"],
        "decoded_healing_actions": ["same_domain", "sew"],
        "sliver_face_count": 0,
        "decoded_sliver_face_count": 0,
        "nonmanifold_edge_count": 0,
        "decoded_nonmanifold_edge_count": 0,
        "input_shape_ids": ["solid:a", "solid:b"],
        "decoded_input_shape_ids": ["solid:a", "solid:b"],
        "output_shape_id": "solid:fused",
        "decoded_output_shape_id": "solid:fused",
        "operation_history": [["solid:a", "solid:fused"], ["solid:b", "solid:fused"]],
        "decoded_operation_history": [["solid:a", "solid:fused"], ["solid:b", "solid:fused"]],
        "input_owner": "boolean:inputs31",
        "decoded_input_owner": "boolean:inputs31",
        "output_owner": "boolean:result31",
        "decoded_output_owner": "boolean:result31",
        "boolean_brep_sha256": "9" * 64,
        "decoded_boolean_brep_sha256": "9" * 64,
    }
    return row


def test_v36_positive_all_four_contracts():
    reference, measured = _public_v36()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v36())["status"] == "ok"


def test_v36_public_fillet_chamfer_edge_selection_radius_tolerance_topology_volume_owner_mismatch():
    reference, measured = _public_v36()
    identity = measured["external_cad"][0]["fillet_chamfer_edge_selection_radius_distance_tolerance_euler_volume_area_owner_brep_result_generation_identity"]
    identity.update({"selection_generation": "fillet-chamfer-contract-230", "result_feature_kind": "chamfer", "result_selected_edge_ids": [99], "result_radius_m": -0.01, "result_distance_m": 0.02, "result_modeling_tolerance_m": 0.1, "result_topology_after_v_e_f": [8, 24, 10], "result_volume_after_m3": 1.2, "result_surface_area_after_m2": -1.0, "result_shape_owner": "stale:part", "accepted_feature_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fillet_chamfer_features_use_current_edges_size_tolerance_topology_volume_area_owner_and_brep"]


def test_v36_public_loft_section_order_orientation_continuity_volume_centroid_owner_mismatch():
    reference, measured = _public_v36()
    identity = measured["external_cad"][0]["loft_section_order_orientation_guide_continuity_volume_centroid_owner_brep_result_generation_identity"]
    identity.update({"section_generation": "loft-contract-230", "result_section_ids": ["wire:z2", "wire:z0"], "result_section_parameters": [1.0, 0.0], "result_section_orientation_signs": [1, -1], "result_guide_correspondence": [[0, 3], [1, 2]], "result_continuity": "C0", "result_loft_volume_m3": -0.75, "result_loft_centroid_m": [0.0, 0.0, 2.0], "result_loft_owner": "stale:loft", "accepted_loft_brep_sha256": "b" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["lofts_use_current_sections_order_orientation_guides_continuity_volume_centroid_owner_and_brep"]


def test_v36_source_step_assembly_hierarchy_transform_name_unit_occurrence_digest_mismatch():
    row = _source_v36(); identity = row["replay_identity"]["step_assembly_hierarchy_occurrence_transform_repeated_part_name_color_unit_owner_file_result_generation_identity"]
    identity.update({"hierarchy_generation": "step-assembly-contract-230", "decoded_occurrence_paths": ["root/old"], "decoded_part_ids": ["old"], "decoded_transform_sha256": ["c" * 64], "decoded_product_names": ["old"], "decoded_colors_rgb": [[0.0, 0.0, 0.0]], "decoded_length_unit": "mm", "decoded_assembly_owner": "stale:assembly", "decoded_step_file_sha256": "d" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["step_assemblies_use_current_hierarchy_occurrences_repeated_parts_transforms_names_colors_units_owner_and_file"]


def test_v36_source_boolean_tolerance_healing_sliver_nonmanifold_history_owner_mismatch():
    row = _source_v36(); identity = row["replay_identity"]["boolean_tolerance_healing_sliver_nonmanifold_operation_history_input_output_owner_brep_result_generation_identity"]
    identity.update({"healing_generation": "boolean-history-contract-230", "decoded_operation": "cut", "decoded_fuzzy_tolerance_m": -1.0, "decoded_healing_actions": ["discard"], "decoded_sliver_face_count": 3, "decoded_nonmanifold_edge_count": 2, "decoded_input_shape_ids": ["solid:a"], "decoded_output_shape_id": "solid:old", "decoded_operation_history": [["solid:x", "solid:old"]], "decoded_input_owner": "stale:inputs", "decoded_output_owner": "stale:result", "decoded_boolean_brep_sha256": "e" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["boolean_roundtrips_use_current_operation_tolerance_healing_slivers_manifold_history_owners_and_brep"]


def test_v36_rejects_self_consistent_fillet_euler_error():
    reference, measured = _public_v36(); identity = measured["external_cad"][0]["fillet_chamfer_edge_selection_radius_distance_tolerance_euler_volume_area_owner_brep_result_generation_identity"]
    identity["topology_after_v_e_f"] = [16, 24, 9]; identity["result_topology_after_v_e_f"] = [16, 24, 9]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v36_rejects_self_consistent_boolean_nonmanifold_output():
    row = _source_v36(); identity = row["replay_identity"]["boolean_tolerance_healing_sliver_nonmanifold_operation_history_input_output_owner_brep_result_generation_identity"]
    identity["nonmanifold_edge_count"] = 1; identity["decoded_nonmanifold_edge_count"] = 1
    assert _source_result(row)["status"] == "needs_attention"
