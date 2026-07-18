from __future__ import annotations

from test_build123d_generalization_v29 import (
    _public_result,
    _public_v29,
    _source_result,
    _source_v29,
)


_PROMOTED_CASE_IDS = (
    "v30_public_helical_sweep_pitch_handedness_profile_frame_self_intersection_volume_mismatch",
    "v30_public_boolean_tolerance_operand_order_volume_centroid_inertia_history_mismatch",
    "v30_source_step_assembly_product_name_color_transform_unit_instance_digest_mismatch",
    "v30_source_stl_facet_normal_winding_watertight_tolerance_volume_digest_mismatch",
)


def _public_v30():
    reference, measured = _public_v29()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            shape_digest = ("1" if index == 0 else "2") * 64
            profile_digest = ("3" if index == 0 else "4") * 64
            generation = "helical-sweep-171"
            volume = 1.2e-6 + index * 0.3e-6
            centroid = [0.0, 0.0, 0.02 + index * 0.005]
            frame = [[1.0, 0.0, 0.0, 0.005], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
            row[
                "helical_sweep_pitch_handedness_profile_frame_turn_self_intersection_volume_centroid_shape_generation_identity"
            ] = {
                "helical_generation": generation,
                "pitch_helical_generation": generation,
                "handedness_helical_generation": generation,
                "profile_helical_generation": generation,
                "turn_helical_generation": generation,
                "intersection_helical_generation": generation,
                "mass_helical_generation": generation,
                "shape_helical_generation": generation,
                "result_helical_generation": generation,
                "pitch_m": 0.01,
                "result_pitch_m": 0.01,
                "handedness": "right",
                "result_handedness": "right",
                "profile_frame_matrix": frame,
                "result_profile_frame_matrix": frame,
                "profile_frame_sha256": profile_digest,
                "result_profile_frame_sha256": profile_digest,
                "turn_count": 4.0,
                "result_turn_count": 4.0,
                "self_intersection": False,
                "result_self_intersection": False,
                "volume_m3": volume,
                "result_volume_m3": volume,
                "centroid_m": centroid,
                "result_centroid_m": centroid,
                "helical_shape_sha256": shape_digest,
                "result_helical_shape_sha256": shape_digest,
            }
            generation = "boolean-history-171"
            boolean_digest = ("5" if index == 0 else "6") * 64
            history_digest = ("7" if index == 0 else "8") * 64
            inertia = [[1.0e-9, 0.0, 0.0], [0.0, 2.0e-9, 0.0], [0.0, 0.0, 3.0e-9]]
            row[
                "boolean_tolerance_operand_order_history_volume_centroid_inertia_shape_generation_identity"
            ] = {
                "boolean_generation": generation,
                "tolerance_boolean_generation": generation,
                "operand_boolean_generation": generation,
                "history_boolean_generation": generation,
                "mass_boolean_generation": generation,
                "inertia_boolean_generation": generation,
                "shape_boolean_generation": generation,
                "result_boolean_generation": generation,
                "operation": "cut",
                "result_operation": "cut",
                "model_tolerance_m": 1.0e-7,
                "result_model_tolerance_m": 1.0e-7,
                "operand_order": ["body-a", "tool-b"],
                "result_operand_order": ["body-a", "tool-b"],
                "subshape_history_sha256": history_digest,
                "result_subshape_history_sha256": history_digest,
                "volume_m3": volume,
                "result_volume_m3": volume,
                "centroid_m": centroid,
                "result_centroid_m": centroid,
                "inertia_tensor_kg_m2": inertia,
                "result_inertia_tensor_kg_m2": inertia,
                "boolean_shape_sha256": boolean_digest,
                "result_boolean_shape_sha256": boolean_digest,
            }
    return reference, measured


def _source_v30():
    row = _source_v29()
    identity = row["replay_identity"]
    generation = "step-assembly-171"
    identity[
        "step_assembly_product_color_instance_transform_unit_hierarchy_file_generation_identity"
    ] = {
        "step_generation": generation,
        "product_step_generation": generation,
        "color_step_generation": generation,
        "instance_step_generation": generation,
        "transform_step_generation": generation,
        "unit_step_generation": generation,
        "hierarchy_step_generation": generation,
        "file_step_generation": generation,
        "result_step_generation": generation,
        "product_names": ["assembly", "housing", "shaft"],
        "decoded_product_names": ["assembly", "housing", "shaft"],
        "product_colors_rgb": [["housing", 0.2, 0.4, 0.8], ["shaft", 0.7, 0.7, 0.7]],
        "decoded_product_colors_rgb": [["housing", 0.2, 0.4, 0.8], ["shaft", 0.7, 0.7, 0.7]],
        "instance_order": ["housing-1", "shaft-1"],
        "decoded_instance_order": ["housing-1", "shaft-1"],
        "instance_transform_sha256": [["housing-1", "9" * 64], ["shaft-1", "a" * 64]],
        "decoded_instance_transform_sha256": [["housing-1", "9" * 64], ["shaft-1", "a" * 64]],
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "assembly_hierarchy": [["assembly", "housing-1"], ["assembly", "shaft-1"]],
        "decoded_assembly_hierarchy": [["assembly", "housing-1"], ["assembly", "shaft-1"]],
        "step_sha256": "b" * 64,
        "decoded_step_sha256": "b" * 64,
        "assembly_shape_sha256": "c" * 64,
        "decoded_assembly_shape_sha256": "c" * 64,
    }
    generation = "stl-solid-171"
    identity[
        "stl_facet_normal_winding_watertight_tolerance_volume_unit_file_generation_identity"
    ] = {
        "stl_generation": generation,
        "normal_stl_generation": generation,
        "winding_stl_generation": generation,
        "watertight_stl_generation": generation,
        "tolerance_stl_generation": generation,
        "volume_stl_generation": generation,
        "unit_stl_generation": generation,
        "file_stl_generation": generation,
        "result_stl_generation": generation,
        "facet_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
        "decoded_facet_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
        "triangle_winding": "outward_counterclockwise",
        "decoded_triangle_winding": "outward_counterclockwise",
        "unmatched_edge_count": 0,
        "decoded_unmatched_edge_count": 0,
        "merge_tolerance_m": 1.0e-7,
        "decoded_merge_tolerance_m": 1.0e-7,
        "signed_volume_m3": 1.5e-6,
        "decoded_signed_volume_m3": 1.5e-6,
        "length_unit": "m",
        "decoded_length_unit": "m",
        "stl_sha256": "d" * 64,
        "decoded_stl_sha256": "d" * 64,
        "stl_solid_sha256": "e" * 64,
        "decoded_stl_solid_sha256": "e" * 64,
    }
    return row


def test_v30_positive_public_and_source_identity():
    reference, measured = _public_v30()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v30())["status"] == "ok"


def test_v30_public_helical_sweep_pitch_handedness_profile_frame_self_intersection_volume_mismatch():
    reference, measured = _public_v30()
    measured["external_cad"][0][
        "helical_sweep_pitch_handedness_profile_frame_turn_self_intersection_volume_centroid_shape_generation_identity"
    ].update(
        {
            "pitch_helical_generation": "helical-sweep-170",
            "profile_helical_generation": "helical-sweep-169",
            "result_helical_generation": "helical-sweep-168",
            "result_pitch_m": 0.02,
            "result_handedness": "left",
            "result_profile_frame_matrix": [[1.0, 0.0, 0.0, 0.0]],
            "result_profile_frame_sha256": "f" * 64,
            "result_turn_count": 3.5,
            "result_self_intersection": True,
            "result_volume_m3": 0.8e-6,
            "result_centroid_m": [0.01, 0.0, 0.01],
            "result_helical_shape_sha256": "0" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "helical_sweeps_use_current_pitch_handedness_profile_frame_turns_intersection_mass_and_shape"
    ]


def test_v30_public_boolean_tolerance_operand_order_volume_centroid_inertia_history_mismatch():
    reference, measured = _public_v30()
    measured["external_cad"][0][
        "boolean_tolerance_operand_order_history_volume_centroid_inertia_shape_generation_identity"
    ].update(
        {
            "tolerance_boolean_generation": "boolean-history-170",
            "history_boolean_generation": "boolean-history-169",
            "result_boolean_generation": "boolean-history-168",
            "result_operation": "fuse",
            "result_model_tolerance_m": 1.0e-3,
            "result_operand_order": ["tool-b", "body-a"],
            "result_subshape_history_sha256": "1" * 64,
            "result_volume_m3": 2.4e-6,
            "result_centroid_m": [0.02, 0.0, 0.0],
            "result_inertia_tensor_kg_m2": [[3.0e-9, 1.0e-9, 0.0]],
            "result_boolean_shape_sha256": "2" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_results_use_current_operation_tolerance_operands_history_mass_inertia_and_shape"
    ]


def test_v30_source_step_assembly_product_name_color_transform_unit_instance_digest_mismatch():
    row = _source_v30()
    row["replay_identity"][
        "step_assembly_product_color_instance_transform_unit_hierarchy_file_generation_identity"
    ].update(
        {
            "product_step_generation": "step-assembly-170",
            "transform_step_generation": "step-assembly-169",
            "result_step_generation": "step-assembly-168",
            "decoded_product_names": ["assembly", "shaft", "housing-old"],
            "decoded_product_colors_rgb": [["housing", 0.8, 0.4, 0.2]],
            "decoded_instance_order": ["shaft-1", "housing-1"],
            "decoded_instance_transform_sha256": [["housing-1", "3" * 64]],
            "decoded_length_unit": "m",
            "decoded_assembly_hierarchy": [["housing-1", "assembly"]],
            "decoded_step_sha256": "4" * 64,
            "decoded_assembly_shape_sha256": "5" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_assemblies_use_current_products_colors_instances_transforms_units_hierarchy_and_digests"
    ]


def test_v30_source_stl_facet_normal_winding_watertight_tolerance_volume_digest_mismatch():
    row = _source_v30()
    row["replay_identity"][
        "stl_facet_normal_winding_watertight_tolerance_volume_unit_file_generation_identity"
    ].update(
        {
            "normal_stl_generation": "stl-solid-170",
            "watertight_stl_generation": "stl-solid-169",
            "result_stl_generation": "stl-solid-168",
            "decoded_facet_normals": [[0.0, 0.0, -1.0]],
            "decoded_triangle_winding": "inward_clockwise",
            "decoded_unmatched_edge_count": 4,
            "decoded_merge_tolerance_m": 1.0e-3,
            "decoded_signed_volume_m3": -1.2e-6,
            "decoded_length_unit": "mm",
            "decoded_stl_sha256": "6" * 64,
            "decoded_stl_solid_sha256": "7" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "stl_solids_use_current_normals_winding_watertight_edges_tolerance_volume_units_and_digests"
    ]
