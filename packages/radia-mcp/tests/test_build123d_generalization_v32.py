from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v31 import _public_v31, _source_v31


_PROMOTED_CASE_IDS = (
    "v32_public_boolean_fuzzy_tolerance_topology_name_face_ancestry_volume_centroid_mismatch",
    "v32_public_sweep_frenet_frame_twist_transition_self_intersection_volume_mismatch",
    "v32_source_step_assembly_instance_transform_unit_color_material_uuid_digest_mismatch",
    "v32_source_stl_repair_tolerance_normal_duplicate_vertex_watertight_volume_digest_mismatch",
)


def _public_v32():
    reference, measured = _public_v31()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            suffix = "1" if index == 0 else "2"
            generation = "boolean-topology-closure-191"
            row[
                "boolean_fuzzy_tolerance_topology_name_face_ancestry_count_volume_centroid_shape_generation_identity"
            ] = {
                "boolean_generation": generation,
                "tolerance_generation": generation,
                "topology_generation": generation,
                "ancestry_generation": generation,
                "mass_generation": generation,
                "shape_generation": generation,
                "result_generation": generation,
                "fuzzy_tolerance_m": 1.0e-7,
                "result_fuzzy_tolerance_m": 1.0e-7,
                "surviving_topology_names": ["body", "cut_face", "outer_shell"],
                "result_surviving_topology_names": ["body", "cut_face", "outer_shell"],
                "face_ancestry": [["box:f1", "result:f3"], ["cylinder:f2", "result:f7"]],
                "result_face_ancestry": [["box:f1", "result:f3"], ["cylinder:f2", "result:f7"]],
                "solid_count": 1,
                "result_solid_count": 1,
                "volume_m3": 9.0e-4,
                "result_volume_m3": 9.0e-4,
                "centroid_m": [0.01, 0.0, 0.0],
                "result_centroid_m": [0.01, 0.0, 0.0],
                "boolean_shape_sha256": suffix * 64,
                "result_boolean_shape_sha256": suffix * 64,
            }
            generation = "sweep-frame-closure-191"
            sweep_suffix = "3" if index == 0 else "4"
            row[
                "sweep_frenet_frame_twist_transition_profile_orientation_self_intersection_volume_owner_shape_generation_identity"
            ] = {
                "sweep_generation": generation,
                "frame_generation": generation,
                "twist_generation": generation,
                "transition_generation": generation,
                "orientation_generation": generation,
                "intersection_generation": generation,
                "mass_generation": generation,
                "owner_generation": generation,
                "result_generation": generation,
                "frame_convention": "corrected_frenet",
                "result_frame_convention": "corrected_frenet",
                "twist_parameters_rad": [0.0, 0.2, 0.4],
                "result_twist_parameters_rad": [0.0, 0.2, 0.4],
                "transition_mode": "right_corner",
                "result_transition_mode": "right_corner",
                "profile_orientation_signs": [1, 1, 1],
                "result_profile_orientation_signs": [1, 1, 1],
                "self_intersection": False,
                "result_self_intersection": False,
                "sweep_volume_m3": 1.2e-3,
                "result_sweep_volume_m3": 1.2e-3,
                "shape_owner": "sweep/body1",
                "result_shape_owner": "sweep/body1",
                "sweep_shape_sha256": sweep_suffix * 64,
                "result_sweep_shape_sha256": sweep_suffix * 64,
            }
    return reference, measured


def _source_v32():
    row = _source_v31()
    identity = row["replay_identity"]
    generation = "step-assembly-closure-191"
    identity[
        "step_assembly_instance_transform_unit_color_material_uuid_component_volume_file_generation_identity"
    ] = {
        "step_generation": generation,
        "instance_generation": generation,
        "transform_generation": generation,
        "unit_generation": generation,
        "color_generation": generation,
        "material_generation": generation,
        "uuid_generation": generation,
        "component_generation": generation,
        "volume_generation": generation,
        "file_generation": generation,
        "result_generation": generation,
        "instance_ids": ["housing-1", "shaft-1"],
        "decoded_instance_ids": ["housing-1", "shaft-1"],
        "instance_transform_sha256": [["housing-1", "5" * 64], ["shaft-1", "6" * 64]],
        "decoded_instance_transform_sha256": [["housing-1", "5" * 64], ["shaft-1", "6" * 64]],
        "length_unit": "m",
        "decoded_length_unit": "m",
        "instance_colors_rgb": [["housing-1", 0.8, 0.8, 0.8], ["shaft-1", 0.3, 0.3, 0.3]],
        "decoded_instance_colors_rgb": [["housing-1", 0.8, 0.8, 0.8], ["shaft-1", 0.3, 0.3, 0.3]],
        "material_labels": [["housing-1", "aluminum"], ["shaft-1", "steel"]],
        "decoded_material_labels": [["housing-1", "aluminum"], ["shaft-1", "steel"]],
        "product_uuids": [["housing-1", "uuid-housing-v4"], ["shaft-1", "uuid-shaft-v6"]],
        "decoded_product_uuids": [["housing-1", "uuid-housing-v4"], ["shaft-1", "uuid-shaft-v6"]],
        "component_shape_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64]],
        "decoded_component_shape_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64]],
        "total_volume_m3": 1.2e-3,
        "decoded_total_volume_m3": 1.2e-3,
        "step_file_sha256": "9" * 64,
        "decoded_step_file_sha256": "9" * 64,
    }
    generation = "stl-repair-closure-191"
    identity[
        "stl_repair_merge_tolerance_normal_duplicate_boundary_watertight_volume_unit_file_generation_identity"
    ] = {
        "repair_generation": generation,
        "tolerance_generation": generation,
        "normal_generation": generation,
        "duplicate_generation": generation,
        "boundary_generation": generation,
        "watertight_generation": generation,
        "volume_generation": generation,
        "unit_generation": generation,
        "file_generation": generation,
        "result_generation": generation,
        "merge_tolerance_m": 1.0e-6,
        "decoded_merge_tolerance_m": 1.0e-6,
        "facet_normal_sha256": "a" * 64,
        "decoded_facet_normal_sha256": "a" * 64,
        "duplicate_vertex_count": 0,
        "decoded_duplicate_vertex_count": 0,
        "boundary_edge_count": 0,
        "decoded_boundary_edge_count": 0,
        "watertight_component_count": 1,
        "decoded_watertight_component_count": 1,
        "repaired_volume_m3": 1.5e-3,
        "decoded_repaired_volume_m3": 1.5e-3,
        "length_unit": "m",
        "decoded_length_unit": "m",
        "source_stl_sha256": "b" * 64,
        "decoded_source_stl_sha256": "b" * 64,
        "repaired_stl_sha256": "c" * 64,
        "decoded_repaired_stl_sha256": "c" * 64,
    }
    return row


def test_v32_positive_public_and_source_identity():
    reference, measured = _public_v32()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v32())["status"] == "ok"


def test_v32_public_boolean_fuzzy_tolerance_topology_name_face_ancestry_volume_centroid_mismatch():
    reference, measured = _public_v32()
    measured["external_cad"][0][
        "boolean_fuzzy_tolerance_topology_name_face_ancestry_count_volume_centroid_shape_generation_identity"
    ].update(
        {
            "tolerance_generation": "boolean-topology-190",
            "ancestry_generation": "boolean-topology-189",
            "result_generation": "boolean-topology-188",
            "result_fuzzy_tolerance_m": 1.0e-3,
            "result_surviving_topology_names": ["body", "face_old"],
            "result_face_ancestry": [["box:f2", "result:f3"]],
            "result_solid_count": 2,
            "result_volume_m3": 8.0e-4,
            "result_centroid_m": [0.02, 0.01, 0.0],
            "result_boolean_shape_sha256": "d" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_results_use_current_fuzzy_tolerance_topology_names_face_ancestry_count_volume_centroid_and_shape"
    ]


def test_v32_public_sweep_frenet_frame_twist_transition_self_intersection_volume_mismatch():
    reference, measured = _public_v32()
    measured["external_cad"][0][
        "sweep_frenet_frame_twist_transition_profile_orientation_self_intersection_volume_owner_shape_generation_identity"
    ].update(
        {
            "frame_generation": "sweep-frame-190",
            "intersection_generation": "sweep-frame-189",
            "result_generation": "sweep-frame-188",
            "result_frame_convention": "fixed_global",
            "result_twist_parameters_rad": [0.4, 0.2, 0.0],
            "result_transition_mode": "transformed",
            "result_profile_orientation_signs": [1, -1, 1],
            "result_self_intersection": True,
            "result_sweep_volume_m3": 0.0,
            "result_shape_owner": "stale/body2",
            "result_sweep_shape_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "swept_solids_use_current_frenet_frame_twist_transition_orientation_intersection_volume_owner_and_shape"
    ]


def test_v32_source_step_assembly_instance_transform_unit_color_material_uuid_digest_mismatch():
    row = _source_v32()
    row["replay_identity"][
        "step_assembly_instance_transform_unit_color_material_uuid_component_volume_file_generation_identity"
    ].update(
        {
            "instance_generation": "step-assembly-190",
            "material_generation": "step-assembly-189",
            "result_generation": "step-assembly-188",
            "decoded_instance_ids": ["shaft-1", "housing-1"],
            "decoded_instance_transform_sha256": [["housing-1", "f" * 64]],
            "decoded_length_unit": "mm",
            "decoded_instance_colors_rgb": [["housing-1", 1.0, 0.0, 0.0]],
            "decoded_material_labels": [["housing-1", "steel"]],
            "decoded_product_uuids": [["housing-1", "uuid-old"]],
            "decoded_component_shape_sha256": [["housing-1", "0" * 64]],
            "decoded_total_volume_m3": 1200.0,
            "decoded_step_file_sha256": "1" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_assemblies_use_current_instances_transforms_units_colors_materials_uuids_components_volume_and_file"
    ]


def test_v32_source_stl_repair_tolerance_normal_duplicate_vertex_watertight_volume_digest_mismatch():
    row = _source_v32()
    row["replay_identity"][
        "stl_repair_merge_tolerance_normal_duplicate_boundary_watertight_volume_unit_file_generation_identity"
    ].update(
        {
            "tolerance_generation": "stl-repair-190",
            "watertight_generation": "stl-repair-189",
            "result_generation": "stl-repair-188",
            "decoded_merge_tolerance_m": 1.0e-3,
            "decoded_facet_normal_sha256": "2" * 64,
            "decoded_duplicate_vertex_count": 12,
            "decoded_boundary_edge_count": 8,
            "decoded_watertight_component_count": 3,
            "decoded_repaired_volume_m3": 1.0e-3,
            "decoded_length_unit": "mm",
            "decoded_source_stl_sha256": "3" * 64,
            "decoded_repaired_stl_sha256": "4" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "stl_repairs_use_current_tolerance_normals_duplicates_boundaries_watertight_volume_unit_and_files"
    ]
