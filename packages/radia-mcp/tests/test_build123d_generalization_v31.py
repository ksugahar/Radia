from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v30 import _public_v30, _source_v30


_PROMOTED_CASE_IDS = (
    "v31_public_assembly_occurrence_location_density_unit_mass_inertia_parallel_axis_mismatch",
    "v31_public_loft_sweep_topology_seam_face_lineage_orientation_volume_mismatch",
    "v31_source_step_ap242_occurrence_transform_length_unit_product_structure_digest_mismatch",
    "v31_source_stl_tessellation_chord_angle_normal_watertight_volume_tolerance_mismatch",
)


def _public_v31():
    reference, measured = _public_v30()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            generation = "assembly-mass-181"
            shape_digest = ("1" if index == 0 else "2") * 64
            locations = [
                [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
                [[1.0, 0.0, 0.0, 0.1], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            ]
            inertia = [
                [0.01, 0.0, 0.0],
                [0.0, 0.0205, 0.0],
                [0.0, 0.0, 0.0205],
            ]
            row[
                "assembly_occurrence_location_density_unit_suppression_mass_center_inertia_parallel_axis_shape_generation_identity"
            ] = {
                "assembly_generation": generation,
                "occurrence_assembly_generation": generation,
                "location_assembly_generation": generation,
                "density_assembly_generation": generation,
                "suppression_assembly_generation": generation,
                "mass_assembly_generation": generation,
                "inertia_assembly_generation": generation,
                "shape_assembly_generation": generation,
                "result_assembly_generation": generation,
                "occurrence_ids": ["housing-1", "shaft-1"],
                "result_occurrence_ids": ["housing-1", "shaft-1"],
                "location_matrices": locations,
                "result_location_matrices": locations,
                "densities_kg_m3": [2700.0, 7800.0],
                "result_densities_kg_m3": [2700.0, 7800.0],
                "density_unit": "kg/m^3",
                "result_density_unit": "kg/m^3",
                "suppressed_occurrence_ids": [],
                "result_suppressed_occurrence_ids": [],
                "part_volumes_m3": [1.0e-3, 2.0e-4],
                "result_part_volumes_m3": [1.0e-3, 2.0e-4],
                "assembly_mass_kg": 4.26,
                "result_assembly_mass_kg": 4.26,
                "center_of_mass_m": [0.03661971830985916, 0.0, 0.0],
                "result_center_of_mass_m": [0.03661971830985916, 0.0, 0.0],
                "inertia_reference_frame": "assembly_global_center_of_mass",
                "result_inertia_reference_frame": "assembly_global_center_of_mass",
                "parallel_axis_applied": True,
                "result_parallel_axis_applied": True,
                "assembly_inertia_kg_m2": inertia,
                "result_assembly_inertia_kg_m2": inertia,
                "assembly_shape_sha256": shape_digest,
                "result_assembly_shape_sha256": shape_digest,
            }
            generation = "loft-lineage-181"
            lineage_digest = ("3" if index == 0 else "4") * 64
            loft_digest = ("5" if index == 0 else "6") * 64
            lineage = [
                ["section-0:e1", "loft:f1"],
                ["section-1:e1", "loft:f1"],
                ["section-2:e1", "loft:f1"],
            ]
            row[
                "loft_sweep_profile_order_seam_guide_orientation_face_lineage_shell_volume_shape_generation_identity"
            ] = {
                "loft_generation": generation,
                "profile_loft_generation": generation,
                "seam_loft_generation": generation,
                "guide_loft_generation": generation,
                "lineage_loft_generation": generation,
                "shell_loft_generation": generation,
                "mass_loft_generation": generation,
                "shape_loft_generation": generation,
                "result_loft_generation": generation,
                "profile_ids": ["section-0", "section-1", "section-2"],
                "result_profile_ids": ["section-0", "section-1", "section-2"],
                "profile_parameters": [0.0, 0.5, 1.0],
                "result_profile_parameters": [0.0, 0.5, 1.0],
                "seam_parameters": [0.0, 0.0, 0.0],
                "result_seam_parameters": [0.0, 0.0, 0.0],
                "guide_orientation": "start_to_end_right_handed",
                "result_guide_orientation": "start_to_end_right_handed",
                "face_lineage_pairs": lineage,
                "result_face_lineage_pairs": lineage,
                "face_lineage_sha256": lineage_digest,
                "result_face_lineage_sha256": lineage_digest,
                "shell_closed": True,
                "result_shell_closed": True,
                "volume_m3": 1.5e-3,
                "result_volume_m3": 1.5e-3,
                "loft_shape_sha256": loft_digest,
                "result_loft_shape_sha256": loft_digest,
            }
    return reference, measured


def _source_v31():
    row = _source_v30()
    identity = row["replay_identity"]
    generation = "step-ap242-context-181"
    identity[
        "step_ap242_representation_context_external_owner_occurrence_transform_mass_product_structure_file_generation_identity"
    ] = {
        "step_generation": generation,
        "schema_step_generation": generation,
        "context_step_generation": generation,
        "owner_step_generation": generation,
        "occurrence_step_generation": generation,
        "transform_step_generation": generation,
        "mass_step_generation": generation,
        "structure_step_generation": generation,
        "file_step_generation": generation,
        "result_step_generation": generation,
        "ap_schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF",
        "decoded_ap_schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF",
        "representation_contexts": [["ctx-mm", "mm", "degree"], ["ctx-m", "m", "radian"]],
        "decoded_representation_contexts": [["ctx-mm", "mm", "degree"], ["ctx-m", "m", "radian"]],
        "part_owners": [["housing", "part-housing-v3"], ["shaft", "part-shaft-v5"]],
        "decoded_part_owners": [["housing", "part-housing-v3"], ["shaft", "part-shaft-v5"]],
        "occurrence_ids": ["housing-1", "shaft-1", "shaft-2"],
        "decoded_occurrence_ids": ["housing-1", "shaft-1", "shaft-2"],
        "occurrence_transform_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64], ["shaft-2", "9" * 64]],
        "decoded_occurrence_transform_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64], ["shaft-2", "9" * 64]],
        "occurrence_mass_properties": [["housing-1", 2.7, 0.0, 0.0, 0.0], ["shaft-1", 0.78, 0.1, 0.0, 0.0], ["shaft-2", 0.78, -0.1, 0.0, 0.0]],
        "decoded_occurrence_mass_properties": [["housing-1", 2.7, 0.0, 0.0, 0.0], ["shaft-1", 0.78, 0.1, 0.0, 0.0], ["shaft-2", 0.78, -0.1, 0.0, 0.0]],
        "product_structure_sha256": "a" * 64,
        "decoded_product_structure_sha256": "a" * 64,
        "step_file_sha256": "b" * 64,
        "decoded_step_file_sha256": "b" * 64,
    }
    generation = "stl-tessellation-error-181"
    identity[
        "stl_tessellation_source_brep_chord_angle_facet_component_deviation_area_volume_digest_generation_identity"
    ] = {
        "stl_generation": generation,
        "source_stl_generation": generation,
        "tolerance_stl_generation": generation,
        "facet_stl_generation": generation,
        "component_stl_generation": generation,
        "deviation_stl_generation": generation,
        "area_stl_generation": generation,
        "volume_stl_generation": generation,
        "file_stl_generation": generation,
        "result_stl_generation": generation,
        "source_brep_sha256": "c" * 64,
        "decoded_source_brep_sha256": "c" * 64,
        "chord_tolerance_m": 1.0e-4,
        "decoded_chord_tolerance_m": 1.0e-4,
        "angular_tolerance_deg": 10.0,
        "decoded_angular_tolerance_deg": 10.0,
        "facet_count": 1200,
        "decoded_facet_count": 1200,
        "connected_component_count": 1,
        "decoded_connected_component_count": 1,
        "maximum_surface_deviation_m": 8.0e-5,
        "decoded_maximum_surface_deviation_m": 8.0e-5,
        "source_surface_area_m2": 0.02,
        "decoded_surface_area_m2": 0.01999,
        "source_volume_m3": 1.5e-3,
        "decoded_volume_m3": 1.4995e-3,
        "relative_area_tolerance": 1.0e-3,
        "relative_volume_tolerance": 1.0e-3,
        "normal_table_sha256": "d" * 64,
        "decoded_normal_table_sha256": "d" * 64,
        "stl_file_sha256": "e" * 64,
        "decoded_stl_file_sha256": "e" * 64,
    }
    return row


def test_v31_positive_public_and_source_identity():
    reference, measured = _public_v31()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v31())["status"] == "ok"


def test_v31_public_assembly_occurrence_location_density_unit_mass_inertia_parallel_axis_mismatch():
    reference, measured = _public_v31()
    measured["external_cad"][0][
        "assembly_occurrence_location_density_unit_suppression_mass_center_inertia_parallel_axis_shape_generation_identity"
    ].update(
        {
            "location_assembly_generation": "assembly-mass-180",
            "density_assembly_generation": "assembly-mass-179",
            "result_assembly_generation": "assembly-mass-178",
            "result_occurrence_ids": ["shaft-1", "housing-1"],
            "result_location_matrices": [[[1.0, 0.0, 0.0, 0.0]]],
            "result_densities_kg_m3": [2.7, 7.8],
            "result_density_unit": "g/cm^3",
            "result_suppressed_occurrence_ids": ["housing-1"],
            "result_part_volumes_m3": [2.0e-4],
            "result_assembly_mass_kg": 1.56,
            "result_center_of_mass_m": [0.1, 0.0, 0.0],
            "result_inertia_reference_frame": "shaft-local-origin",
            "result_parallel_axis_applied": False,
            "result_assembly_inertia_kg_m2": [[0.01]],
            "result_assembly_shape_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "assembly_occurrences_use_current_locations_densities_units_suppression_mass_center_parallel_axis_inertia_and_shape"
    ]


def test_v31_public_loft_sweep_topology_seam_face_lineage_orientation_volume_mismatch():
    reference, measured = _public_v31()
    measured["external_cad"][0][
        "loft_sweep_profile_order_seam_guide_orientation_face_lineage_shell_volume_shape_generation_identity"
    ].update(
        {
            "profile_loft_generation": "loft-lineage-180",
            "lineage_loft_generation": "loft-lineage-179",
            "result_loft_generation": "loft-lineage-178",
            "result_profile_ids": ["section-2", "section-1", "section-0"],
            "result_profile_parameters": [0.0, 0.8, 0.5],
            "result_seam_parameters": [0.5, 0.0, 0.25],
            "result_guide_orientation": "end_to_start_left_handed",
            "result_face_lineage_pairs": [["section-0:e2", "loft:f2"]],
            "result_face_lineage_sha256": "0" * 64,
            "result_shell_closed": False,
            "result_volume_m3": 0.0,
            "result_loft_shape_sha256": "1" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "lofts_use_current_profile_order_seams_guide_face_lineage_shell_volume_and_shape"
    ]


def test_v31_source_step_ap242_occurrence_transform_length_unit_product_structure_digest_mismatch():
    row = _source_v31()
    row["replay_identity"][
        "step_ap242_representation_context_external_owner_occurrence_transform_mass_product_structure_file_generation_identity"
    ].update(
        {
            "context_step_generation": "step-ap242-context-180",
            "owner_step_generation": "step-ap242-context-179",
            "result_step_generation": "step-ap242-context-178",
            "decoded_ap_schema": "AP203",
            "decoded_representation_contexts": [["ctx-inch", "inch", "degree"]],
            "decoded_part_owners": [["housing", "part-old"]],
            "decoded_occurrence_ids": ["shaft-2", "housing-1"],
            "decoded_occurrence_transform_sha256": [["housing-1", "2" * 64]],
            "decoded_occurrence_mass_properties": [["housing-1", 2700.0, 0.0, 0.0, 0.0]],
            "decoded_product_structure_sha256": "3" * 64,
            "decoded_step_file_sha256": "4" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_ap242_occurrences_use_current_schema_context_owners_transforms_mass_structure_and_file"
    ]


def test_v31_source_stl_tessellation_chord_angle_normal_watertight_volume_tolerance_mismatch():
    row = _source_v31()
    row["replay_identity"][
        "stl_tessellation_source_brep_chord_angle_facet_component_deviation_area_volume_digest_generation_identity"
    ].update(
        {
            "source_stl_generation": "stl-tessellation-error-180",
            "deviation_stl_generation": "stl-tessellation-error-179",
            "result_stl_generation": "stl-tessellation-error-178",
            "decoded_source_brep_sha256": "5" * 64,
            "decoded_chord_tolerance_m": 1.0e-2,
            "decoded_angular_tolerance_deg": 45.0,
            "decoded_facet_count": 120,
            "decoded_connected_component_count": 3,
            "decoded_maximum_surface_deviation_m": 2.0e-2,
            "decoded_surface_area_m2": 0.03,
            "decoded_volume_m3": 1.0e-3,
            "decoded_normal_table_sha256": "6" * 64,
            "decoded_stl_file_sha256": "7" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "stl_tessellations_use_current_brep_chord_angle_facets_components_deviation_area_volume_and_digests"
    ]
