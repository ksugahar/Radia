from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v32 import _public_v32, _source_v32


_PROMOTED_CASE_IDS = (
    "v33_public_loft_section_guide_parameterization_seam_orientation_smooth_volume_mismatch",
    "v33_public_mass_property_inertia_origin_density_unit_principal_axis_parallel_axis_mismatch",
    "v33_source_brep_face_pcurve_orientation_location_tolerance_serialization_digest_mismatch",
    "v33_source_gltf_node_transform_triangle_material_unit_tessellation_volume_digest_mismatch",
)


def _public_v33():
    reference, measured = _public_v32()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            suffix = "1" if index == 0 else "2"
            generation = "guided-loft-201"
            row[
                "loft_section_guide_parameterization_seam_orientation_mode_intersection_volume_shape_generation_identity"
            ] = {
                "loft_generation": generation,
                "section_generation": generation,
                "parameter_generation": generation,
                "guide_generation": generation,
                "seam_generation": generation,
                "mode_generation": generation,
                "intersection_generation": generation,
                "mass_generation": generation,
                "shape_generation": generation,
                "result_generation": generation,
                "section_order": ["section-0", "section-1", "section-2"],
                "result_section_order": ["section-0", "section-1", "section-2"],
                "wire_parameterization_sha256": "3" * 64,
                "result_wire_parameterization_sha256": "3" * 64,
                "guide_intersections": [
                    ["guide-0", "section-0"],
                    ["guide-0", "section-1"],
                    ["guide-0", "section-2"],
                ],
                "result_guide_intersections": [
                    ["guide-0", "section-0"],
                    ["guide-0", "section-1"],
                    ["guide-0", "section-2"],
                ],
                "seam_orientation_signs": [1, 1, 1],
                "result_seam_orientation_signs": [1, 1, 1],
                "loft_mode": "smooth",
                "result_loft_mode": "smooth",
                "self_intersection": False,
                "result_self_intersection": False,
                "loft_volume_m3": 8.0e-4,
                "result_loft_volume_m3": 8.0e-4,
                "loft_shape_sha256": suffix * 64,
                "result_loft_shape_sha256": suffix * 64,
            }
            generation = "mass-inertia-201"
            mass_suffix = "4" if index == 0 else "5"
            row[
                "mass_property_density_unit_origin_center_principal_axis_degeneracy_parallel_axis_owner_shape_generation_identity"
            ] = {
                "mass_property_generation": generation,
                "density_generation": generation,
                "origin_generation": generation,
                "center_generation": generation,
                "principal_generation": generation,
                "axis_generation": generation,
                "parallel_axis_generation": generation,
                "owner_generation": generation,
                "shape_generation": generation,
                "result_generation": generation,
                "density_kg_m3": 7800.0,
                "result_density_kg_m3": 7800.0,
                "density_unit": "kg/m^3",
                "result_density_unit": "kg/m^3",
                "mass_kg": 7.8,
                "result_mass_kg": 7.8,
                "inertia_origin_m": [0.0, 0.0, 0.0],
                "result_inertia_origin_m": [0.0, 0.0, 0.0],
                "center_of_mass_m": [0.1, 0.0, 0.0],
                "result_center_of_mass_m": [0.1, 0.0, 0.0],
                "inertia_at_com_kg_m2": [
                    [0.02, 0.0, 0.0],
                    [0.0, 0.03, 0.0],
                    [0.0, 0.0, 0.04],
                ],
                "result_inertia_at_com_kg_m2": [
                    [0.02, 0.0, 0.0],
                    [0.0, 0.03, 0.0],
                    [0.0, 0.0, 0.04],
                ],
                "inertia_at_origin_kg_m2": [
                    [0.02, 0.0, 0.0],
                    [0.0, 0.108, 0.0],
                    [0.0, 0.0, 0.118],
                ],
                "result_inertia_at_origin_kg_m2": [
                    [0.02, 0.0, 0.0],
                    [0.0, 0.108, 0.0],
                    [0.0, 0.0, 0.118],
                ],
                "principal_moments_kg_m2": [0.02, 0.03, 0.04],
                "result_principal_moments_kg_m2": [0.02, 0.03, 0.04],
                "principal_axes": [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                "result_principal_axes": [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                "degeneracy_convention": "right_handed_sorted_moments",
                "result_degeneracy_convention": "right_handed_sorted_moments",
                "shape_owner": "solid/body-1",
                "result_shape_owner": "solid/body-1",
                "mass_shape_sha256": mass_suffix * 64,
                "result_mass_shape_sha256": mass_suffix * 64,
            }
    return reference, measured


def _source_v33():
    row = _source_v32()
    identity = row["replay_identity"]
    generation = "brep-semantic-roundtrip-201"
    identity[
        "brep_face_pcurve_wire_orientation_location_tolerance_surface_serializer_shape_generation_identity"
    ] = {
        "brep_generation": generation,
        "pcurve_generation": generation,
        "wire_generation": generation,
        "location_generation": generation,
        "tolerance_generation": generation,
        "surface_generation": generation,
        "serializer_generation": generation,
        "shape_generation": generation,
        "result_generation": generation,
        "face_pcurve_sha256": [[1, "6" * 64], [2, "7" * 64]],
        "decoded_face_pcurve_sha256": [[1, "6" * 64], [2, "7" * 64]],
        "wire_orientation_signs": [[1, 1], [2, -1]],
        "decoded_wire_orientation_signs": [[1, 1], [2, -1]],
        "nested_location_sha256": [
            ["root/part-1", "8" * 64],
            ["root/part-1/face-2", "9" * 64],
        ],
        "decoded_nested_location_sha256": [
            ["root/part-1", "8" * 64],
            ["root/part-1/face-2", "9" * 64],
        ],
        "edge_tolerances_m": [[11, 1.0e-7], [12, 2.0e-7]],
        "decoded_edge_tolerances_m": [[11, 1.0e-7], [12, 2.0e-7]],
        "surface_types": [[1, "plane"], [2, "cylinder"]],
        "decoded_surface_types": [[1, "plane"], [2, "cylinder"]],
        "serializer_version": "occt-brep-v3",
        "decoded_serializer_version": "occt-brep-v3",
        "brep_shape_sha256": "a" * 64,
        "decoded_brep_shape_sha256": "a" * 64,
    }
    generation = "gltf-roundtrip-201"
    identity[
        "gltf_node_hierarchy_transform_winding_material_unit_tessellation_volume_file_generation_identity"
    ] = {
        "gltf_generation": generation,
        "node_generation": generation,
        "transform_generation": generation,
        "winding_generation": generation,
        "material_generation": generation,
        "unit_generation": generation,
        "tessellation_generation": generation,
        "volume_generation": generation,
        "file_generation": generation,
        "result_generation": generation,
        "node_hierarchy": [["root", "part-1"], ["part-1", "mesh-1"]],
        "decoded_node_hierarchy": [["root", "part-1"], ["part-1", "mesh-1"]],
        "instance_transform_sha256": [["part-1", "b" * 64]],
        "decoded_instance_transform_sha256": [["part-1", "b" * 64]],
        "triangle_winding_sha256": "c" * 64,
        "decoded_triangle_winding_sha256": "c" * 64,
        "material_assignments": [["mesh-1", "steel-painted"]],
        "decoded_material_assignments": [["mesh-1", "steel-painted"]],
        "length_unit": "m",
        "decoded_length_unit": "m",
        "linear_deflection_m": 1.0e-4,
        "decoded_linear_deflection_m": 1.0e-4,
        "angular_deflection_rad": 0.1,
        "decoded_angular_deflection_rad": 0.1,
        "triangle_count": 512,
        "decoded_triangle_count": 512,
        "enclosed_volume_m3": 1.2e-3,
        "decoded_enclosed_volume_m3": 1.2e-3,
        "gltf_file_sha256": "d" * 64,
        "decoded_gltf_file_sha256": "d" * 64,
    }
    return row


def test_v33_positive_public_and_source_identity():
    reference, measured = _public_v33()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v33())["status"] == "ok"


def test_v33_public_loft_section_guide_parameterization_seam_orientation_smooth_volume_mismatch():
    reference, measured = _public_v33()
    measured["external_cad"][0][
        "loft_section_guide_parameterization_seam_orientation_mode_intersection_volume_shape_generation_identity"
    ].update(
        {
            "section_generation": "guided-loft-200",
            "guide_generation": "guided-loft-199",
            "result_generation": "guided-loft-198",
            "result_section_order": ["section-0", "section-2", "section-1"],
            "result_wire_parameterization_sha256": "e" * 64,
            "result_guide_intersections": [["guide-0", "section-0"]],
            "result_seam_orientation_signs": [1, -1, 1],
            "result_loft_mode": "ruled",
            "result_self_intersection": True,
            "result_loft_volume_m3": 0.0,
            "result_loft_shape_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "guided_lofts_use_current_sections_parameterization_guides_seams_mode_intersection_volume_and_shape"
    ]


def test_v33_public_mass_property_inertia_origin_density_unit_principal_axis_parallel_axis_mismatch():
    reference, measured = _public_v33()
    measured["external_cad"][0][
        "mass_property_density_unit_origin_center_principal_axis_degeneracy_parallel_axis_owner_shape_generation_identity"
    ].update(
        {
            "density_generation": "mass-inertia-200",
            "parallel_axis_generation": "mass-inertia-199",
            "result_generation": "mass-inertia-198",
            "result_density_kg_m3": 7.8,
            "result_density_unit": "g/cm^3",
            "result_mass_kg": 0.0078,
            "result_inertia_origin_m": [0.1, 0.0, 0.0],
            "result_center_of_mass_m": [0.0, 0.0, 0.0],
            "result_principal_moments_kg_m2": [0.04, 0.03, 0.02],
            "result_principal_axes": [
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            "result_degeneracy_convention": "unsorted_left_handed",
            "result_shape_owner": "stale/body-2",
            "result_mass_shape_sha256": "0" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mass_properties_use_current_density_units_origin_center_principal_axes_parallel_axis_owner_and_shape"
    ]


def test_v33_source_brep_face_pcurve_orientation_location_tolerance_serialization_digest_mismatch():
    row = _source_v33()
    row["replay_identity"][
        "brep_face_pcurve_wire_orientation_location_tolerance_surface_serializer_shape_generation_identity"
    ].update(
        {
            "pcurve_generation": "brep-semantic-roundtrip-200",
            "location_generation": "brep-semantic-roundtrip-199",
            "result_generation": "brep-semantic-roundtrip-198",
            "decoded_face_pcurve_sha256": [[1, "1" * 64]],
            "decoded_wire_orientation_signs": [[1, -1], [2, 1]],
            "decoded_nested_location_sha256": [["root/part-1", "2" * 64]],
            "decoded_edge_tolerances_m": [[11, 1.0e-3]],
            "decoded_surface_types": [[1, "bspline"], [2, "plane"]],
            "decoded_serializer_version": "occt-brep-v2",
            "decoded_brep_shape_sha256": "3" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "brep_roundtrips_use_current_pcurves_wire_orientations_locations_tolerances_surfaces_serializer_and_shape"
    ]


def test_v33_source_gltf_node_transform_triangle_material_unit_tessellation_volume_digest_mismatch():
    row = _source_v33()
    row["replay_identity"][
        "gltf_node_hierarchy_transform_winding_material_unit_tessellation_volume_file_generation_identity"
    ].update(
        {
            "node_generation": "gltf-roundtrip-200",
            "tessellation_generation": "gltf-roundtrip-199",
            "result_generation": "gltf-roundtrip-198",
            "decoded_node_hierarchy": [["root", "mesh-1"]],
            "decoded_instance_transform_sha256": [["part-1", "4" * 64]],
            "decoded_triangle_winding_sha256": "5" * 64,
            "decoded_material_assignments": [["mesh-1", "default"]],
            "decoded_length_unit": "mm",
            "decoded_linear_deflection_m": 1.0e-2,
            "decoded_angular_deflection_rad": 0.5,
            "decoded_triangle_count": 128,
            "decoded_enclosed_volume_m3": 1.2e6,
            "decoded_gltf_file_sha256": "6" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "gltf_roundtrips_use_current_hierarchy_transforms_winding_materials_units_tessellation_volume_and_file"
    ]
