from __future__ import annotations

from test_build123d_generalization_v26 import (
    _public_result,
    _public_v26,
    _source_result,
    _source_v26,
)


def _public_v27():
    reference, measured = _public_v26()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["boolean_imprint_interface_owner_topology_name_tolerance_mass_generation_identity"] = {
                "boolean_generation": "boolean-imprint-141", "interface_boolean_generation": "boolean-imprint-141",
                "owner_boolean_generation": "boolean-imprint-141", "topology_boolean_generation": "boolean-imprint-141",
                "tolerance_boolean_generation": "boolean-imprint-141", "mass_boolean_generation": "boolean-imprint-141",
                "result_boolean_generation": "boolean-imprint-141", "operation": "imprint", "result_operation": "imprint",
                "interface_face_names": ["imprint-a-b"], "result_interface_face_names": ["imprint-a-b"],
                "interface_owner_pairs": [["body-a", "body-b"]],
                "result_interface_owner_pairs": [["body-a", "body-b"]],
                "persistent_topology_names": ["body-a", "body-b", "imprint-a-b"],
                "result_persistent_topology_names": ["body-a", "body-b", "imprint-a-b"],
                "linear_tolerance_m": 1.0e-7, "result_linear_tolerance_m": 1.0e-7,
                "solid_count": 2, "result_solid_count": 2,
                "total_volume_m3": 0.0015, "result_total_volume_m3": 0.0015,
                "interface_map_sha256": "3" * 64, "result_interface_map_sha256": "3" * 64,
                "mass_property_sha256": digest, "result_mass_property_sha256": digest,
            }
            row["loft_section_wire_seam_continuity_solid_volume_generation_identity"] = {
                "loft_generation": "loft-141", "section_loft_generation": "loft-141",
                "wire_loft_generation": "loft-141", "seam_loft_generation": "loft-141",
                "continuity_loft_generation": "loft-141", "solid_loft_generation": "loft-141",
                "volume_loft_generation": "loft-141", "result_loft_generation": "loft-141",
                "section_names": ["section-z0", "section-z1", "section-z2"],
                "result_section_names": ["section-z0", "section-z1", "section-z2"],
                "wire_orientation_signs": [1, 1, 1], "result_wire_orientation_signs": [1, 1, 1],
                "seam_vertex_names": ["seam-0", "seam-1", "seam-2"],
                "result_seam_vertex_names": ["seam-0", "seam-1", "seam-2"],
                "continuity": "C1", "result_continuity": "C1",
                "is_valid_solid": True, "result_is_valid_solid": True,
                "volume_m3": 0.0008, "result_volume_m3": 0.0008,
                "loft_shape_sha256": digest, "result_loft_shape_sha256": digest,
            }
    return reference, measured


def _source_v27():
    row = _source_v26()
    identity = row["replay_identity"]
    identity["step_import_unit_hierarchy_placement_color_shape_checksum_generation_identity"] = {
        "step_generation": "step-import-141", "unit_step_generation": "step-import-141",
        "hierarchy_step_generation": "step-import-141", "placement_step_generation": "step-import-141",
        "color_step_generation": "step-import-141", "shape_step_generation": "step-import-141",
        "result_step_generation": "step-import-141", "length_unit": "mm", "decoded_length_unit": "mm",
        "product_hierarchy": [["assembly", "base"], ["assembly", "cover"]],
        "decoded_product_hierarchy": [["assembly", "base"], ["assembly", "cover"]],
        "placements": [["base", 0.0, 0.0, 0.0], ["cover", 0.0, 0.0, 10.0]],
        "decoded_placements": [["base", 0.0, 0.0, 0.0], ["cover", 0.0, 0.0, 10.0]],
        "colors_rgb": [["base", 0.8, 0.1, 0.1], ["cover", 0.1, 0.1, 0.8]],
        "decoded_colors_rgb": [["base", 0.8, 0.1, 0.1], ["cover", 0.1, 0.1, 0.8]],
        "shape_ids": ["base", "cover"], "decoded_shape_ids": ["base", "cover"],
        "step_sha256": "4" * 64, "decoded_step_sha256": "4" * 64,
        "shape_map_sha256": "5" * 64, "decoded_shape_map_sha256": "5" * 64,
    }
    identity["brep_serialization_kernel_tolerance_location_cache_generation_identity"] = {
        "brep_generation": "brep-cache-141", "kernel_brep_generation": "brep-cache-141",
        "tolerance_brep_generation": "brep-cache-141", "location_brep_generation": "brep-cache-141",
        "shape_brep_generation": "brep-cache-141", "cache_brep_generation": "brep-cache-141",
        "result_brep_generation": "brep-cache-141", "kernel_version": "OCCT-7.8.1",
        "decoded_kernel_version": "OCCT-7.8.1", "modeling_tolerance_m": 1.0e-7,
        "decoded_modeling_tolerance_m": 1.0e-7,
        "location_transform": [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03]],
        "decoded_location_transform": [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03]],
        "shape_generation": "shape-141", "decoded_shape_generation": "shape-141",
        "shape_sha256": "6" * 64, "decoded_shape_sha256": "6" * 64,
        "brep_cache_sha256": "7" * 64, "decoded_brep_cache_sha256": "7" * 64,
    }
    return row


def test_v27_positive_public_and_source_identity():
    reference, measured = _public_v27()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v27())["status"] == "ok"


def test_v27_public_boolean_imprint_interface_owner_topology_name_tolerance_mass_property_mismatch():
    reference, measured = _public_v27()
    measured["external_cad"][0]["boolean_imprint_interface_owner_topology_name_tolerance_mass_generation_identity"].update(
        {"interface_boolean_generation": "boolean-imprint-140", "result_operation": "cut",
         "result_interface_face_names": ["stale-face"],
         "result_interface_owner_pairs": [["body-a", "body-c"]],
         "result_linear_tolerance_m": 1.0e-3, "result_total_volume_m3": 0.0012,
         "result_mass_property_sha256": "9" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["boolean_imprints_use_current_interfaces_owners_topology_names_tolerance_and_mass"]


def test_v27_public_loft_section_order_wire_orientation_seam_continuity_solid_volume_mismatch():
    reference, measured = _public_v27()
    measured["external_cad"][0]["loft_section_wire_seam_continuity_solid_volume_generation_identity"].update(
        {"section_loft_generation": "loft-140",
         "result_section_names": ["section-z0", "section-z2", "section-z1"],
         "result_wire_orientation_signs": [1, -1, 1], "result_continuity": "C0",
         "result_is_valid_solid": False, "result_volume_m3": 0.0006}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["lofts_use_current_sections_wire_orientation_seams_continuity_solid_and_volume"]


def test_v27_source_step_import_unit_product_hierarchy_placement_color_shape_checksum_mismatch():
    row = _source_v27()
    row["replay_identity"]["step_import_unit_hierarchy_placement_color_shape_checksum_generation_identity"].update(
        {"unit_step_generation": "step-import-140", "decoded_length_unit": "m",
         "decoded_product_hierarchy": [["assembly", "cover"], ["assembly", "base-old"]],
         "decoded_shape_ids": ["cover", "base-old"], "decoded_step_sha256": "b" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["step_import_uses_current_units_hierarchy_placements_colors_shapes_and_checksums"]


def test_v27_source_brep_serialization_kernel_version_tolerance_location_cache_digest_mismatch():
    row = _source_v27()
    row["replay_identity"]["brep_serialization_kernel_tolerance_location_cache_generation_identity"].update(
        {"kernel_brep_generation": "brep-cache-140", "decoded_kernel_version": "OCCT-7.7.0",
         "decoded_modeling_tolerance_m": 1.0e-4, "decoded_shape_generation": "shape-140",
         "decoded_brep_cache_sha256": "e" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["brep_cache_uses_current_kernel_tolerance_location_shape_and_digest"]
