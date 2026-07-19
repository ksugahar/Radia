from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v39 import _generation
from test_build123d_generalization_v40 import _public_v40, _source_v40


_LOFT = (
    "loft_section_order_parameterization_seam_twist_closure_topology_volume_"
    "centroid_owner_brep_generation_identity"
)
_SHELL = (
    "shell_offset_removedface_direction_thickness_join_selfintersection_"
    "validity_mass_owner_brep_generation_identity"
)
_ASSEMBLY = (
    "assembly_interference_part_transform_unit_contact_overlap_clearance_"
    "owner_result_generation_identity"
)
_STEP = (
    "step_name_color_layer_hierarchy_subshape_label_source_import_owner_"
    "result_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v41_public_loft_section_order_parameterization_seam_twist_solid_volume_centroid_mismatch",
    "v41_public_shell_offset_thickness_joinmode_selfintersection_massproperty_brep_mismatch",
    "v41_source_assembly_interference_contactpair_overlapvolume_transform_unit_owner_mismatch",
    "v41_source_step_name_color_layer_metadata_subshape_label_roundtrip_owner_mismatch",
)


def _public_v41():
    reference, measured = _public_v40()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "loft-solid-724"
            row[_LOFT] = {
                "loft_generation": generation,
                **_generation(
                    generation,
                    "section_generation",
                    "parameter_generation",
                    "seam_generation",
                    "twist_generation",
                    "closure_generation",
                    "topology_generation",
                    "mass_generation",
                    "owner_generation",
                    "brep_generation",
                    "result_generation",
                ),
                "section_ids": [1, 2, 3],
                "result_section_ids": [1, 2, 3],
                "section_parameters": [0.0, 0.5, 1.0],
                "result_section_parameters": [0.0, 0.5, 1.0],
                "seam_vertex_ids": [10, 20, 30],
                "result_seam_vertex_ids": [10, 20, 30],
                "section_twist_deg": [0.0, 10.0, 20.0],
                "result_section_twist_deg": [0.0, 10.0, 20.0],
                "closed_profile": True,
                "result_closed_profile": True,
                "solid_valid": True,
                "result_solid_valid": True,
                "topology_signature": {
                    "solid": 1,
                    "shell": 1,
                    "face": 8,
                    "edge": 18,
                    "vertex": 12,
                },
                "result_topology_signature": {
                    "solid": 1,
                    "shell": 1,
                    "face": 8,
                    "edge": 18,
                    "vertex": 12,
                },
                "volume_m3": 2.5e-3,
                "result_volume_m3": 2.5e-3,
                "centroid_m": [0.0, 0.0, 0.05],
                "result_centroid_m": [0.0, 0.0, 0.05],
                "shape_owner": "part:loft-solid-724",
                "result_shape_owner": "part:loft-solid-724",
                "loft_brep_sha256": suffix * 64,
                "accepted_loft_brep_sha256": suffix * 64,
            }

            generation = "shell-offset-724"
            row[_SHELL] = {
                "shell_generation": generation,
                **_generation(
                    generation,
                    "face_generation",
                    "offset_generation",
                    "thickness_generation",
                    "join_generation",
                    "intersection_generation",
                    "validity_generation",
                    "mass_generation",
                    "owner_generation",
                    "brep_generation",
                    "result_generation",
                ),
                "removed_face_ids": [1],
                "result_removed_face_ids": [1],
                "offset_direction": "inward",
                "result_offset_direction": "inward",
                "signed_offset_m": -2.0e-3,
                "result_signed_offset_m": -2.0e-3,
                "wall_thickness_m": 2.0e-3,
                "result_wall_thickness_m": 2.0e-3,
                "join_mode": "arc",
                "result_join_mode": "arc",
                "self_intersection_free": True,
                "result_self_intersection_free": True,
                "solid_valid": True,
                "result_solid_valid": True,
                "volume_m3": 4.0e-4,
                "result_volume_m3": 4.0e-4,
                "surface_area_m2": 0.22,
                "result_surface_area_m2": 0.22,
                "centroid_m": [0.0, 0.0, 0.04],
                "result_centroid_m": [0.0, 0.0, 0.04],
                "shape_owner": "part:shell-offset-724",
                "result_shape_owner": "part:shell-offset-724",
                "shell_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_shell_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v41():
    row = _source_v40()
    identity = row["replay_identity"]
    generation = "assembly-interference-724"
    base_transform = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
    slider_transform = [
        [1, 0, 0, 0.1],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
    identity[_ASSEMBLY] = {
        "interference_generation": generation,
        **_generation(
            generation,
            "part_generation",
            "transform_generation",
            "unit_generation",
            "contact_generation",
            "overlap_generation",
            "clearance_generation",
            "owner_generation",
            "result_generation",
        ),
        "part_names": ["base", "slider"],
        "replayed_part_names": ["base", "slider"],
        "part_transforms": {"base": base_transform, "slider": slider_transform},
        "replayed_part_transforms": {"base": base_transform, "slider": slider_transform},
        "length_unit": "m",
        "replayed_length_unit": "m",
        "contact_pairs": [["base:face1", "slider:face2"]],
        "replayed_contact_pairs": [["base:face1", "slider:face2"]],
        "overlap_volumes_m3": [0.0],
        "replayed_overlap_volumes_m3": [0.0],
        "minimum_clearance_m": 1.0e-3,
        "replayed_minimum_clearance_m": 1.0e-3,
        "assembly_generation_id": 724,
        "replayed_assembly_generation_id": 724,
        "shape_owners": {"base": "headless:base-724", "slider": "headless:slider-724"},
        "replayed_shape_owners": {
            "base": "headless:base-724",
            "slider": "headless:slider-724",
        },
        "assembly_result_sha256": "5" * 64,
        "accepted_assembly_result_sha256": "5" * 64,
    }

    generation = "step-metadata-724"
    identity[_STEP] = {
        "step_generation": generation,
        **_generation(
            generation,
            "name_generation",
            "color_generation",
            "layer_generation",
            "hierarchy_generation",
            "label_generation",
            "source_generation",
            "import_generation",
            "result_generation",
        ),
        "product_names": {"part:1": "rotor", "part:2": "shaft"},
        "replayed_product_names": {"part:1": "rotor", "part:2": "shaft"},
        "colors_rgb": {"part:1": [1.0, 0.0, 0.0], "part:2": [0.5, 0.5, 0.5]},
        "replayed_colors_rgb": {
            "part:1": [1.0, 0.0, 0.0],
            "part:2": [0.5, 0.5, 0.5],
        },
        "layers": {"part:1": "magnet", "part:2": "mechanical"},
        "replayed_layers": {"part:1": "magnet", "part:2": "mechanical"},
        "assembly_hierarchy": {"assembly:1": ["part:1", "part:2"]},
        "replayed_assembly_hierarchy": {"assembly:1": ["part:1", "part:2"]},
        "subshape_labels": {"part:1/face:1": "airgap", "part:2/face:2": "bearing"},
        "replayed_subshape_labels": {
            "part:1/face:1": "airgap",
            "part:2/face:2": "bearing",
        },
        "source_shape_owner": "headless:step-source-724",
        "replayed_source_shape_owner": "headless:step-source-724",
        "imported_shape_owner": "headless:step-import-724",
        "replayed_imported_shape_owner": "headless:step-import-724",
        "source_brep_sha256": "6" * 64,
        "replayed_source_brep_sha256": "6" * 64,
        "step_file_sha256": "7" * 64,
        "replayed_step_file_sha256": "7" * 64,
        "metadata_result_sha256": "8" * 64,
        "accepted_metadata_result_sha256": "8" * 64,
    }
    return row


def test_v41_positive_contracts():
    reference, measured = _public_v41()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v41())["status"] == "ok"


def test_v41_public_loft_mismatch():
    reference, measured = _public_v41()
    measured["external_cad"][0][_LOFT].update(
        {
            "section_generation": "loft-solid-723",
            "result_section_ids": [3, 2, 1],
            "result_section_parameters": [0.0, 0.8, 0.7],
            "result_seam_vertex_ids": [30, 10, 20],
            "result_section_twist_deg": [0.0, 200.0],
            "result_closed_profile": False,
            "result_solid_valid": False,
            "result_volume_m3": -1.0,
            "result_shape_owner": "stale:loft",
            "accepted_loft_brep_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "lofts_use_current_sections_parameters_seams_twist_closure_topology_mass_owner_and_brep"
    ]


def test_v41_public_shell_mismatch():
    reference, measured = _public_v41()
    measured["external_cad"][0][_SHELL].update(
        {
            "offset_generation": "shell-offset-723",
            "result_removed_face_ids": [99],
            "result_offset_direction": "outward",
            "result_signed_offset_m": 2.0e-3,
            "result_wall_thickness_m": -2.0e-3,
            "result_join_mode": "bad",
            "result_self_intersection_free": False,
            "result_solid_valid": False,
            "result_volume_m3": -1.0,
            "result_shape_owner": "stale:shell",
            "accepted_shell_brep_sha256": "b" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "shells_use_current_faces_offset_thickness_join_intersection_validity_mass_owner_and_brep"
    ]


def test_v41_source_assembly_interference_mismatch():
    row = _source_v41()
    row["replay_identity"][_ASSEMBLY].update(
        {
            "transform_generation": "assembly-interference-723",
            "replayed_part_names": ["old"],
            "replayed_part_transforms": {},
            "replayed_length_unit": "mm",
            "replayed_contact_pairs": [],
            "replayed_overlap_volumes_m3": [1.0],
            "replayed_minimum_clearance_m": -1.0,
            "replayed_assembly_generation_id": 723,
            "replayed_shape_owners": {},
            "accepted_assembly_result_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "assembly_interference_replays_use_current_parts_transforms_units_contacts_overlap_clearance_owner_and_result"
    ]


def test_v41_source_step_metadata_mismatch():
    row = _source_v41()
    row["replay_identity"][_STEP].update(
        {
            "color_generation": "step-metadata-723",
            "replayed_product_names": {"part:1": "old"},
            "replayed_colors_rgb": {"part:1": [2.0, -1.0, 0.0]},
            "replayed_layers": {},
            "replayed_assembly_hierarchy": {},
            "replayed_subshape_labels": {},
            "replayed_source_shape_owner": "stale:source",
            "replayed_imported_shape_owner": "stale:import",
            "replayed_step_file_sha256": "d" * 64,
            "accepted_metadata_result_sha256": "e" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_metadata_replays_use_current_names_colors_layers_hierarchy_labels_owners_and_result"
    ]


def test_v41_rejects_self_consistent_nonmonotone_loft_parameterization():
    reference, measured = _public_v41()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_LOFT]["section_parameters"] = [0.0, 0.8, 0.7]
            row[_LOFT]["result_section_parameters"] = [0.0, 0.8, 0.7]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v41_rejects_self_consistent_shell_direction_sign_error():
    reference, measured = _public_v41()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_SHELL]["offset_direction"] = "outward"
            row[_SHELL]["result_offset_direction"] = "outward"
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v41_rejects_self_consistent_nonrigid_interfering_assembly():
    row = _source_v41()
    transform = [
        [2.0, 0.0, 0.0, 0.1],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    value = row["replay_identity"][_ASSEMBLY]
    value["part_transforms"]["slider"] = transform
    value["replayed_part_transforms"]["slider"] = transform
    value["overlap_volumes_m3"] = [1.0e-6]
    value["replayed_overlap_volumes_m3"] = [1.0e-6]
    assert _source_result(row)["status"] == "needs_attention"


def test_v41_rejects_self_consistent_step_color_and_label_reference_error():
    row = _source_v41()
    value = row["replay_identity"][_STEP]
    value["colors_rgb"]["part:1"] = [2.0, 0.0, 0.0]
    value["replayed_colors_rgb"]["part:1"] = [2.0, 0.0, 0.0]
    value["subshape_labels"]["part:9/face:1"] = "ghost"
    value["replayed_subshape_labels"]["part:9/face:1"] = "ghost"
    assert _source_result(row)["status"] == "needs_attention"
