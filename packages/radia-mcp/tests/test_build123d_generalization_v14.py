from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v13 import _public_v13, _source_v13


def _public_v14():
    reference, measured = _public_v13()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("6" if index == 0 else "7") * 64
            row["center_of_mass_density_length_unit_identity"] = {
                "mass_property_generation": "mass-property-46",
                "geometry_generation": "mass-property-46",
                "density_generation": "mass-property-46",
                "geometry_length_unit": "mm",
                "geometry_length_scale_to_m": 1.0e-3,
                "center_of_mass_length_unit": "mm",
                "center_of_mass_length_scale_to_m": 1.0e-3,
                "density_unit": "kg/m^3",
                "density_scale_to_kg_per_m3": 1.0,
                "volume_length_unit": "mm",
                "volume_length_scale_to_m": 1.0e-3,
                "reported_mass_unit": "kg",
            }
            row["periodic_face_selector_fillet_topology_identity"] = {
                "final_fillet_generation": "fillet-generation-46",
                "final_topology_generation": "topology-generation-46",
                "selector_topology_generation": "topology-generation-46",
                "periodic_pair_topology_generation": "topology-generation-46",
                "source_face_ids": [11, 12],
                "selected_source_face_ids": [11, 12],
                "target_face_ids": [21, 22],
                "selected_target_face_ids": [21, 22],
                "final_shape_sha256": digest,
                "selector_parent_shape_sha256": digest,
            }
    return reference, measured


def _source_v14():
    row = _source_v13()
    identity = row["replay_identity"]
    identity["step_assembly_placement_unit_identity"] = {
        "step_import_generation": "step-import-46",
        "part_geometry_generation": "step-import-46",
        "placement_metadata_generation": "step-import-46",
        "part_geometry_length_unit": "mm",
        "part_geometry_scale_to_m": 1.0e-3,
        "placement_translation_unit": "mm",
        "placement_translation_scale_to_m": 1.0e-3,
        "placement_transform_sha256": "8" * 64,
        "applied_transform_sha256": "8" * 64,
    }
    identity["brep_serialization_tolerance_kernel_identity"] = {
        "cache_generation": "brep-cache-46",
        "active_kernel_generation": "occt-kernel-46",
        "serialization_kernel_generation": "occt-kernel-46",
        "modeling_tolerance_value": 1.0e-6,
        "modeling_tolerance_unit": "m",
        "cache_tolerance_value": 1.0e-6,
        "cache_tolerance_unit": "m",
        "shape_sha256": "9" * 64,
        "cached_shape_sha256": "9" * 64,
    }
    return row


def test_v14_public_center_of_mass_density_length_unit_covariance_mismatch():
    reference, measured = _public_v14()
    measured["external_cad"][0]["center_of_mass_density_length_unit_identity"].update(
        {
            "volume_length_unit": "m",
            "volume_length_scale_to_m": 1.0,
            "density_scale_to_kg_per_m3": 1.0e9,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["center_of_mass_density_and_volume_share_length_unit_covariance"] is False


def test_v14_public_periodic_face_selector_after_fillet_topology_mismatch():
    reference, measured = _public_v14()
    measured["external_cad"][0]["periodic_face_selector_fillet_topology_identity"].update(
        {
            "selector_topology_generation": "topology-generation-45",
            "selected_source_face_ids": [10, 12],
            "selector_parent_shape_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["periodic_face_selectors_follow_final_fillet_topology"] is False


def test_v14_source_step_assembly_placement_unit_transform_mismatch():
    row = _source_v14()
    row["replay_identity"]["step_assembly_placement_unit_identity"].update(
        {
            "placement_metadata_generation": "step-import-45",
            "placement_translation_unit": "m",
            "placement_translation_scale_to_m": 1.0,
            "applied_transform_sha256": "a" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["step_assembly_placement_uses_current_unit_transform_generation"] is False


def test_v14_source_brep_serialization_tolerance_kernel_generation_mismatch():
    row = _source_v14()
    row["replay_identity"]["brep_serialization_tolerance_kernel_identity"].update(
        {
            "serialization_kernel_generation": "occt-kernel-45",
            "cache_tolerance_value": 1.0e-4,
            "cached_shape_sha256": "a" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["brep_cache_uses_current_kernel_tolerance_and_shape"] is False
