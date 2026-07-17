from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v19 import _public_v19, _source_v19


def _mass_properties_identity(digest: str):
    return {
        "assembly_generation": "assembly-54",
        "mass_assembly_generation": "assembly-54",
        "center_of_mass_assembly_generation": "assembly-54",
        "inertia_assembly_generation": "assembly-54",
        "density_mapping_generation": "density-54",
        "mass_density_mapping_generation": "density-54",
        "center_of_mass_density_mapping_generation": "density-54",
        "inertia_density_mapping_generation": "density-54",
        "part_location_generation": "location-54",
        "mass_part_location_generation": "location-54",
        "center_of_mass_part_location_generation": "location-54",
        "inertia_part_location_generation": "location-54",
        "density_unit": "kg/m^3",
        "mass_density_unit": "kg/m^3",
        "center_of_mass_density_unit": "kg/m^3",
        "inertia_density_unit": "kg/m^3",
        "part_names": ["rotor", "housing"],
        "resolved_part_names": ["rotor", "housing"],
        "part_location_sha256": ["1" * 64, "2" * 64],
        "resolved_part_location_sha256": ["1" * 64, "2" * 64],
        "mass_property_table_sha256": digest,
        "resolved_mass_property_table_sha256": digest,
    }


def _sweep_identity(digest: str):
    return {
        "sweep_generation": "sweep-54",
        "solid_sweep_generation": "sweep-54",
        "path_generation": "path-54",
        "frame_path_generation": "path-54",
        "twist_path_generation": "path-54",
        "profile_generation": "profile-54",
        "orientation_profile_generation": "profile-54",
        "solid_profile_generation": "profile-54",
        "path_parameters": [0.0, 0.5, 1.0],
        "frame_path_parameters": [0.0, 0.5, 1.0],
        "twist_degrees": [0.0, 45.0, 90.0],
        "solid_twist_degrees": [0.0, 45.0, 90.0],
        "path_frame_sha256": ["3" * 64, "4" * 64, "5" * 64],
        "solid_path_frame_sha256": ["3" * 64, "4" * 64, "5" * 64],
        "profile_orientation_sha256": "7" * 64,
        "solid_profile_orientation_sha256": "7" * 64,
        "swept_solid_sha256": digest,
        "resolved_swept_solid_sha256": digest,
    }


def _public_v20():
    reference, measured = _public_v19()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("5" if index == 0 else "6") * 64
            row["mass_properties_density_unit_location_generation_identity"] = (
                _mass_properties_identity(digest)
            )
            row[
                "sweep_path_frame_twist_profile_orientation_generation_identity"
            ] = _sweep_identity(digest)
    return reference, measured


def _source_v20():
    row = _source_v19()
    identity = row["replay_identity"]
    identity["brep_serialization_shape_digest_occt_location_generation_identity"] = {
        "serialization_generation": "serialization-54",
        "deserialization_serialization_generation": "serialization-54",
        "shape_generation": "shape-54",
        "serialized_shape_generation": "shape-54",
        "deserialized_shape_generation": "shape-54",
        "kernel_version": "OCCT-7.8",
        "serialized_kernel_version": "OCCT-7.8",
        "deserialized_kernel_version": "OCCT-7.8",
        "location_generation": "location-54",
        "serialized_location_generation": "location-54",
        "deserialized_location_generation": "location-54",
        "shape_sha256": "8" * 64,
        "serialized_shape_sha256": "8" * 64,
        "deserialized_shape_sha256": "8" * 64,
        "top_level_location_sha256": "9" * 64,
        "serialized_top_level_location_sha256": "9" * 64,
        "deserialized_top_level_location_sha256": "9" * 64,
        "brep_payload_sha256": "a" * 64,
        "deserialized_brep_payload_sha256": "a" * 64,
    }
    identity["dxf_wire_plane_orientation_layer_generation_identity"] = {
        "dxf_import_generation": "dxf-import-54",
        "wire_import_generation": "dxf-import-54",
        "plane_import_generation": "dxf-import-54",
        "layer_import_generation": "dxf-import-54",
        "closure_import_generation": "dxf-import-54",
        "plane_generation": "plane-54",
        "wire_plane_generation": "plane-54",
        "extrusion_plane_generation": "plane-54",
        "layer_generation": "layer-54",
        "wire_layer_generation": "layer-54",
        "closure_generation": "closure-54",
        "wire_closure_generation": "closure-54",
        "wire_ids": [101, 102],
        "imported_wire_ids": [101, 102],
        "wire_layers": ["outline", "holes"],
        "imported_wire_layers": ["outline", "holes"],
        "wire_closed": [True, True],
        "imported_wire_closed": [True, True],
        "plane_orientation_sha256": "b" * 64,
        "imported_plane_orientation_sha256": "b" * 64,
        "wire_table_sha256": "c" * 64,
        "imported_wire_table_sha256": "c" * 64,
    }
    return row


def test_v20_positive_public_and_source_cad_identity():
    reference, measured = _public_v20()
    public = _public_result(reference, measured)
    assert public["status"] == "ok"
    assert public["checks"][
        "assembly_mass_properties_share_density_unit_and_part_locations"
    ]
    assert public["checks"][
        "swept_solid_uses_current_path_frames_twist_and_profile_orientation"
    ]
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v20()))
    )
    assert source["status"] == "ok"
    assert source["checks"][
        "brep_deserialization_uses_current_shape_kernel_and_location"
    ]
    assert source["checks"][
        "dxf_wires_use_current_plane_layer_and_closure_generations"
    ]


def test_v20_public_mass_properties_density_unit_location_generation_mismatch():
    reference, measured = _public_v20()
    measured["external_cad"][0][
        "mass_properties_density_unit_location_generation_identity"
    ].update(
        {
            "center_of_mass_density_mapping_generation": "density-53",
            "inertia_part_location_generation": "location-53",
            "center_of_mass_density_unit": "g/cm^3",
            "resolved_part_location_sha256": ["2" * 64, "1" * 64],
            "resolved_mass_property_table_sha256": "d" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "assembly_mass_properties_share_density_unit_and_part_locations"
    ] is False


def test_v20_public_sweep_path_frame_twist_profile_orientation_generation_mismatch():
    reference, measured = _public_v20()
    measured["external_cad"][0][
        "sweep_path_frame_twist_profile_orientation_generation_identity"
    ].update(
        {
            "frame_path_generation": "path-53",
            "twist_path_generation": "path-53",
            "orientation_profile_generation": "profile-53",
            "solid_path_frame_sha256": ["5" * 64, "4" * 64, "3" * 64],
            "solid_twist_degrees": [0.0, 30.0, 60.0],
            "solid_profile_orientation_sha256": "d" * 64,
            "resolved_swept_solid_sha256": "d" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "swept_solid_uses_current_path_frames_twist_and_profile_orientation"
    ] is False


def test_v20_source_brep_serialization_shape_digest_occt_location_generation_mismatch():
    row = _source_v20()
    row["replay_identity"][
        "brep_serialization_shape_digest_occt_location_generation_identity"
    ].update(
        {
            "deserialization_serialization_generation": "serialization-53",
            "deserialized_shape_generation": "shape-53",
            "deserialized_kernel_version": "OCCT-7.7",
            "deserialized_location_generation": "location-53",
            "deserialized_shape_sha256": "d" * 64,
            "deserialized_top_level_location_sha256": "d" * 64,
            "deserialized_brep_payload_sha256": "d" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "brep_deserialization_uses_current_shape_kernel_and_location"
    ] is False


def test_v20_source_dxf_wire_plane_orientation_layer_generation_mismatch():
    row = _source_v20()
    row["replay_identity"][
        "dxf_wire_plane_orientation_layer_generation_identity"
    ].update(
        {
            "plane_import_generation": "dxf-import-53",
            "layer_import_generation": "dxf-import-53",
            "wire_plane_generation": "plane-53",
            "wire_layer_generation": "layer-53",
            "imported_wire_ids": [102, 101],
            "imported_wire_layers": ["holes", "outline"],
            "imported_wire_closed": [True, False],
            "imported_plane_orientation_sha256": "d" * 64,
            "imported_wire_table_sha256": "d" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "dxf_wires_use_current_plane_layer_and_closure_generations"
    ] is False
