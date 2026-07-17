from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v22 import _public_v22, _source_v22


def _public_v23():
    reference, measured = _public_v22()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row[
                "assembly_location_boolean_operand_revision_mass_property_generation_identity"
            ] = {
                "assembly_generation": "assembly-91",
                "location_assembly_generation": "assembly-91",
                "boolean_operand_assembly_generation": "assembly-91",
                "compound_membership_assembly_generation": "assembly-91",
                "density_map_assembly_generation": "assembly-91",
                "mass_property_assembly_generation": "assembly-91",
                "part_ids": ["base", "rotor", "housing"],
                "evaluated_part_ids": ["base", "rotor", "housing"],
                "part_location_sha256": ["3" * 64, "4" * 64, "5" * 64],
                "evaluated_part_location_sha256": ["3" * 64, "4" * 64, "5" * 64],
                "boolean_operand_revisions": ["operand-a-91", "operand-b-91"],
                "evaluated_boolean_operand_revisions": [
                    "operand-a-91",
                    "operand-b-91",
                ],
                "compound_member_ids": [101, 102, 103],
                "evaluated_compound_member_ids": [101, 102, 103],
                "density_map_sha256": "6" * 64,
                "evaluated_density_map_sha256": "6" * 64,
                "mass_property_sha256": digest,
                "evaluated_mass_property_sha256": digest,
            }
            row[
                "loft_spline_tessellation_watertight_volume_generation_identity"
            ] = {
                "shape_generation": "loft-shape-91",
                "spline_shape_generation": "loft-shape-91",
                "tessellation_shape_generation": "loft-shape-91",
                "watertight_shape_generation": "loft-shape-91",
                "volume_shape_generation": "loft-shape-91",
                "spline_sha256": "7" * 64,
                "evaluated_spline_sha256": "7" * 64,
                "chord_tolerance": 1.0e-4,
                "evaluated_chord_tolerance": 1.0e-4,
                "angular_tolerance_deg": 10.0,
                "evaluated_angular_tolerance_deg": 10.0,
                "length_unit": "m",
                "evaluated_length_unit": "m",
                "watertight": True,
                "evaluated_watertight": True,
                "tessellated_shell_sha256": "8" * 64,
                "evaluated_tessellated_shell_sha256": "8" * 64,
                "volume": 0.0125,
                "evaluated_volume": 0.0125,
                "volume_result_sha256": digest,
                "evaluated_volume_result_sha256": digest,
            }
    return reference, measured


def _source_v23():
    row = _source_v22()
    identity = row["replay_identity"]
    identity["step_import_label_color_unit_topology_generation_identity"] = {
        "import_generation": "step-import-91",
        "source_content_import_generation": "step-import-91",
        "label_import_generation": "step-import-91",
        "color_import_generation": "step-import-91",
        "unit_import_generation": "step-import-91",
        "topology_import_generation": "step-import-91",
        "source_content_sha256": "9" * 64,
        "imported_source_content_sha256": "9" * 64,
        "shape_labels": ["assembly", "rotor", "housing"],
        "decoded_shape_labels": ["assembly", "rotor", "housing"],
        "shape_colors_rgb": [[180, 180, 180], [220, 40, 40], [80, 100, 140]],
        "decoded_shape_colors_rgb": [
            [180, 180, 180],
            [220, 40, 40],
            [80, 100, 140],
        ],
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "unit_scale_to_m": 0.001,
        "decoded_unit_scale_to_m": 0.001,
        "brep_topology_sha256": "a" * 64,
        "decoded_brep_topology_sha256": "a" * 64,
    }
    identity[
        "mesh_export_facet_normal_tolerance_shape_digest_generation_identity"
    ] = {
        "mesh_export_generation": "mesh-export-91",
        "shape_mesh_export_generation": "mesh-export-91",
        "facet_mesh_export_generation": "mesh-export-91",
        "normal_mesh_export_generation": "mesh-export-91",
        "tolerance_mesh_export_generation": "mesh-export-91",
        "source_shape_sha256": "b" * 64,
        "exported_source_shape_sha256": "b" * 64,
        "facet_ids": [201, 202, 203],
        "exported_facet_ids": [201, 202, 203],
        "facet_normals": [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
        "exported_facet_normals": [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
        "chord_tolerance": 1.0e-4,
        "exported_chord_tolerance": 1.0e-4,
        "angular_tolerance_deg": 10.0,
        "exported_angular_tolerance_deg": 10.0,
        "facet_topology_sha256": "c" * 64,
        "exported_facet_topology_sha256": "c" * 64,
    }
    return row


def test_v23_positive_public_and_source_identity():
    reference, measured = _public_v23()
    assert _public_result(reference, measured)["status"] == "ok"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v23()))
    )
    assert result["status"] == "ok"


def test_v23_public_assembly_location_boolean_operand_revision_mass_property_mismatch():
    reference, measured = _public_v23()
    measured["external_cad"][0][
        "assembly_location_boolean_operand_revision_mass_property_generation_identity"
    ].update(
        {
            "location_assembly_generation": "assembly-90",
            "boolean_operand_assembly_generation": "assembly-89",
            "compound_membership_assembly_generation": "assembly-88",
            "density_map_assembly_generation": "assembly-87",
            "mass_property_assembly_generation": "assembly-86",
            "evaluated_part_ids": ["base", "housing", "rotor"],
            "evaluated_part_location_sha256": ["3" * 64, "5" * 64, "4" * 64],
            "evaluated_boolean_operand_revisions": [
                "operand-a-90",
                "operand-b-90",
            ],
            "evaluated_compound_member_ids": [101, 103, 104],
            "evaluated_density_map_sha256": "d" * 64,
            "evaluated_mass_property_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "assembly_mass_properties_use_current_locations_operands_members_and_density"
    ]


def test_v23_public_loft_spline_tessellation_tolerance_watertight_volume_generation_mismatch():
    reference, measured = _public_v23()
    measured["external_cad"][0][
        "loft_spline_tessellation_watertight_volume_generation_identity"
    ].update(
        {
            "spline_shape_generation": "loft-shape-90",
            "tessellation_shape_generation": "loft-shape-89",
            "watertight_shape_generation": "loft-shape-88",
            "volume_shape_generation": "loft-shape-87",
            "evaluated_spline_sha256": "f" * 64,
            "evaluated_chord_tolerance": 0.1,
            "evaluated_angular_tolerance_deg": 25.0,
            "evaluated_length_unit": "mm",
            "evaluated_watertight": False,
            "evaluated_tessellated_shell_sha256": "0" * 64,
            "evaluated_volume": 0.011,
            "evaluated_volume_result_sha256": "1" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "loft_volume_uses_current_spline_tolerances_and_watertight_shell"
    ]


def test_v23_source_step_import_label_color_unit_topology_hash_generation_mismatch():
    row = _source_v23()
    row["replay_identity"][
        "step_import_label_color_unit_topology_generation_identity"
    ].update(
        {
            "source_content_import_generation": "step-import-90",
            "label_import_generation": "step-import-89",
            "color_import_generation": "step-import-88",
            "unit_import_generation": "step-import-87",
            "topology_import_generation": "step-import-86",
            "imported_source_content_sha256": "2" * 64,
            "decoded_shape_labels": ["assembly", "housing", "old_rotor"],
            "decoded_shape_colors_rgb": [
                [180, 180, 180],
                [80, 100, 140],
                [220, 40, 40],
            ],
            "decoded_length_unit": "m",
            "decoded_unit_scale_to_m": 1.0,
            "decoded_brep_topology_sha256": "3" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_import_uses_current_source_labels_colors_units_and_topology"
    ]


def test_v23_source_mesh_export_facet_normal_tolerance_shape_digest_generation_mismatch():
    row = _source_v23()
    row["replay_identity"][
        "mesh_export_facet_normal_tolerance_shape_digest_generation_identity"
    ].update(
        {
            "shape_mesh_export_generation": "mesh-export-90",
            "facet_mesh_export_generation": "mesh-export-89",
            "normal_mesh_export_generation": "mesh-export-88",
            "tolerance_mesh_export_generation": "mesh-export-87",
            "exported_source_shape_sha256": "4" * 64,
            "exported_facet_ids": [203, 202, 204],
            "exported_facet_normals": [
                [0.0, 0.0, -1.0],
                [0.0, -1.0, 0.0],
                [-1.0, 0.0, 0.0],
            ],
            "exported_chord_tolerance": 0.1,
            "exported_angular_tolerance_deg": 30.0,
            "exported_facet_topology_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mesh_export_uses_current_shape_facets_normals_and_tolerances"
    ]
