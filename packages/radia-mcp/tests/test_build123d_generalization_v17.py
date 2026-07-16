from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v16 import _public_v16, _source_v16


def _public_v17():
    reference, measured = _public_v16()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("6" if index == 0 else "7") * 64
            row["boolean_tolerance_model_length_unit_generation_identity"] = {
                "model_length_unit_generation": "model-unit-51",
                "boolean_tolerance_unit_generation": "model-unit-51",
                "boolean_result_unit_generation": "model-unit-51",
                "model_length_unit": "mm",
                "tolerance_length_unit": "mm",
                "boolean_result_length_unit": "mm",
                "boolean_tolerance_value": 1.0e-6,
                "boolean_tolerance_si_m": 1.0e-9,
                "boolean_result_tolerance_si_m": 1.0e-9,
                "boolean_tolerance_sha256": digest,
                "boolean_result_tolerance_sha256": digest,
            }
            row["assembly_center_of_mass_part_density_mapping_identity"] = {
                "assembly_configuration_generation": "assembly-config-51",
                "part_density_mapping_generation": "assembly-config-51",
                "center_of_mass_configuration_generation": "assembly-config-51",
                "part_names": ["housing", "shaft"],
                "density_part_names": ["housing", "shaft"],
                "part_densities_kg_m3": [2700.0, 7850.0],
                "center_of_mass_density_values_kg_m3": [2700.0, 7850.0],
                "density_mapping_sha256": digest,
                "center_of_mass_density_mapping_sha256": digest,
            }
    return reference, measured


def _source_v17():
    row = _source_v16()
    identity = row["replay_identity"]
    identity["step_occurrence_name_color_hierarchy_identity"] = {
        "step_import_generation": "step-import-51",
        "occurrence_metadata_import_generation": "step-import-51",
        "assembly_hierarchy_generation": "assembly-hierarchy-51",
        "occurrence_hierarchy_generation": "assembly-hierarchy-51",
        "occurrence_ids": ["root/housing", "root/shaft"],
        "metadata_occurrence_ids": ["root/housing", "root/shaft"],
        "occurrence_names": ["housing", "shaft"],
        "occurrence_colors_rgb": [[0.8, 0.8, 0.8], [0.2, 0.2, 0.2]],
        "occurrence_parent_paths": ["root", "root"],
        "hierarchy_metadata_sha256": "8" * 64,
        "imported_hierarchy_metadata_sha256": "8" * 64,
    }
    identity["brep_edge_tolerance_shape_fix_topology_identity"] = {
        "shape_fix_generation": "shape-fix-51",
        "edge_tolerance_shape_fix_generation": "shape-fix-51",
        "topology_digest_shape_fix_generation": "shape-fix-51",
        "edge_ids": [101, 102, 103],
        "edge_tolerance_edge_ids": [101, 102, 103],
        "edge_tolerances_m": [1.0e-8, 2.0e-8, 1.0e-8],
        "topology_edge_count": 3,
        "topology_sha256": "9" * 64,
        "edge_tolerance_topology_sha256": "9" * 64,
    }
    return row


def test_v17_positive_public_and_source_lineage():
    reference, measured = _public_v17()
    assert _public_result(reference, measured)["status"] == "ok"
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v17()))
    )
    assert source["status"] == "ok"


def test_v17_public_boolean_tolerance_model_length_unit_generation_mismatch():
    reference, measured = _public_v17()
    measured["external_cad"][0][
        "boolean_tolerance_model_length_unit_generation_identity"
    ].update(
        {
            "boolean_tolerance_unit_generation": "model-unit-50",
            "tolerance_length_unit": "m",
            "boolean_tolerance_si_m": 1.0e-6,
            "boolean_result_tolerance_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "boolean_tolerance_uses_current_model_length_unit_generation"
        ]
        is False
    )


def test_v17_public_assembly_center_of_mass_part_density_mapping_generation_mismatch():
    reference, measured = _public_v17()
    measured["external_cad"][0][
        "assembly_center_of_mass_part_density_mapping_identity"
    ].update(
        {
            "part_density_mapping_generation": "assembly-config-50",
            "center_of_mass_density_values_kg_m3": [7850.0, 2700.0],
            "center_of_mass_density_mapping_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "assembly_center_of_mass_uses_current_part_density_mapping"
        ]
        is False
    )


def test_v17_source_step_occurrence_name_color_hierarchy_generation_mismatch():
    row = _source_v17()
    row["replay_identity"][
        "step_occurrence_name_color_hierarchy_identity"
    ].update(
        {
            "occurrence_metadata_import_generation": "step-import-50",
            "occurrence_hierarchy_generation": "assembly-hierarchy-50",
            "occurrence_parent_paths": ["root/old", "root/old"],
            "imported_hierarchy_metadata_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "step_occurrence_metadata_uses_current_assembly_hierarchy"
        ]
        is False
    )


def test_v17_source_brep_edge_tolerance_topology_digest_after_shape_fix_mismatch():
    row = _source_v17()
    row["replay_identity"][
        "brep_edge_tolerance_shape_fix_topology_identity"
    ].update(
        {
            "edge_tolerance_shape_fix_generation": "shape-fix-50",
            "edge_tolerance_edge_ids": [101, 103, 104],
            "edge_tolerance_topology_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["brep_edge_tolerances_follow_final_shape_fix_topology"]
        is False
    )
