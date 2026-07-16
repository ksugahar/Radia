from __future__ import annotations

import json

from radia_mcp.build123d.server import (
    build123d_jointed_assembly_source_replay_gate,
)
from test_build123d_generalization_v11 import (
    _public_result,
    _public_v11,
    _source_v11,
)


def _public_v12():
    reference, measured = _public_v11()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            placement_digest = ("b" if index == 0 else "c") * 64
            final_digest = row["shape_healing_identity"]["final_brep_sha256"]
            pre_heal_digest = row["shape_healing_identity"][
                "pre_heal_brep_sha256"
            ]
            row["assembly_mass_property_coordinate_identity"] = {
                "assembly_generation": "assembly-generation-44",
                "placement_matrix_generation": "placement-generation-44",
                "centroid_transform_generation": "placement-generation-44",
                "inertia_transform_generation": "placement-generation-44",
                "coordinate_frame_id": "assembly-global-frame",
                "centroid_coordinate_frame_id": "assembly-global-frame",
                "inertia_coordinate_frame_id": "assembly-global-frame",
                "placement_matrix_sha256": placement_digest,
                "centroid_placement_matrix_sha256": placement_digest,
                "inertia_placement_matrix_sha256": placement_digest,
            }
            row["boolean_final_shape_identity"] = {
                "boolean_result_generation": "boolean-generation-44",
                "healing_generation": "healing-generation-44",
                "final_shape_generation": "shape-generation-43",
                "mass_property_shape_generation": "shape-generation-43",
                "validity_shape_generation": "shape-generation-43",
                "topology_shape_generation": "shape-generation-43",
                "pre_heal_brep_sha256": pre_heal_digest,
                "final_brep_sha256": final_digest,
                "mass_property_brep_sha256": final_digest,
                "validity_brep_sha256": final_digest,
                "topology_brep_sha256": final_digest,
            }
    return reference, measured


def _source_v12():
    row = _source_v11()
    identity = row["replay_identity"]
    identity["assembly_mass_property_coordinate_identity"] = {
        "assembly_generation": "assembly-generation-44",
        "placement_matrix_generation": "placement-generation-44",
        "centroid_transform_generation": "placement-generation-44",
        "inertia_transform_generation": "placement-generation-44",
        "coordinate_frame_id": "assembly-global-frame",
        "centroid_coordinate_frame_id": "assembly-global-frame",
        "inertia_coordinate_frame_id": "assembly-global-frame",
        "placement_matrix_sha256": "d" * 64,
        "centroid_placement_matrix_sha256": "d" * 64,
        "inertia_placement_matrix_sha256": "d" * 64,
    }
    identity["boolean_final_shape_report_identity"] = {
        "boolean_result_generation": "boolean-generation-44",
        "healing_generation": "healing-generation-44",
        "final_shape_generation": "shape-generation-44",
        "mass_property_shape_generation": "shape-generation-44",
        "validity_shape_generation": "shape-generation-44",
        "topology_shape_generation": "shape-generation-44",
        "pre_heal_brep_sha256": "e" * 64,
        "final_brep_sha256": "f" * 64,
        "mass_property_brep_sha256": "f" * 64,
        "validity_brep_sha256": "f" * 64,
        "topology_brep_sha256": "f" * 64,
    }
    return row


def test_v12_public_assembly_mass_properties_coordinate_frame_mismatch() -> None:
    reference, measured = _public_v12()
    measured["external_cad"][0]["assembly_mass_property_coordinate_identity"].update(
        {
            "centroid_transform_generation": "placement-generation-43",
            "centroid_coordinate_frame_id": "part-local-frame",
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["assembly_mass_properties_use_final_coordinate_frame"] is False


def test_v12_public_boolean_final_shape_healing_generation_mismatch() -> None:
    reference, measured = _public_v12()
    identity = measured["external_cad"][0]["boolean_final_shape_identity"]
    identity.update(
        {
            "mass_property_shape_generation": "shape-generation-42",
            "mass_property_brep_sha256": identity["pre_heal_brep_sha256"],
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "boolean_validity_topology_and_mass_share_final_healed_shape"
        ]
        is False
    )


def test_v12_source_assembly_mass_properties_coordinate_frame_mismatch() -> None:
    row = _source_v12()
    row["replay_identity"]["assembly_mass_property_coordinate_identity"].update(
        {
            "placement_matrix_generation": "placement-generation-43",
            "placement_matrix_sha256": "1" * 64,
        }
    )
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "assembly_mass_property_report_uses_current_placement_frame"
        ]
        is False
    )


def test_v12_source_boolean_final_shape_healing_generation_mismatch() -> None:
    row = _source_v12()
    identity = row["replay_identity"]["boolean_final_shape_report_identity"]
    identity.update(
        {
            "mass_property_shape_generation": "shape-generation-43",
            "mass_property_brep_sha256": identity["pre_heal_brep_sha256"],
        }
    )
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["final_shape_report_uses_one_post_heal_generation"] is False
