from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v24 import _public_v24, _source_v24


def _public_v25():
    reference, measured = _public_v24()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["boolean_tolerance_healing_topology_volume_generation_identity"] = {
                "boolean_generation": "boolean-111",
                "operand_boolean_generation": "boolean-111",
                "tolerance_boolean_generation": "boolean-111",
                "healing_boolean_generation": "boolean-111",
                "topology_boolean_generation": "boolean-111",
                "volume_boolean_generation": "boolean-111",
                "result_boolean_generation": "boolean-111",
                "operand_ids": ["base", "tool"],
                "result_operand_ids": ["base", "tool"],
                "linear_tolerance": 1.0e-7,
                "result_linear_tolerance": 1.0e-7,
                "healing_policy": "sew_then_fix_small_edges",
                "result_healing_policy": "sew_then_fix_small_edges",
                "topology_signature": {
                    "solids": 1,
                    "shells": 1,
                    "faces": 14,
                    "edges": 28,
                },
                "result_topology_signature": {
                    "solids": 1,
                    "shells": 1,
                    "faces": 14,
                    "edges": 28,
                },
                "volume_m3": 0.0125,
                "result_volume_m3": 0.0125,
                "operand_shape_sha256": "3" * 64,
                "result_operand_shape_sha256": "3" * 64,
                "boolean_shape_sha256": digest,
                "result_boolean_shape_sha256": digest,
            }
            row["assembly_mate_transform_dof_loop_closure_generation_identity"] = {
                "assembly_generation": "mate-111",
                "mate_assembly_generation": "mate-111",
                "transform_assembly_generation": "mate-111",
                "dof_assembly_generation": "mate-111",
                "closure_assembly_generation": "mate-111",
                "solver_assembly_generation": "mate-111",
                "result_assembly_generation": "mate-111",
                "mate_ids": ["fixed-base", "revolute-rotor", "coincident-cover"],
                "result_mate_ids": [
                    "fixed-base",
                    "revolute-rotor",
                    "coincident-cover",
                ],
                "part_transform_sha256": ["4" * 64, "5" * 64, "6" * 64],
                "result_part_transform_sha256": ["4" * 64, "5" * 64, "6" * 64],
                "remaining_dof": 1,
                "result_remaining_dof": 1,
                "loop_closure_residual_m": 2.0e-12,
                "result_loop_closure_residual_m": 2.0e-12,
                "kinematic_solver_sha256": "7" * 64,
                "result_kinematic_solver_sha256": "7" * 64,
                "assembly_pose_sha256": digest,
                "result_assembly_pose_sha256": digest,
            }
    return reference, measured


def _source_v25():
    row = _source_v24()
    identity = row["replay_identity"]
    identity["step_label_color_unit_hierarchy_shape_roundtrip_identity"] = {
        "roundtrip_generation": "step-meta-111",
        "label_roundtrip_generation": "step-meta-111",
        "color_roundtrip_generation": "step-meta-111",
        "unit_roundtrip_generation": "step-meta-111",
        "hierarchy_roundtrip_generation": "step-meta-111",
        "shape_roundtrip_generation": "step-meta-111",
        "result_roundtrip_generation": "step-meta-111",
        "part_labels": ["base", "rotor", "cover"],
        "decoded_part_labels": ["base", "rotor", "cover"],
        "part_colors_rgb": [
            [0.3, 0.3, 0.3],
            [0.8, 0.2, 0.2],
            [0.2, 0.4, 0.8],
        ],
        "decoded_part_colors_rgb": [
            [0.3, 0.3, 0.3],
            [0.8, 0.2, 0.2],
            [0.2, 0.4, 0.8],
        ],
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "assembly_hierarchy": [
            ["root", "base"],
            ["root", "rotor"],
            ["root", "cover"],
        ],
        "decoded_assembly_hierarchy": [
            ["root", "base"],
            ["root", "rotor"],
            ["root", "cover"],
        ],
        "part_shape_sha256": ["8" * 64, "9" * 64, "a" * 64],
        "decoded_part_shape_sha256": ["8" * 64, "9" * 64, "a" * 64],
        "step_sha256": "b" * 64,
        "decoded_step_sha256": "b" * 64,
    }
    identity[
        "occ_version_tolerance_tessellation_cache_build_fingerprint_identity"
    ] = {
        "build_generation": "occ-build-111",
        "occ_build_generation": "occ-build-111",
        "tolerance_build_generation": "occ-build-111",
        "tessellation_build_generation": "occ-build-111",
        "cache_build_generation": "occ-build-111",
        "fingerprint_build_generation": "occ-build-111",
        "result_build_generation": "occ-build-111",
        "occ_version": "7.8.1",
        "result_occ_version": "7.8.1",
        "linear_tolerance": 1.0e-7,
        "result_linear_tolerance": 1.0e-7,
        "tessellation_linear_deflection": 1.0e-3,
        "result_tessellation_linear_deflection": 1.0e-3,
        "tessellation_angular_deflection_rad": 0.1,
        "result_tessellation_angular_deflection_rad": 0.1,
        "module_cache_fingerprint_sha256": "c" * 64,
        "result_module_cache_fingerprint_sha256": "c" * 64,
        "build_fingerprint_sha256": "d" * 64,
        "result_build_fingerprint_sha256": "d" * 64,
        "tessellation_sha256": "e" * 64,
        "result_tessellation_sha256": "e" * 64,
    }
    return row


def _source_result(row):
    return json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))


def test_v25_positive_public_and_source_identity():
    reference, measured = _public_v25()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v25())["status"] == "ok"


def test_v25_public_boolean_near_tolerance_shape_healing_topology_volume_generation_mismatch():
    reference, measured = _public_v25()
    measured["external_cad"][0][
        "boolean_tolerance_healing_topology_volume_generation_identity"
    ].update(
        {
            "operand_boolean_generation": "boolean-110",
            "tolerance_boolean_generation": "boolean-109",
            "healing_boolean_generation": "boolean-108",
            "topology_boolean_generation": "boolean-107",
            "volume_boolean_generation": "boolean-106",
            "result_operand_ids": ["tool", "base-old"],
            "result_linear_tolerance": 1.0e-3,
            "result_healing_policy": "none",
            "result_topology_signature": {
                "solids": 2,
                "shells": 2,
                "faces": 16,
                "edges": 31,
            },
            "result_volume_m3": 0.0118,
            "result_operand_shape_sha256": "f" * 64,
            "result_boolean_shape_sha256": "0" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_result_uses_current_operands_tolerance_healing_topology_and_volume"
    ]


def test_v25_public_assembly_mate_kinematic_transform_dof_loop_closure_generation_mismatch():
    reference, measured = _public_v25()
    measured["external_cad"][0][
        "assembly_mate_transform_dof_loop_closure_generation_identity"
    ].update(
        {
            "mate_assembly_generation": "mate-110",
            "transform_assembly_generation": "mate-109",
            "dof_assembly_generation": "mate-108",
            "closure_assembly_generation": "mate-107",
            "solver_assembly_generation": "mate-106",
            "result_mate_ids": [
                "fixed-base",
                "revolute-old",
                "coincident-cover",
            ],
            "result_part_transform_sha256": ["4" * 64, "1" * 64, "6" * 64],
            "result_remaining_dof": 3,
            "result_loop_closure_residual_m": 2.0e-2,
            "result_kinematic_solver_sha256": "2" * 64,
            "result_assembly_pose_sha256": "3" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "assembly_mates_use_current_transforms_dof_solver_and_loop_closure"
    ]


def test_v25_source_step_label_color_unit_hierarchy_shape_roundtrip_digest_mismatch():
    row = _source_v25()
    row["replay_identity"][
        "step_label_color_unit_hierarchy_shape_roundtrip_identity"
    ].update(
        {
            "label_roundtrip_generation": "step-meta-110",
            "color_roundtrip_generation": "step-meta-109",
            "unit_roundtrip_generation": "step-meta-108",
            "hierarchy_roundtrip_generation": "step-meta-107",
            "shape_roundtrip_generation": "step-meta-106",
            "decoded_part_labels": ["base", "cover", "rotor-old"],
            "decoded_part_colors_rgb": [
                [0.3, 0.3, 0.3],
                [0.2, 0.4, 0.8],
                [0.8, 0.2, 0.2],
            ],
            "decoded_length_unit": "m",
            "decoded_assembly_hierarchy": [
                ["root", "base"],
                ["sub", "rotor"],
                ["root", "cover"],
            ],
            "decoded_part_shape_sha256": ["8" * 64, "a" * 64, "4" * 64],
            "decoded_step_sha256": "5" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_roundtrip_preserves_labels_colors_units_hierarchy_and_shapes"
    ]


def test_v25_source_occ_version_tolerance_tessellation_cache_build_fingerprint_mismatch():
    row = _source_v25()
    row["replay_identity"][
        "occ_version_tolerance_tessellation_cache_build_fingerprint_identity"
    ].update(
        {
            "occ_build_generation": "occ-build-110",
            "tolerance_build_generation": "occ-build-109",
            "tessellation_build_generation": "occ-build-108",
            "cache_build_generation": "occ-build-107",
            "fingerprint_build_generation": "occ-build-106",
            "result_occ_version": "7.7.2",
            "result_linear_tolerance": 1.0e-4,
            "result_tessellation_linear_deflection": 1.0e-2,
            "result_tessellation_angular_deflection_rad": 0.5,
            "result_module_cache_fingerprint_sha256": "6" * 64,
            "result_build_fingerprint_sha256": "7" * 64,
            "result_tessellation_sha256": "8" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "occ_build_uses_current_version_tolerances_tessellation_cache_and_fingerprint"
    ]
