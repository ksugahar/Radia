from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v23 import _public_v23, _source_v23


def _public_v24():
    reference, measured = _public_v23()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row[
                "transformed_assembly_com_inertia_axis_density_unit_generation_identity"
            ] = {
                "assembly_generation": "inertia-101",
                "transform_assembly_generation": "inertia-101",
                "density_assembly_generation": "inertia-101",
                "unit_assembly_generation": "inertia-101",
                "mass_property_assembly_generation": "inertia-101",
                "result_assembly_generation": "inertia-101",
                "part_ids": ["base", "rotor", "housing"],
                "result_part_ids": ["base", "rotor", "housing"],
                "local_to_global_transform_sha256": ["3" * 64, "4" * 64, "5" * 64],
                "result_local_to_global_transform_sha256": [
                    "3" * 64,
                    "4" * 64,
                    "5" * 64,
                ],
                "density_kg_m3": [7800.0, 7600.0, 2700.0],
                "result_density_kg_m3": [7800.0, 7600.0, 2700.0],
                "length_unit": "m",
                "result_length_unit": "m",
                "center_of_mass_m": [0.01, -0.02, 0.03],
                "result_center_of_mass_m": [0.01, -0.02, 0.03],
                "inertia_tensor_kg_m2": [
                    [0.5, 0.01, 0.0],
                    [0.01, 0.7, 0.02],
                    [0.0, 0.02, 0.9],
                ],
                "result_inertia_tensor_kg_m2": [
                    [0.5, 0.01, 0.0],
                    [0.01, 0.7, 0.02],
                    [0.0, 0.02, 0.9],
                ],
                "principal_axes_sha256": "6" * 64,
                "result_principal_axes_sha256": "6" * 64,
                "mass_property_sha256": digest,
                "result_mass_property_sha256": digest,
            }
            row[
                "fillet_chamfer_topology_naming_edge_selection_fingerprint_identity"
            ] = {
                "build_generation": "topology-101",
                "selection_build_generation": "topology-101",
                "fillet_build_generation": "topology-101",
                "chamfer_build_generation": "topology-101",
                "naming_build_generation": "topology-101",
                "result_build_generation": "topology-101",
                "operation_order": ["fillet", "chamfer"],
                "result_operation_order": ["fillet", "chamfer"],
                "selected_edge_ids": [11, 14, 18],
                "result_selected_edge_ids": [11, 14, 18],
                "persistent_edge_names": ["rim-a", "rim-b", "key-edge"],
                "result_persistent_edge_names": ["rim-a", "rim-b", "key-edge"],
                "pre_operation_topology_sha256": "7" * 64,
                "result_pre_operation_topology_sha256": "7" * 64,
                "final_topology_sha256": "8" * 64,
                "result_final_topology_sha256": "8" * 64,
                "build_fingerprint_sha256": digest,
                "result_build_fingerprint_sha256": digest,
            }
    return reference, measured


def _source_v24():
    row = _source_v23()
    identity = row["replay_identity"]
    identity[
        "brep_step_roundtrip_tolerance_orientation_volume_generation_identity"
    ] = {
        "roundtrip_generation": "roundtrip-101",
        "source_roundtrip_generation": "roundtrip-101",
        "tolerance_roundtrip_generation": "roundtrip-101",
        "orientation_roundtrip_generation": "roundtrip-101",
        "volume_roundtrip_generation": "roundtrip-101",
        "topology_roundtrip_generation": "roundtrip-101",
        "result_roundtrip_generation": "roundtrip-101",
        "source_format": "STEP-AP214",
        "decoded_source_format": "STEP-AP214",
        "linear_tolerance": 1.0e-7,
        "decoded_linear_tolerance": 1.0e-7,
        "angular_tolerance_deg": 0.1,
        "decoded_angular_tolerance_deg": 0.1,
        "shell_orientation": "outward",
        "decoded_shell_orientation": "outward",
        "volume": 0.0125,
        "decoded_volume": 0.0125,
        "source_shape_sha256": "9" * 64,
        "decoded_source_shape_sha256": "9" * 64,
        "topology_sha256": "a" * 64,
        "decoded_topology_sha256": "a" * 64,
    }
    identity[
        "fresh_subprocess_timeout_exception_cache_output_generation_identity"
    ] = {
        "run_generation": "subprocess-101",
        "process_run_generation": "subprocess-101",
        "interpreter_run_generation": "subprocess-101",
        "cache_run_generation": "subprocess-101",
        "temporary_output_run_generation": "subprocess-101",
        "result_run_generation": "subprocess-101",
        "fresh_interpreter": True,
        "timed_out": False,
        "exception_raised": False,
        "module_cache_preloaded": False,
        "temporary_directory_unique": True,
        "owned_process_count_after": 0,
        "source_script_sha256": "b" * 64,
        "executed_source_script_sha256": "b" * 64,
        "output_shape_sha256": "c" * 64,
        "accepted_output_shape_sha256": "c" * 64,
        "process_log_sha256": "d" * 64,
        "accepted_process_log_sha256": "d" * 64,
    }
    return row


def test_v24_positive_public_and_source_identity():
    reference, measured = _public_v24()
    assert _public_result(reference, measured)["status"] == "ok"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v24()))
    )
    assert result["status"] == "ok"


def test_v24_public_transformed_assembly_com_inertia_principal_axis_density_unit_mismatch():
    reference, measured = _public_v24()
    measured["external_cad"][0][
        "transformed_assembly_com_inertia_axis_density_unit_generation_identity"
    ].update(
        {
            "transform_assembly_generation": "inertia-100",
            "density_assembly_generation": "inertia-99",
            "unit_assembly_generation": "inertia-98",
            "mass_property_assembly_generation": "inertia-97",
            "result_part_ids": ["base", "housing", "rotor"],
            "result_local_to_global_transform_sha256": ["3" * 64, "5" * 64, "4" * 64],
            "result_density_kg_m3": [7800.0, 2700.0, 7600.0],
            "result_length_unit": "mm",
            "result_center_of_mass_m": [10.0, -20.0, 30.0],
            "result_inertia_tensor_kg_m2": [
                [0.9, 0.02, 0.0],
                [0.02, 0.7, 0.01],
                [0.0, 0.01, 0.5],
            ],
            "result_principal_axes_sha256": "e" * 64,
            "result_mass_property_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "transformed_assembly_mass_properties_use_current_transforms_density_units_and_axes"
    ]


def test_v24_public_fillet_chamfer_topology_naming_edge_selection_build_fingerprint_mismatch():
    reference, measured = _public_v24()
    measured["external_cad"][0][
        "fillet_chamfer_topology_naming_edge_selection_fingerprint_identity"
    ].update(
        {
            "selection_build_generation": "topology-100",
            "fillet_build_generation": "topology-99",
            "chamfer_build_generation": "topology-98",
            "naming_build_generation": "topology-97",
            "result_operation_order": ["chamfer", "fillet"],
            "result_selected_edge_ids": [11, 15, 19],
            "result_persistent_edge_names": ["rim-b", "rim-a", "stale-edge"],
            "result_pre_operation_topology_sha256": "0" * 64,
            "result_final_topology_sha256": "1" * 64,
            "result_build_fingerprint_sha256": "2" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fillet_chamfer_topology_uses_current_selection_names_order_and_fingerprint"
    ]


def test_v24_source_brep_step_roundtrip_tolerance_orientation_volume_digest_mismatch():
    row = _source_v24()
    row["replay_identity"][
        "brep_step_roundtrip_tolerance_orientation_volume_generation_identity"
    ].update(
        {
            "source_roundtrip_generation": "roundtrip-100",
            "tolerance_roundtrip_generation": "roundtrip-99",
            "orientation_roundtrip_generation": "roundtrip-98",
            "volume_roundtrip_generation": "roundtrip-97",
            "topology_roundtrip_generation": "roundtrip-96",
            "decoded_source_format": "BREP",
            "decoded_linear_tolerance": 1.0e-3,
            "decoded_angular_tolerance_deg": 5.0,
            "decoded_shell_orientation": "inward",
            "decoded_volume": 0.011,
            "decoded_source_shape_sha256": "3" * 64,
            "decoded_topology_sha256": "4" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "brep_step_roundtrip_uses_current_tolerances_orientation_volume_and_topology"
    ]


def test_v24_source_fresh_subprocess_timeout_exception_cache_output_generation_mismatch():
    row = _source_v24()
    row["replay_identity"][
        "fresh_subprocess_timeout_exception_cache_output_generation_identity"
    ].update(
        {
            "process_run_generation": "subprocess-100",
            "interpreter_run_generation": "subprocess-99",
            "cache_run_generation": "subprocess-98",
            "temporary_output_run_generation": "subprocess-97",
            "fresh_interpreter": False,
            "timed_out": True,
            "exception_raised": True,
            "module_cache_preloaded": True,
            "temporary_directory_unique": False,
            "owned_process_count_after": 1,
            "executed_source_script_sha256": "5" * 64,
            "accepted_output_shape_sha256": "6" * 64,
            "accepted_process_log_sha256": "7" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fresh_subprocess_rejects_timeout_exception_cache_and_stale_output"
    ]
