"""Payload builders for test_build123d_generalization.py (not collected by pytest)."""

from __future__ import annotations

import copy
import json
import math
from copy import deepcopy

from radia_mcp.build123d.assembly_exchange_identity_v53 import (
    HISTORY,
    MASS,
    MATE as MATE_V53,
    STEP as STEP_V53,
)
from radia_mcp.build123d.assembly_replay_identity_v49 import (
    ASSEMBLY as ASSEMBLY_V49,
    BOOLEAN as BOOLEAN_V49,
    SKETCH as SKETCH_V49,
    STEP as STEP_V49,
)
from radia_mcp.build123d.assembly_tessellation_identity_v54 import (
    ASSEMBLY as ASSEMBLY_V54,
    LOFT as LOFT_V54,
    STEP as STEP_V54,
    TESSELLATION,
)
from radia_mcp.build123d.cross_artifact_cad_lineage_v47 import (
    ASSEMBLY as ASSEMBLY_V47,
    COMPOUND,
    EXTERNAL,
    ROUNDTRIP,
)
from radia_mcp.build123d.feature_replay_identity_v50 import (
    FEATURE,
    HEALING,
    MATE as MATE_V50,
    STL,
)
from radia_mcp.build123d.gear_pipe_exchange_identity_v56 import (
    GEAR,
    MESH,
    PIPE,
    STEP as STEP_V56,
)
from radia_mcp.build123d.selector_exchange_identity_v52 import (
    BREP as BREP_V52,
    GLTF,
    SELECTOR as SELECTOR_V52,
    WORKPLANE,
)
from radia_mcp.build123d.semantic_cad_identity_v48 import (
    IMPORT,
    LOFT as LOFT_V48,
    SELECTOR as SELECTOR_V48,
    SUPPRESSION,
)
from radia_mcp.build123d.server import (
    build123d_jointed_assembly_source_replay_gate,
    build123d_mass_property_crosscheck,
)
from radia_mcp.build123d.thread_sheet_identity_v55 import (
    BREP as BREP_V55,
    PMI,
    SHEET,
    THREAD,
)


def _public_v7():
    reference = [_box("frame", (2.0, 3.0, 4.0)), _box("insert", (1.0, 2.0, 5.0))]
    return reference, {"external_cad": copy.deepcopy(reference)}


def _source_v7():
    return {
        "source_kind": "upstream_source_native_example_with_display_stub_only",
        "source_sha256": "a" * 64,
        "source_url": "https://example.invalid/project/blob/v0.10.0/examples/model.py",
        "source_preserved": True,
        "display_stubbed_only": True,
        "components": [
            {"name": "frame", "joint_names": ["frame_joint"]},
            {"name": "insert", "joint_names": ["insert_joint"]},
        ],
        "joint_connections": [
            {"from": "frame_joint", "to": "insert_joint", "kind": "rigid"}
        ],
        "external_execution": {
            "mode": "python_api_headless_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "diagnosis_gate_status": "ok",
        "diagnosis": "component_solid_closure_loss",
        "solver_ready": False,
        "timing_breakdown_s": {
            "source_replay": 0.2,
            "neutral_cad_export": 0.1,
            "external_replay": 0.3,
            "identity_validation": 0.05,
        },
        "replay_identity": {
            "source_commit": "b" * 40,
            "replayed_source_commit": "b" * 40,
            "cad_artifacts": [
                {
                    "name": "assembly.step",
                    "sha256": "c" * 64,
                    "fresh": True,
                    "source_commit": "b" * 40,
                }
            ],
            "external_kernel": {
                "name": "OCCT",
                "claimed_version": "7.8.1",
                "replay_versions": ["7.8.1", "7.8.1"],
            },
        },
    }


def _box(name: str, size: tuple[float, float, float]) -> dict:
    x, y, z = size
    return {
        "name": name,
        "type": "Solid",
        "is_valid": True,
        "volume": x * y * z,
        "area": 2.0 * (x * y + y * z + x * z),
        "faces": 6,
        "edges": 12,
        "vertices": 8,
        "solids": 1,
        "bounding_box": {
            "min": [0.0, 0.0, 0.0],
            "max": [x, y, z],
            "center": [x / 2.0, y / 2.0, z / 2.0],
            "size": [x, y, z],
            "diagonal": (x * x + y * y + z * z) ** 0.5,
        },
    }


def _public() -> tuple[list[dict], dict[str, list[dict]]]:
    reference = [
        _box("frame", (2.0, 3.0, 4.0)),
        _box("insert", (1.0, 2.0, 5.0)),
    ]
    measured = copy.deepcopy(reference)
    child_revisions = {"frame": "child-frame-5", "insert": "child-insert-3"}
    for rows in (reference, measured):
        for row in rows:
            name = row["name"]
            revision = f"brep-{name}-42"
            digest = ("d" if name == "frame" else "e") * 64
            row["brep_identity"] = {"revision": revision, "sha256": digest}
            row["mass_property_identity"] = {
                "brep_revision": revision,
                "brep_sha256": digest,
            }
            row["assembly_identity"] = {
                "generation": "assembly-generation-42",
                "child_revisions": copy.deepcopy(child_revisions),
            }
            row["mass_property_frame_identity"] = {
                "frame_id": "assembly-global-frame-42",
                "transform_generation": "assembly-transform-42",
            }
            row["topology_identity"] = {
                "brep_revision": revision,
                "brep_sha256": digest,
                "face_adjacency_sha256": ("1" if name == "frame" else "2") * 64,
            }
            row["compound_volume_identity"] = {
                "topology_kind": "physical_union_solid",
                "reported_volume_basis": "physical_union",
                "overlap_volume": 0.0,
                "topology_generation": "compound-generation-42",
                "volume_generation": "compound-generation-42",
            }
            row["placement_transform_identity"] = {
                "center_of_mass_frame": "assembly-global-frame-42",
                "center_of_mass_transform_generation": "assembly-transform-42",
                "final_placement_transform_generation": "assembly-transform-42",
            }
    return reference, {"external_cad": measured}


def _public_result(reference: list[dict], measured: dict[str, list[dict]]) -> dict:
    return json.loads(
        build123d_mass_property_crosscheck(
            json.dumps(reference),
            json.dumps(measured),
            rtol=1.0e-10,
            bbox_atol=1.0e-10,
        )
    )


def _source() -> dict:
    return {
        "source_kind": "upstream_source_native_example_with_display_stub_only",
        "source_sha256": "a" * 64,
        "source_url": "https://example.invalid/project/blob/v0.10.0/examples/model.py",
        "source_preserved": True,
        "display_stubbed_only": True,
        "components": [
            {"name": "frame", "joint_names": ["frame_joint"]},
            {"name": "insert", "joint_names": ["insert_joint"]},
        ],
        "joint_connections": [
            {"from": "frame_joint", "to": "insert_joint", "kind": "rigid"}
        ],
        "external_execution": {
            "mode": "python_api_headless_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "diagnosis_gate_status": "ok",
        "diagnosis": "component_solid_closure_loss",
        "solver_ready": False,
        "timing_breakdown_s": {
            "source_replay": 0.2,
            "neutral_cad_export": 0.1,
            "external_replay": 0.3,
            "identity_validation": 0.05,
        },
        "replay_identity": {
            "source_commit": "b" * 40,
            "replayed_source_commit": "b" * 40,
            "source_replay_started_utc": "2026-07-16T02:00:00Z",
            "cad_artifacts": [
                {
                    "name": "assembly.step",
                    "sha256": "c" * 64,
                    "fresh": True,
                    "source_commit": "b" * 40,
                    "export_completed_utc": "2026-07-16T02:00:01Z",
                }
            ],
            "external_kernel": {
                "name": "OCCT",
                "claimed_version": "7.8.1",
                "replay_versions": ["7.8.1", "7.8.1"],
                "claimed_session_generation": "occt-session-42",
                "replay_sessions": [
                    {
                        "session_generation": "occt-session-42",
                        "process_start_utc": "2026-07-16T01:59:00Z",
                    },
                    {
                        "session_generation": "occt-session-42",
                        "process_start_utc": "2026-07-16T01:59:00Z",
                    },
                ],
            },
            "topology_replay_identity": {
                "source_topology_sha256": "3" * 64,
                "imports": [
                    {"mode": "noheal", "topology_sha256": "3" * 64},
                    {"mode": "heal", "topology_sha256": "3" * 64},
                ],
            },
            "unit_conversion_identity": {
                "source_geometry_unit": "mm",
                "target_geometry_unit": "m",
                "length_scale_to_target": 0.001,
                "external_measurement_stage": "after_unit_conversion",
                "external_volume_unit": "m^3",
                "declared_volume_scale_to_target": 1.0,
            },
            "boolean_clean_identity": {
                "boolean_result_sha256": "4" * 64,
                "shape_clean_input_sha256": "4" * 64,
                "cleaned_topology_sha256": "5" * 64,
                "export_topology_sha256": "5" * 64,
                "shape_generation": "shape-generation-42",
                "export_shape_generation": "shape-generation-42",
            },
            "tessellation_identity": {
                "shape_generation": "shape-generation-42",
                "tolerance_shape_generation": "shape-generation-42",
                "linear_deflection": 0.01,
                "angular_tolerance_rad": 0.1,
                "tessellation_generation": "tessellation-generation-42",
            },
        },
    }


def _public_v11():
    reference, measured = _public()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            final_digest = ("6" if index == 0 else "7") * 64
            row["shape_healing_identity"] = {
                "pre_heal_brep_sha256": ("4" if index == 0 else "5") * 64,
                "healed_brep_sha256": final_digest,
                "final_brep_sha256": final_digest,
                "mass_property_brep_sha256": final_digest,
                "final_shape_generation": "shape-generation-43",
                "mass_property_shape_generation": "shape-generation-43",
            }
            frame = row["mass_property_frame_identity"]["frame_id"]
            transform = row["placement_transform_identity"][
                "final_placement_transform_generation"
            ]
            row["inertia_tensor_identity"] = {
                "tensor_frame_id": frame,
                "center_of_mass_frame_id": frame,
                "tensor_transform_generation": transform,
                "final_placement_transform_generation": transform,
                "mirror_transform_applied": True,
                "mirror_transform_determinant": -1.0,
                "tensor_basis_handedness": "right_handed",
            }
    return reference, measured


def _source_v11():
    row = _source()
    identity = row["replay_identity"]
    identity["step_export_tolerance_identity"] = {
        "kernel_session_generation": "occt-session-42",
        "tolerance_kernel_session_generation": "occt-session-42",
        "shape_generation": "shape-generation-42",
        "tolerance_shape_generation": "shape-generation-42",
        "sewing_tolerance": 1.0e-7,
        "brep_tolerance": 1.0e-7,
        "export_artifact_sha256": "8" * 64,
    }
    identity["assembly_replacement_identity"] = {
        "assembly_generation": "assembly-generation-43",
        "components": [
            {
                "slot_id": "insert-slot",
                "replacement_generation": "replacement-generation-43",
                "removed_instance_uuid": "instance-uuid-42",
                "current_instance_uuid": "instance-uuid-43",
                "removed_shape_sha256": "9" * 64,
                "current_shape_sha256": "a" * 64,
                "placement_shape_sha256": "a" * 64,
                "placement_assembly_generation": "assembly-generation-43",
            }
        ],
    }
    return row


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


def _public_v13():
    reference, measured = _public_v12()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["tessellation_tolerance_unit_identity"] = {
                "tessellation_generation": "tessellation-generation-45", "surface_area_generation": "tessellation-generation-45",
                "linear_deflection_value": 0.05, "linear_deflection_unit": "mm", "linear_deflection_scale_to_m": 1.0e-3,
                "area_evaluation_deflection_unit": "mm", "area_evaluation_deflection_scale_to_m": 1.0e-3,
            }
            row["compound_label_topology_identity"] = {
                "boolean_generation": "boolean-generation-45", "label_table_boolean_generation": "boolean-generation-45", "selector_boolean_generation": "boolean-generation-45",
                "label": "mounting_face", "topology_index": 7, "label_topology_index": 7,
                "final_shape_sha256": digest, "selected_subshape_parent_sha256": digest,
            }
    return reference, measured


def _source_v13():
    row = _source_v12()
    identity = row["replay_identity"]
    identity["step_geometry_unit_scale_identity"] = {
        "step_import_generation": "step-import-45", "geometry_coordinate_generation": "step-import-45", "metadata_generation": "step-import-45",
        "geometry_length_unit": "m", "metadata_length_unit": "m", "geometry_scale_to_m": 1.0, "metadata_scale_to_m": 1.0,
    }
    identity["selector_cache_shape_identity"] = {
        "active_shape_generation": "shape-generation-45", "selector_cache_shape_generation": "shape-generation-45", "selected_face_shape_generation": "shape-generation-45",
        "selector_query_sha256": "3" * 64, "cached_selector_query_sha256": "3" * 64, "selected_face_ids": [7, 8], "live_face_ids": [7, 8],
    }
    return row


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


def _public_v15():
    reference, measured = _public_v14()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("b" if index == 0 else "c") * 64
            row["boolean_tolerance_length_unit_identity"] = {
                "boolean_generation": "boolean-47",
                "result_geometry_generation": "boolean-47",
                "model_length_unit": "mm",
                "tolerance_value": 1.0e-6,
                "tolerance_unit": "mm",
                "tolerance_scale_to_m": 1.0e-3,
                "kernel_tolerance_m": 1.0e-9,
                "input_shape_sha256": digest,
                "boolean_input_shape_sha256": digest,
            }
            row["nested_assembly_placement_order_identity"] = {
                "assembly_generation": "assembly-47",
                "parent_placement_generation": "assembly-47",
                "child_placement_generation": "assembly-47",
                "world_placement_generation": "assembly-47",
                "multiplication_order": "parent_then_child",
                "applied_multiplication_order": "parent_then_child",
                "placement_chain_sha256": digest,
                "world_transform_sha256": digest,
            }
    return reference, measured


def _source_v15():
    row = _source_v14()
    identity = row["replay_identity"]
    identity["step_color_label_topology_identity"] = {
        "step_import_generation": "step-import-47",
        "topology_generation": "step-topology-47",
        "attribute_map_topology_generation": "step-topology-47",
        "face_ids": [31, 32, 33],
        "attribute_face_ids": [31, 32, 33],
        "labels": ["housing", "shaft", "terminal"],
        "colors_rgb": [[0.8, 0.8, 0.8], [0.2, 0.2, 0.2], [0.8, 0.1, 0.1]],
        "attribute_map_sha256": "d" * 64,
        "resolved_attribute_map_sha256": "d" * 64,
    }
    identity["brep_surface_parameter_orientation_identity"] = {
        "serialization_generation": "brep-serialization-47",
        "surface_parameter_generation": "brep-serialization-47",
        "surface_ids": [41, 42],
        "exported_surface_ids": [41, 42],
        "parameter_orientation": "u_cross_v_outward",
        "exported_parameter_orientation": "u_cross_v_outward",
        "parameter_range_sha256": "e" * 64,
        "exported_parameter_range_sha256": "e" * 64,
    }
    return row


def _public_v16():
    reference, measured = _public_v15()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["mass_inertia_reference_frame_placement_identity"] = {
                "shape_generation": "shape-50",
                "mass_property_shape_generation": "shape-50",
                "placement_generation": "placement-50",
                "mass_property_placement_generation": "placement-50",
                "mass_reference_frame": "world",
                "inertia_reference_frame": "world",
                "center_of_mass_reference_frame": "world",
                "placed_shape_sha256": digest,
                "mass_property_shape_sha256": digest,
            }
            row["loft_wire_correspondence_seam_identity"] = {
                "loft_generation": "loft-50",
                "section_wire_loft_generation": "loft-50",
                "seam_normalization_generation": "seam-50",
                "wire_correspondence_seam_generation": "seam-50",
                "section_wire_ids": [11, 12, 13],
                "loft_section_wire_ids": [11, 12, 13],
                "wire_correspondence_sha256": digest,
                "loft_wire_correspondence_sha256": digest,
            }
    return reference, measured


def _source_v16():
    row = _source_v15()
    identity = row["replay_identity"]
    identity["step_assembly_child_parent_unit_identity"] = {
        "step_import_generation": "step-import-50",
        "child_placement_import_generation": "step-import-50",
        "parent_placement_import_generation": "step-import-50",
        "assembly_length_unit": "mm",
        "child_placement_length_unit": "mm",
        "parent_placement_length_unit": "mm",
        "assembly_scale_to_m": 1.0e-3,
        "child_placement_scale_to_m": 1.0e-3,
        "parent_placement_scale_to_m": 1.0e-3,
        "assembly_placement_sha256": "3" * 64,
        "resolved_assembly_placement_sha256": "3" * 64,
    }
    identity["selector_normal_world_frame_identity"] = {
        "shape_generation": "selector-shape-50",
        "selector_shape_generation": "selector-shape-50",
        "placement_generation": "selector-placement-50",
        "selector_placement_generation": "selector-placement-50",
        "normal_predicate_frame": "world",
        "evaluated_normal_frame": "world",
        "selected_face_ids": [21, 22],
        "resolved_face_ids": [21, 22],
        "normal_table_sha256": "4" * 64,
        "evaluated_normal_table_sha256": "4" * 64,
    }
    return row


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


def _public_v18():
    reference, measured = _public_v17()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("6" if index == 0 else "7") * 64
            row["nested_assembly_location_transform_composition_identity"] = {
                "assembly_generation": "assembly-52",
                "location_tree_assembly_generation": "assembly-52",
                "composition_assembly_generation": "assembly-52",
                "transform_order_generation": "transform-order-52",
                "location_transform_order_generation": "transform-order-52",
                "child_paths": ["root/frame", "root/frame/shaft"],
                "location_child_paths": ["root/frame", "root/frame/shaft"],
                "local_transform_sha256": ["1" * 64, "2" * 64],
                "composed_transform_sha256": ["3" * 64, "4" * 64],
                "resolved_composed_transform_sha256": ["3" * 64, "4" * 64],
                "composition_order": "parent_then_child",
                "resolved_composition_order": "parent_then_child",
                "location_tree_sha256": digest,
                "resolved_location_tree_sha256": digest,
            }
            row["boolean_retained_face_name_history_refine_identity"] = {
                "boolean_generation": "boolean-52",
                "refine_generation": "refine-52",
                "retained_name_boolean_generation": "boolean-52",
                "retained_name_refine_generation": "refine-52",
                "topology_history_refine_generation": "refine-52",
                "retained_face_names": ["inlet", "outlet"],
                "resolved_face_names": ["inlet", "outlet"],
                "retained_face_ids": [41, 42],
                "resolved_face_ids": [41, 42],
                "topology_history_sha256": digest,
                "resolved_topology_history_sha256": digest,
            }
    return reference, measured


def _source_v18():
    row = _source_v17()
    identity = row["replay_identity"]
    identity["step_occurrence_color_material_inheritance_identity"] = {
        "step_import_generation": "step-import-52",
        "assembly_generation": "assembly-52",
        "occurrence_metadata_import_generation": "step-import-52",
        "occurrence_metadata_assembly_generation": "assembly-52",
        "color_inheritance_assembly_generation": "assembly-52",
        "material_inheritance_assembly_generation": "assembly-52",
        "occurrence_ids": ["root/frame", "root/frame/shaft"],
        "metadata_occurrence_ids": ["root/frame", "root/frame/shaft"],
        "parent_occurrence_ids": ["root", "root/frame"],
        "metadata_parent_occurrence_ids": ["root", "root/frame"],
        "inherited_colors_rgb": [[0.7, 0.7, 0.7], [0.2, 0.2, 0.2]],
        "imported_colors_rgb": [[0.7, 0.7, 0.7], [0.2, 0.2, 0.2]],
        "inherited_material_names": ["aluminum", "steel"],
        "imported_material_names": ["aluminum", "steel"],
        "occurrence_metadata_sha256": "8" * 64,
        "imported_occurrence_metadata_sha256": "8" * 64,
    }
    identity["stl_tolerance_model_length_unit_generation_identity"] = {
        "model_length_unit_generation": "model-unit-52",
        "tessellation_model_unit_generation": "model-unit-52",
        "tolerance_conversion_model_unit_generation": "model-unit-52",
        "model_length_unit": "mm",
        "chordal_tolerance_length_unit": "mm",
        "angular_tolerance_unit": "deg",
        "chordal_tolerance_value": 0.01,
        "chordal_tolerance_si_m": 1.0e-5,
        "tessellator_chordal_tolerance_si_m": 1.0e-5,
        "angular_tolerance_value": 5.0,
        "angular_tolerance_rad": 0.08726646259971647,
        "tessellator_angular_tolerance_rad": 0.08726646259971647,
        "tolerance_contract_sha256": "9" * 64,
        "tessellator_tolerance_contract_sha256": "9" * 64,
    }
    return row


def _public_v19():
    reference, measured = _public_v18()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("a" if index == 0 else "b") * 64
            row["boolean_history_subshape_label_fillet_order_identity"] = {
                "boolean_generation": "boolean-53",
                "fillet_generation": "fillet-53",
                "history_boolean_generation": "boolean-53",
                "label_fillet_generation": "fillet-53",
                "edge_order_fillet_generation": "fillet-53",
                "subshape_labels": ["inlet", "outlet"],
                "resolved_subshape_labels": ["inlet", "outlet"],
                "history_face_ids": [51, 52],
                "resolved_face_ids": [51, 52],
                "fillet_edge_ids": [71, 72, 73],
                "resolved_fillet_edge_ids": [71, 72, 73],
                "subshape_history_sha256": digest,
                "resolved_subshape_history_sha256": digest,
            }
            row["assembly_mate_frame_unit_location_generation_identity"] = {
                "assembly_generation": "assembly-53",
                "mate_frame_assembly_generation": "assembly-53",
                "parent_location_assembly_generation": "assembly-53",
                "unit_assembly_generation": "assembly-53",
                "length_unit": "mm",
                "mate_frame_length_unit": "mm",
                "parent_location_length_unit": "mm",
                "mate_names": ["shaft_axis", "housing_axis"],
                "resolved_mate_names": ["shaft_axis", "housing_axis"],
                "local_frame_sha256": ["1" * 64, "2" * 64],
                "resolved_local_frame_sha256": ["1" * 64, "2" * 64],
                "parent_location_sha256": ["3" * 64, "4" * 64],
                "resolved_parent_location_sha256": ["3" * 64, "4" * 64],
                "mate_resolution_sha256": digest,
                "resolved_mate_resolution_sha256": digest,
            }
    return reference, measured


def _source_v19():
    row = _source_v18()
    identity = row["replay_identity"]
    identity["step_import_tolerance_unit_healing_generation_identity"] = {
        "step_import_generation": "step-import-53",
        "healing_generation": "healing-53",
        "tolerance_import_generation": "step-import-53",
        "tolerance_healing_generation": "healing-53",
        "healed_edge_map_import_generation": "step-import-53",
        "healed_edge_map_healing_generation": "healing-53",
        "source_length_unit": "mm",
        "tolerance_length_unit": "mm",
        "tolerance_value": 1.0e-5,
        "tolerance_si_m": 1.0e-8,
        "healed_edge_ids": [101, 102, 103],
        "imported_healed_edge_ids": [101, 102, 103],
        "healed_edge_map_sha256": "c" * 64,
        "imported_healed_edge_map_sha256": "c" * 64,
    }
    identity["tessellation_vertex_index_normal_transform_generation_identity"] = {
        "shape_generation": "shape-53",
        "location_transform_generation": "location-53",
        "vertex_shape_generation": "shape-53",
        "index_shape_generation": "shape-53",
        "normal_shape_generation": "shape-53",
        "vertex_location_transform_generation": "location-53",
        "index_location_transform_generation": "location-53",
        "normal_location_transform_generation": "location-53",
        "vertex_count": 4,
        "triangle_indices": [[0, 1, 2], [0, 2, 3]],
        "normal_count": 4,
        "transformed_triangle_indices": [[0, 1, 2], [0, 2, 3]],
        "tessellation_sha256": "d" * 64,
        "rendered_tessellation_sha256": "d" * 64,
    }
    return row


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


def _public_v21():
    reference, measured = _public_v20()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["boolean_result_solid_orientation_location_label_generation_identity"] = {
                "boolean_generation": "boolean-71",
                "result_boolean_generation": "boolean-71",
                "orientation_boolean_generation": "boolean-71",
                "location_boolean_generation": "boolean-71",
                "label_boolean_generation": "boolean-71",
                "operand_generations": ["operand-a-71", "operand-b-71"],
                "result_operand_generations": ["operand-a-71", "operand-b-71"],
                "solid_orientation": "forward",
                "resolved_solid_orientation": "forward",
                "semantic_labels": ["body", "interface"],
                "resolved_semantic_labels": ["body", "interface"],
                "result_location_sha256": "3" * 64,
                "resolved_result_location_sha256": "3" * 64,
                "boolean_result_sha256": digest,
                "resolved_boolean_result_sha256": digest,
            }
            row["tessellation_chord_angle_unit_location_generation_identity"] = {
                "shape_generation": "shape-71",
                "tessellation_shape_generation": "shape-71",
                "tessellation_generation": "tessellation-71",
                "metric_tessellation_generation": "tessellation-71",
                "location_tessellation_generation": "tessellation-71",
                "chord_tolerance": 1.0e-4,
                "evaluated_chord_tolerance": 1.0e-4,
                "angular_tolerance_deg": 12.0,
                "evaluated_angular_tolerance_deg": 12.0,
                "length_unit": "m",
                "evaluated_length_unit": "m",
                "object_location_sha256": "4" * 64,
                "evaluated_object_location_sha256": "4" * 64,
                "tessellation_sha256": digest,
                "evaluated_tessellation_sha256": digest,
            }
    return reference, measured


def _source_v21():
    row = _source_v20()
    identity = row["replay_identity"]
    identity["step_assembly_product_id_color_location_generation_identity"] = {
        "step_export_generation": "step-export-71",
        "decoder_step_export_generation": "step-export-71",
        "assembly_generation": "assembly-71",
        "product_id_assembly_generation": "assembly-71",
        "color_assembly_generation": "assembly-71",
        "hierarchy_assembly_generation": "assembly-71",
        "location_assembly_generation": "assembly-71",
        "product_ids": ["root", "rotor", "housing"],
        "decoded_product_ids": ["root", "rotor", "housing"],
        "parent_product_ids": ["", "root", "root"],
        "decoded_parent_product_ids": ["", "root", "root"],
        "colors_rgb": [[180, 180, 180], [220, 40, 40], [80, 100, 140]],
        "decoded_colors_rgb": [[180, 180, 180], [220, 40, 40], [80, 100, 140]],
        "component_location_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "decoded_component_location_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "assembly_metadata_sha256": "8" * 64,
        "decoded_assembly_metadata_sha256": "8" * 64,
    }
    identity["sketch_constraint_entity_id_solver_order_generation_identity"] = {
        "sketch_generation": "sketch-71",
        "entity_table_sketch_generation": "sketch-71",
        "constraint_table_sketch_generation": "sketch-71",
        "solver_order_sketch_generation": "sketch-71",
        "entity_ids": [101, 102, 103],
        "replay_entity_ids": [101, 102, 103],
        "constraint_ids": [201, 202],
        "solver_constraint_order": [201, 202],
        "replay_solver_constraint_order": [201, 202],
        "constraint_entity_ids": [[101, 102], [102, 103]],
        "replay_constraint_entity_ids": [[101, 102], [102, 103]],
        "entity_table_sha256": "9" * 64,
        "replay_entity_table_sha256": "9" * 64,
        "constraint_table_sha256": "a" * 64,
        "replay_constraint_table_sha256": "a" * 64,
    }
    return row


def _public_v22():
    reference, measured = _public_v21()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["joint_connector_frame_labeled_face_subshape_generation_identity"] = {
                "shape_generation": "shape-81",
                "label_table_shape_generation": "shape-81",
                "connector_shape_generation": "shape-81",
                "location_shape_generation": "shape-81",
                "labeled_face_subshape_id": "face:mount_a",
                "resolved_labeled_face_subshape_id": "face:mount_a",
                "labeled_face_geometry_sha256": digest,
                "resolved_labeled_face_geometry_sha256": digest,
                "connector_origin": [0.0, 0.0, 0.0],
                "evaluated_connector_origin": [0.0, 0.0, 0.0],
                "connector_axis": [0.0, 0.0, 1.0],
                "evaluated_connector_axis": [0.0, 0.0, 1.0],
                "parent_location_sha256": "3" * 64,
                "evaluated_parent_location_sha256": "3" * 64,
                "connector_frame_sha256": "4" * 64,
                "evaluated_connector_frame_sha256": "4" * 64,
            }
            row[
                "inertia_tensor_principal_axes_density_unit_location_generation_identity"
            ] = {
                "shape_generation": "shape-81",
                "density_shape_generation": "shape-81",
                "mass_property_shape_generation": "shape-81",
                "location_shape_generation": "shape-81",
                "principal_axis_shape_generation": "shape-81",
                "density_value": 7800.0,
                "evaluated_density_value": 7800.0,
                "density_unit": "kg/m^3",
                "evaluated_density_unit": "kg/m^3",
                "shape_location_sha256": "5" * 64,
                "evaluated_shape_location_sha256": "5" * 64,
                "center_of_mass": [0.01, 0.02, 0.03],
                "evaluated_center_of_mass": [0.01, 0.02, 0.03],
                "inertia_tensor_sha256": "6" * 64,
                "evaluated_inertia_tensor_sha256": "6" * 64,
                "principal_axes_sha256": "7" * 64,
                "evaluated_principal_axes_sha256": "7" * 64,
                "mass_property_sha256": digest,
                "evaluated_mass_property_sha256": digest,
            }
    return reference, measured


def _source_v22():
    row = _source_v21()
    identity = row["replay_identity"]
    identity["step_ap242_component_transform_name_product_generation_identity"] = {
        "step_export_generation": "step-export-81",
        "decoder_step_export_generation": "step-export-81",
        "assembly_generation": "assembly-81",
        "product_assembly_generation": "assembly-81",
        "name_assembly_generation": "assembly-81",
        "transform_assembly_generation": "assembly-81",
        "component_product_ids": ["root", "rotor", "housing"],
        "decoded_component_product_ids": ["root", "rotor", "housing"],
        "component_names": ["machine", "rotor", "housing"],
        "decoded_component_names": ["machine", "rotor", "housing"],
        "nested_transform_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "decoded_nested_transform_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "ap242_product_map_sha256": "4" * 64,
        "decoded_ap242_product_map_sha256": "4" * 64,
    }
    identity["curved_mesh_export_edge_chord_surface_label_generation_identity"] = {
        "shape_generation": "shape-81",
        "mesh_shape_generation": "shape-81",
        "edge_curve_shape_generation": "shape-81",
        "surface_label_shape_generation": "shape-81",
        "mesh_export_generation": "mesh-export-81",
        "metric_mesh_export_generation": "mesh-export-81",
        "edge_curve_sha256": ["5" * 64, "6" * 64],
        "exported_edge_curve_sha256": ["5" * 64, "6" * 64],
        "chordal_tolerance": 1.0e-4,
        "evaluated_chordal_tolerance": 1.0e-4,
        "length_unit": "m",
        "evaluated_length_unit": "m",
        "boundary_surface_labels": ["outer", "interface"],
        "exported_boundary_surface_labels": ["outer", "interface"],
        "curved_mesh_sha256": "7" * 64,
        "exported_curved_mesh_sha256": "7" * 64,
    }
    return row


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


def _public_v26():
    reference, measured = _public_v25()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_identity"] = {
                "feature_generation": "fillet-131", "selector_feature_generation": "fillet-131",
                "topology_feature_generation": "fillet-131", "tolerance_feature_generation": "fillet-131",
                "shape_feature_generation": "fillet-131", "result_feature_generation": "fillet-131",
                "feature_type": "fillet", "result_feature_type": "fillet",
                "edge_selector_names": ["outer-top-1", "outer-top-2"],
                "result_edge_selector_names": ["outer-top-1", "outer-top-2"],
                "persistent_edge_ids": [101, 102], "result_persistent_edge_ids": [101, 102],
                "feature_radius_m": 0.002, "result_feature_radius_m": 0.002,
                "linear_tolerance_m": 1.0e-7, "result_linear_tolerance_m": 1.0e-7,
                "topology_signature": {"solids": 1, "faces": 18, "edges": 36},
                "result_topology_signature": {"solids": 1, "faces": 18, "edges": 36},
                "input_shape_sha256": "3" * 64, "result_input_shape_sha256": "3" * 64,
                "feature_shape_sha256": digest, "result_feature_shape_sha256": digest,
            }
            row["mass_density_center_inertia_reference_frame_assembly_generation_identity"] = {
                "mass_generation": "mass-131", "density_mass_generation": "mass-131",
                "center_mass_generation": "mass-131", "inertia_mass_generation": "mass-131",
                "frame_mass_generation": "mass-131", "assembly_mass_generation": "mass-131",
                "result_mass_generation": "mass-131", "density_kg_m3": 7800.0,
                "result_density_kg_m3": 7800.0, "volume_m3": 0.001, "result_volume_m3": 0.001,
                "mass_kg": 7.8, "result_mass_kg": 7.8,
                "center_of_mass_m": [0.01, 0.02, 0.03], "result_center_of_mass_m": [0.01, 0.02, 0.03],
                "inertia_tensor_kg_m2": [[0.2, 0.01, 0.0], [0.01, 0.3, 0.02], [0.0, 0.02, 0.4]],
                "result_inertia_tensor_kg_m2": [[0.2, 0.01, 0.0], [0.01, 0.3, 0.02], [0.0, 0.02, 0.4]],
                "reference_frame": "assembly-root", "result_reference_frame": "assembly-root",
                "assembly_transform_sha256": "4" * 64, "result_assembly_transform_sha256": "4" * 64,
                "mass_property_sha256": digest, "result_mass_property_sha256": digest,
            }
    return reference, measured


def _source_v26():
    row = _source_v25()
    identity = row["replay_identity"]
    identity["stl_chord_tolerance_triangle_normal_orientation_component_digest_generation_identity"] = {
        "stl_generation": "stl-131", "tolerance_stl_generation": "stl-131",
        "triangle_stl_generation": "stl-131", "normal_stl_generation": "stl-131",
        "component_stl_generation": "stl-131", "result_stl_generation": "stl-131",
        "linear_deflection_m": 0.0001, "decoded_linear_deflection_m": 0.0001,
        "angular_deflection_rad": 0.1, "decoded_angular_deflection_rad": 0.1,
        "triangle_count": 240, "decoded_triangle_count": 240,
        "normal_orientation": "outward", "decoded_normal_orientation": "outward",
        "component_ids": ["base", "rotor", "cover"], "decoded_component_ids": ["base", "rotor", "cover"],
        "source_shape_sha256": "5" * 64, "tessellated_source_shape_sha256": "5" * 64,
        "stl_sha256": "6" * 64, "decoded_stl_sha256": "6" * 64,
    }
    identity["builder_context_workplane_local_frame_part_identity_cache_generation_identity"] = {
        "context_generation": "builder-131", "stack_context_generation": "builder-131",
        "workplane_context_generation": "builder-131", "frame_context_generation": "builder-131",
        "part_context_generation": "builder-131", "cache_context_generation": "builder-131",
        "result_context_generation": "builder-131", "context_stack": ["BuildPart", "BuildSketch", "Locations"],
        "result_context_stack": ["BuildPart", "BuildSketch", "Locations"],
        "workplane_origin_m": [0.0, 0.0, 0.01], "result_workplane_origin_m": [0.0, 0.0, 0.01],
        "workplane_normal": [0.0, 0.0, 1.0], "result_workplane_normal": [0.0, 0.0, 1.0],
        "local_frame_transform": [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.01]],
        "result_local_frame_transform": [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.01]],
        "part_ids": ["base", "rotor", "cover"], "result_part_ids": ["base", "rotor", "cover"],
        "builder_cache_sha256": "7" * 64, "result_builder_cache_sha256": "7" * 64,
        "builder_result_sha256": "8" * 64, "result_builder_result_sha256": "8" * 64,
    }
    return row


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


def _public_v28():
    reference, measured = _public_v27()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["shell_offset_face_normal_thickness_join_self_intersection_mass_generation_identity"] = {
                "shell_generation": "shell-offset-151", "face_shell_generation": "shell-offset-151",
                "normal_shell_generation": "shell-offset-151", "thickness_shell_generation": "shell-offset-151",
                "join_shell_generation": "shell-offset-151", "intersection_shell_generation": "shell-offset-151",
                "shape_shell_generation": "shell-offset-151", "mass_shell_generation": "shell-offset-151",
                "result_shell_generation": "shell-offset-151",
                "selected_face_names": ["top-face", "bottom-face"],
                "result_selected_face_names": ["top-face", "bottom-face"],
                "normal_direction": "outward", "result_normal_direction": "outward",
                "thickness_m": 0.002, "result_thickness_m": 0.002,
                "join_mode": "arc", "result_join_mode": "arc",
                "self_intersection": False, "result_self_intersection": False,
                "shape_generation": "shell-shape-151", "result_shape_generation": "shell-shape-151",
                "is_valid_solid": True, "result_is_valid_solid": True,
                "volume_m3": 0.00042, "result_volume_m3": 0.00042,
                "mass_property_sha256": digest, "result_mass_property_sha256": digest,
            }
            row["path_sweep_frame_transition_profile_orientation_solid_volume_generation_identity"] = {
                "sweep_generation": "path-sweep-151", "path_sweep_generation": "path-sweep-151",
                "frame_sweep_generation": "path-sweep-151", "transition_sweep_generation": "path-sweep-151",
                "profile_sweep_generation": "path-sweep-151", "orientation_sweep_generation": "path-sweep-151",
                "solid_sweep_generation": "path-sweep-151", "volume_sweep_generation": "path-sweep-151",
                "result_sweep_generation": "path-sweep-151",
                "path_edge_names": ["path-0", "path-1", "path-2"],
                "result_path_edge_names": ["path-0", "path-1", "path-2"],
                "moving_frame": "parallel_transport", "result_moving_frame": "parallel_transport",
                "transition_mode": "round", "result_transition_mode": "round",
                "profile_wire_names": ["profile-outer"], "result_profile_wire_names": ["profile-outer"],
                "profile_orientation_deg": [0.0, 12.0, 24.0],
                "result_profile_orientation_deg": [0.0, 12.0, 24.0],
                "is_valid_solid": True, "result_is_valid_solid": True,
                "volume_m3": 0.00031, "result_volume_m3": 0.00031,
                "sweep_shape_sha256": digest, "result_sweep_shape_sha256": digest,
            }
    return reference, measured


def _source_v28():
    row = _source_v27()
    identity = row["replay_identity"]
    identity["step_ap242_context_product_uuid_unit_color_placement_shape_file_generation_identity"] = {
        "step_generation": "ap242-import-151", "context_step_generation": "ap242-import-151",
        "product_step_generation": "ap242-import-151", "unit_step_generation": "ap242-import-151",
        "color_step_generation": "ap242-import-151", "placement_step_generation": "ap242-import-151",
        "shape_step_generation": "ap242-import-151", "file_step_generation": "ap242-import-151",
        "result_step_generation": "ap242-import-151", "representation_context": "mechanical_design_3d",
        "decoded_representation_context": "mechanical_design_3d",
        "product_uuid_map": [["base", "11111111-1111-4111-8111-111111111111"], ["cover", "22222222-2222-4222-8222-222222222222"]],
        "decoded_product_uuid_map": [["base", "11111111-1111-4111-8111-111111111111"], ["cover", "22222222-2222-4222-8222-222222222222"]],
        "length_unit": "mm", "decoded_length_unit": "mm",
        "colors_rgb": [["base", 0.8, 0.1, 0.1], ["cover", 0.1, 0.1, 0.8]],
        "decoded_colors_rgb": [["base", 0.8, 0.1, 0.1], ["cover", 0.1, 0.1, 0.8]],
        "placements": [["base", 0.0, 0.0, 0.0], ["cover", 0.0, 0.0, 10.0]],
        "decoded_placements": [["base", 0.0, 0.0, 0.0], ["cover", 0.0, 0.0, 10.0]],
        "shape_owner_map": [["base-shape", "base"], ["cover-shape", "cover"]],
        "decoded_shape_owner_map": [["base-shape", "base"], ["cover-shape", "cover"]],
        "step_sha256": "3" * 64, "decoded_step_sha256": "3" * 64,
        "shape_map_sha256": "4" * 64, "decoded_shape_map_sha256": "4" * 64,
    }
    identity["occt_kernel_shape_location_tolerance_triangulation_serialization_cache_generation_identity"] = {
        "cache_generation": "occt-cache-151", "kernel_cache_generation": "occt-cache-151",
        "shape_cache_generation": "occt-cache-151", "location_cache_generation": "occt-cache-151",
        "tolerance_cache_generation": "occt-cache-151", "triangulation_cache_generation": "occt-cache-151",
        "serialization_cache_generation": "occt-cache-151", "result_cache_generation": "occt-cache-151",
        "kernel_version": "OCCT-7.8.1", "decoded_kernel_version": "OCCT-7.8.1",
        "shape_sha256": "5" * 64, "decoded_shape_sha256": "5" * 64,
        "location_transform": [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03]],
        "decoded_location_transform": [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03]],
        "modeling_tolerance_m": 1.0e-7, "decoded_modeling_tolerance_m": 1.0e-7,
        "triangulation_parameters": [0.001, 0.2, 1.0],
        "decoded_triangulation_parameters": [0.001, 0.2, 1.0],
        "triangulation_sha256": "6" * 64, "decoded_triangulation_sha256": "6" * 64,
        "serialization_format": "brep-binary-v1", "decoded_serialization_format": "brep-binary-v1",
        "serialization_sha256": "7" * 64, "decoded_serialization_sha256": "7" * 64,
        "cache_sha256": "8" * 64, "decoded_cache_sha256": "8" * 64,
    }
    return row


def _public_v29():
    reference, measured = _public_v28()
    neutral_radius = 0.003 + 0.42 * 0.0015
    bend_allowance = math.radians(90.0) * neutral_radius
    flat_area = 0.05 * (0.1 + 0.08 + bend_allowance)
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            sheet_digest = ("1" if index == 0 else "2") * 64
            joint_digest = ("3" if index == 0 else "4") * 64
            generation = "sheet-metal-161"
            row[
                "sheet_metal_bend_allowance_kfactor_neutral_axis_relief_thickness_flat_pattern_area_generation_identity"
            ] = {
                "sheet_generation": generation,
                "bend_sheet_generation": generation,
                "neutral_axis_sheet_generation": generation,
                "relief_sheet_generation": generation,
                "thickness_sheet_generation": generation,
                "pattern_sheet_generation": generation,
                "area_sheet_generation": generation,
                "result_sheet_generation": generation,
                "bend_radius_m": 0.003,
                "result_bend_radius_m": 0.003,
                "bend_angle_deg": 90.0,
                "result_bend_angle_deg": 90.0,
                "k_factor": 0.42,
                "result_k_factor": 0.42,
                "neutral_axis_radius_m": neutral_radius,
                "result_neutral_axis_radius_m": neutral_radius,
                "bend_allowance_m": bend_allowance,
                "result_bend_allowance_m": bend_allowance,
                "relief_type": "rectangular",
                "result_relief_type": "rectangular",
                "relief_width_m": 0.002,
                "result_relief_width_m": 0.002,
                "thickness_m": 0.0015,
                "result_thickness_m": 0.0015,
                "strip_width_m": 0.05,
                "result_strip_width_m": 0.05,
                "straight_lengths_m": [0.1, 0.08],
                "result_straight_lengths_m": [0.1, 0.08],
                "flat_pattern_area_m2": flat_area,
                "result_flat_pattern_area_m2": flat_area,
                "flat_pattern_wire_closed": True,
                "result_flat_pattern_wire_closed": True,
                "flat_pattern_sha256": sheet_digest,
                "result_flat_pattern_sha256": sheet_digest,
            }
            generation = "joint-loop-161"
            row[
                "joint_kinematic_loop_graph_dof_limit_connector_frame_closure_configuration_swept_volume_generation_identity"
            ] = {
                "joint_generation": generation,
                "graph_joint_generation": generation,
                "dof_joint_generation": generation,
                "limit_joint_generation": generation,
                "frame_joint_generation": generation,
                "closure_joint_generation": generation,
                "configuration_joint_generation": generation,
                "swept_joint_generation": generation,
                "result_joint_generation": generation,
                "joint_graph_edges": [
                    ["ground", "revolute-a", "link-a"],
                    ["link-a", "prismatic-b", "slider-b"],
                    ["slider-b", "fixed-c", "ground"],
                ],
                "result_joint_graph_edges": [
                    ["ground", "revolute-a", "link-a"],
                    ["link-a", "prismatic-b", "slider-b"],
                    ["slider-b", "fixed-c", "ground"],
                ],
                "dof_names": ["theta_a_deg", "travel_b_m"],
                "result_dof_names": ["theta_a_deg", "travel_b_m"],
                "dof_types": ["revolute", "prismatic"],
                "result_dof_types": ["revolute", "prismatic"],
                "lower_limits": [-30.0, 0.0],
                "result_lower_limits": [-30.0, 0.0],
                "upper_limits": [30.0, 0.02],
                "result_upper_limits": [30.0, 0.02],
                "configuration_values": [10.0, 0.012],
                "result_configuration_values": [10.0, 0.012],
                "connector_frame_sha256": "5" * 64,
                "result_connector_frame_sha256": "5" * 64,
                "loop_closure_error_m": 2.0e-10,
                "result_loop_closure_error_m": 2.0e-10,
                "loop_closure_tolerance_m": 1.0e-8,
                "result_loop_closure_tolerance_m": 1.0e-8,
                "configuration_id": "pose-10deg-12mm",
                "result_configuration_id": "pose-10deg-12mm",
                "swept_volume_m3": 0.00073,
                "result_swept_volume_m3": 0.00073,
                "swept_shape_sha256": joint_digest,
                "result_swept_shape_sha256": joint_digest,
            }
    return reference, measured


def _source_v29():
    row = _source_v28()
    identity = row["replay_identity"]
    generation = "dxf-face-161"
    identity[
        "dxf_arc_spline_layer_plane_unit_closed_wire_orientation_face_digest_generation_identity"
    ] = {
        "dxf_generation": generation,
        "arc_dxf_generation": generation,
        "spline_dxf_generation": generation,
        "layer_dxf_generation": generation,
        "plane_dxf_generation": generation,
        "unit_dxf_generation": generation,
        "wire_dxf_generation": generation,
        "face_dxf_generation": generation,
        "result_dxf_generation": generation,
        "arc_parameters": [["arc-1", 0.0, 0.0, 10.0, 0.0, 180.0]],
        "decoded_arc_parameters": [["arc-1", 0.0, 0.0, 10.0, 0.0, 180.0]],
        "spline_parameters": [["spline-1", 3, 4, "6" * 64]],
        "decoded_spline_parameters": [["spline-1", 3, 4, "6" * 64]],
        "entity_layer_map": [["arc-1", "profile"], ["spline-1", "profile"]],
        "decoded_entity_layer_map": [["arc-1", "profile"], ["spline-1", "profile"]],
        "workplane_matrix": [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "decoded_workplane_matrix": [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "wire_closed": True,
        "decoded_wire_closed": True,
        "wire_orientation": "counterclockwise",
        "decoded_wire_orientation": "counterclockwise",
        "face_area_mm2": 314.1592653589793,
        "decoded_face_area_mm2": 314.1592653589793,
        "dxf_sha256": "7" * 64,
        "decoded_dxf_sha256": "7" * 64,
        "face_sha256": "8" * 64,
        "decoded_face_sha256": "8" * 64,
    }
    generation = "3mf-import-161"
    identity[
        "three_mf_component_transform_triangle_winding_material_watertight_volume_unit_digest_generation_identity"
    ] = {
        "three_mf_generation": generation,
        "component_three_mf_generation": generation,
        "transform_three_mf_generation": generation,
        "winding_three_mf_generation": generation,
        "material_three_mf_generation": generation,
        "watertight_three_mf_generation": generation,
        "volume_three_mf_generation": generation,
        "file_three_mf_generation": generation,
        "result_three_mf_generation": generation,
        "component_names": ["housing", "insert"],
        "decoded_component_names": ["housing", "insert"],
        "component_transform_sha256": [["housing", "9" * 64], ["insert", "a" * 64]],
        "decoded_component_transform_sha256": [["housing", "9" * 64], ["insert", "a" * 64]],
        "triangle_winding": "outward_counterclockwise",
        "decoded_triangle_winding": "outward_counterclockwise",
        "material_id_map": [["housing", 1], ["insert", 2]],
        "decoded_material_id_map": [["housing", 1], ["insert", 2]],
        "watertight_component_names": ["housing", "insert"],
        "decoded_watertight_component_names": ["housing", "insert"],
        "component_volumes_mm3": [["housing", 1200.0], ["insert", 300.0]],
        "decoded_component_volumes_mm3": [["housing", 1200.0], ["insert", 300.0]],
        "total_volume_mm3": 1500.0,
        "decoded_total_volume_mm3": 1500.0,
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "three_mf_sha256": "b" * 64,
        "decoded_three_mf_sha256": "b" * 64,
        "triangle_mesh_sha256": "c" * 64,
        "decoded_triangle_mesh_sha256": "c" * 64,
    }
    return row


def _public_v30():
    reference, measured = _public_v29()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            shape_digest = ("1" if index == 0 else "2") * 64
            profile_digest = ("3" if index == 0 else "4") * 64
            generation = "helical-sweep-171"
            volume = 1.2e-6 + index * 0.3e-6
            centroid = [0.0, 0.0, 0.02 + index * 0.005]
            frame = [[1.0, 0.0, 0.0, 0.005], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
            row[
                "helical_sweep_pitch_handedness_profile_frame_turn_self_intersection_volume_centroid_shape_generation_identity"
            ] = {
                "helical_generation": generation,
                "pitch_helical_generation": generation,
                "handedness_helical_generation": generation,
                "profile_helical_generation": generation,
                "turn_helical_generation": generation,
                "intersection_helical_generation": generation,
                "mass_helical_generation": generation,
                "shape_helical_generation": generation,
                "result_helical_generation": generation,
                "pitch_m": 0.01,
                "result_pitch_m": 0.01,
                "handedness": "right",
                "result_handedness": "right",
                "profile_frame_matrix": frame,
                "result_profile_frame_matrix": frame,
                "profile_frame_sha256": profile_digest,
                "result_profile_frame_sha256": profile_digest,
                "turn_count": 4.0,
                "result_turn_count": 4.0,
                "self_intersection": False,
                "result_self_intersection": False,
                "volume_m3": volume,
                "result_volume_m3": volume,
                "centroid_m": centroid,
                "result_centroid_m": centroid,
                "helical_shape_sha256": shape_digest,
                "result_helical_shape_sha256": shape_digest,
            }
            generation = "boolean-history-171"
            boolean_digest = ("5" if index == 0 else "6") * 64
            history_digest = ("7" if index == 0 else "8") * 64
            inertia = [[1.0e-9, 0.0, 0.0], [0.0, 2.0e-9, 0.0], [0.0, 0.0, 3.0e-9]]
            row[
                "boolean_tolerance_operand_order_history_volume_centroid_inertia_shape_generation_identity"
            ] = {
                "boolean_generation": generation,
                "tolerance_boolean_generation": generation,
                "operand_boolean_generation": generation,
                "history_boolean_generation": generation,
                "mass_boolean_generation": generation,
                "inertia_boolean_generation": generation,
                "shape_boolean_generation": generation,
                "result_boolean_generation": generation,
                "operation": "cut",
                "result_operation": "cut",
                "model_tolerance_m": 1.0e-7,
                "result_model_tolerance_m": 1.0e-7,
                "operand_order": ["body-a", "tool-b"],
                "result_operand_order": ["body-a", "tool-b"],
                "subshape_history_sha256": history_digest,
                "result_subshape_history_sha256": history_digest,
                "volume_m3": volume,
                "result_volume_m3": volume,
                "centroid_m": centroid,
                "result_centroid_m": centroid,
                "inertia_tensor_kg_m2": inertia,
                "result_inertia_tensor_kg_m2": inertia,
                "boolean_shape_sha256": boolean_digest,
                "result_boolean_shape_sha256": boolean_digest,
            }
    return reference, measured


def _source_v30():
    row = _source_v29()
    identity = row["replay_identity"]
    generation = "step-assembly-171"
    identity[
        "step_assembly_product_color_instance_transform_unit_hierarchy_file_generation_identity"
    ] = {
        "step_generation": generation,
        "product_step_generation": generation,
        "color_step_generation": generation,
        "instance_step_generation": generation,
        "transform_step_generation": generation,
        "unit_step_generation": generation,
        "hierarchy_step_generation": generation,
        "file_step_generation": generation,
        "result_step_generation": generation,
        "product_names": ["assembly", "housing", "shaft"],
        "decoded_product_names": ["assembly", "housing", "shaft"],
        "product_colors_rgb": [["housing", 0.2, 0.4, 0.8], ["shaft", 0.7, 0.7, 0.7]],
        "decoded_product_colors_rgb": [["housing", 0.2, 0.4, 0.8], ["shaft", 0.7, 0.7, 0.7]],
        "instance_order": ["housing-1", "shaft-1"],
        "decoded_instance_order": ["housing-1", "shaft-1"],
        "instance_transform_sha256": [["housing-1", "9" * 64], ["shaft-1", "a" * 64]],
        "decoded_instance_transform_sha256": [["housing-1", "9" * 64], ["shaft-1", "a" * 64]],
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "assembly_hierarchy": [["assembly", "housing-1"], ["assembly", "shaft-1"]],
        "decoded_assembly_hierarchy": [["assembly", "housing-1"], ["assembly", "shaft-1"]],
        "step_sha256": "b" * 64,
        "decoded_step_sha256": "b" * 64,
        "assembly_shape_sha256": "c" * 64,
        "decoded_assembly_shape_sha256": "c" * 64,
    }
    generation = "stl-solid-171"
    identity[
        "stl_facet_normal_winding_watertight_tolerance_volume_unit_file_generation_identity"
    ] = {
        "stl_generation": generation,
        "normal_stl_generation": generation,
        "winding_stl_generation": generation,
        "watertight_stl_generation": generation,
        "tolerance_stl_generation": generation,
        "volume_stl_generation": generation,
        "unit_stl_generation": generation,
        "file_stl_generation": generation,
        "result_stl_generation": generation,
        "facet_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
        "decoded_facet_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
        "triangle_winding": "outward_counterclockwise",
        "decoded_triangle_winding": "outward_counterclockwise",
        "unmatched_edge_count": 0,
        "decoded_unmatched_edge_count": 0,
        "merge_tolerance_m": 1.0e-7,
        "decoded_merge_tolerance_m": 1.0e-7,
        "signed_volume_m3": 1.5e-6,
        "decoded_signed_volume_m3": 1.5e-6,
        "length_unit": "m",
        "decoded_length_unit": "m",
        "stl_sha256": "d" * 64,
        "decoded_stl_sha256": "d" * 64,
        "stl_solid_sha256": "e" * 64,
        "decoded_stl_solid_sha256": "e" * 64,
    }
    return row


def _public_v31():
    reference, measured = _public_v30()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            generation = "assembly-mass-181"
            shape_digest = ("1" if index == 0 else "2") * 64
            locations = [
                [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
                [[1.0, 0.0, 0.0, 0.1], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            ]
            inertia = [
                [0.01, 0.0, 0.0],
                [0.0, 0.0205, 0.0],
                [0.0, 0.0, 0.0205],
            ]
            row[
                "assembly_occurrence_location_density_unit_suppression_mass_center_inertia_parallel_axis_shape_generation_identity"
            ] = {
                "assembly_generation": generation,
                "occurrence_assembly_generation": generation,
                "location_assembly_generation": generation,
                "density_assembly_generation": generation,
                "suppression_assembly_generation": generation,
                "mass_assembly_generation": generation,
                "inertia_assembly_generation": generation,
                "shape_assembly_generation": generation,
                "result_assembly_generation": generation,
                "occurrence_ids": ["housing-1", "shaft-1"],
                "result_occurrence_ids": ["housing-1", "shaft-1"],
                "location_matrices": locations,
                "result_location_matrices": locations,
                "densities_kg_m3": [2700.0, 7800.0],
                "result_densities_kg_m3": [2700.0, 7800.0],
                "density_unit": "kg/m^3",
                "result_density_unit": "kg/m^3",
                "suppressed_occurrence_ids": [],
                "result_suppressed_occurrence_ids": [],
                "part_volumes_m3": [1.0e-3, 2.0e-4],
                "result_part_volumes_m3": [1.0e-3, 2.0e-4],
                "assembly_mass_kg": 4.26,
                "result_assembly_mass_kg": 4.26,
                "center_of_mass_m": [0.03661971830985916, 0.0, 0.0],
                "result_center_of_mass_m": [0.03661971830985916, 0.0, 0.0],
                "inertia_reference_frame": "assembly_global_center_of_mass",
                "result_inertia_reference_frame": "assembly_global_center_of_mass",
                "parallel_axis_applied": True,
                "result_parallel_axis_applied": True,
                "assembly_inertia_kg_m2": inertia,
                "result_assembly_inertia_kg_m2": inertia,
                "assembly_shape_sha256": shape_digest,
                "result_assembly_shape_sha256": shape_digest,
            }
            generation = "loft-lineage-181"
            lineage_digest = ("3" if index == 0 else "4") * 64
            loft_digest = ("5" if index == 0 else "6") * 64
            lineage = [
                ["section-0:e1", "loft:f1"],
                ["section-1:e1", "loft:f1"],
                ["section-2:e1", "loft:f1"],
            ]
            row[
                "loft_sweep_profile_order_seam_guide_orientation_face_lineage_shell_volume_shape_generation_identity"
            ] = {
                "loft_generation": generation,
                "profile_loft_generation": generation,
                "seam_loft_generation": generation,
                "guide_loft_generation": generation,
                "lineage_loft_generation": generation,
                "shell_loft_generation": generation,
                "mass_loft_generation": generation,
                "shape_loft_generation": generation,
                "result_loft_generation": generation,
                "profile_ids": ["section-0", "section-1", "section-2"],
                "result_profile_ids": ["section-0", "section-1", "section-2"],
                "profile_parameters": [0.0, 0.5, 1.0],
                "result_profile_parameters": [0.0, 0.5, 1.0],
                "seam_parameters": [0.0, 0.0, 0.0],
                "result_seam_parameters": [0.0, 0.0, 0.0],
                "guide_orientation": "start_to_end_right_handed",
                "result_guide_orientation": "start_to_end_right_handed",
                "face_lineage_pairs": lineage,
                "result_face_lineage_pairs": lineage,
                "face_lineage_sha256": lineage_digest,
                "result_face_lineage_sha256": lineage_digest,
                "shell_closed": True,
                "result_shell_closed": True,
                "volume_m3": 1.5e-3,
                "result_volume_m3": 1.5e-3,
                "loft_shape_sha256": loft_digest,
                "result_loft_shape_sha256": loft_digest,
            }
    return reference, measured


def _source_v31():
    row = _source_v30()
    identity = row["replay_identity"]
    generation = "step-ap242-context-181"
    identity[
        "step_ap242_representation_context_external_owner_occurrence_transform_mass_product_structure_file_generation_identity"
    ] = {
        "step_generation": generation,
        "schema_step_generation": generation,
        "context_step_generation": generation,
        "owner_step_generation": generation,
        "occurrence_step_generation": generation,
        "transform_step_generation": generation,
        "mass_step_generation": generation,
        "structure_step_generation": generation,
        "file_step_generation": generation,
        "result_step_generation": generation,
        "ap_schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF",
        "decoded_ap_schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF",
        "representation_contexts": [["ctx-mm", "mm", "degree"], ["ctx-m", "m", "radian"]],
        "decoded_representation_contexts": [["ctx-mm", "mm", "degree"], ["ctx-m", "m", "radian"]],
        "part_owners": [["housing", "part-housing-v3"], ["shaft", "part-shaft-v5"]],
        "decoded_part_owners": [["housing", "part-housing-v3"], ["shaft", "part-shaft-v5"]],
        "occurrence_ids": ["housing-1", "shaft-1", "shaft-2"],
        "decoded_occurrence_ids": ["housing-1", "shaft-1", "shaft-2"],
        "occurrence_transform_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64], ["shaft-2", "9" * 64]],
        "decoded_occurrence_transform_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64], ["shaft-2", "9" * 64]],
        "occurrence_mass_properties": [["housing-1", 2.7, 0.0, 0.0, 0.0], ["shaft-1", 0.78, 0.1, 0.0, 0.0], ["shaft-2", 0.78, -0.1, 0.0, 0.0]],
        "decoded_occurrence_mass_properties": [["housing-1", 2.7, 0.0, 0.0, 0.0], ["shaft-1", 0.78, 0.1, 0.0, 0.0], ["shaft-2", 0.78, -0.1, 0.0, 0.0]],
        "product_structure_sha256": "a" * 64,
        "decoded_product_structure_sha256": "a" * 64,
        "step_file_sha256": "b" * 64,
        "decoded_step_file_sha256": "b" * 64,
    }
    generation = "stl-tessellation-error-181"
    identity[
        "stl_tessellation_source_brep_chord_angle_facet_component_deviation_area_volume_digest_generation_identity"
    ] = {
        "stl_generation": generation,
        "source_stl_generation": generation,
        "tolerance_stl_generation": generation,
        "facet_stl_generation": generation,
        "component_stl_generation": generation,
        "deviation_stl_generation": generation,
        "area_stl_generation": generation,
        "volume_stl_generation": generation,
        "file_stl_generation": generation,
        "result_stl_generation": generation,
        "source_brep_sha256": "c" * 64,
        "decoded_source_brep_sha256": "c" * 64,
        "chord_tolerance_m": 1.0e-4,
        "decoded_chord_tolerance_m": 1.0e-4,
        "angular_tolerance_deg": 10.0,
        "decoded_angular_tolerance_deg": 10.0,
        "facet_count": 1200,
        "decoded_facet_count": 1200,
        "connected_component_count": 1,
        "decoded_connected_component_count": 1,
        "maximum_surface_deviation_m": 8.0e-5,
        "decoded_maximum_surface_deviation_m": 8.0e-5,
        "source_surface_area_m2": 0.02,
        "decoded_surface_area_m2": 0.01999,
        "source_volume_m3": 1.5e-3,
        "decoded_volume_m3": 1.4995e-3,
        "relative_area_tolerance": 1.0e-3,
        "relative_volume_tolerance": 1.0e-3,
        "normal_table_sha256": "d" * 64,
        "decoded_normal_table_sha256": "d" * 64,
        "stl_file_sha256": "e" * 64,
        "decoded_stl_file_sha256": "e" * 64,
    }
    return row


def _public_v32():
    reference, measured = _public_v31()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            suffix = "1" if index == 0 else "2"
            generation = "boolean-topology-closure-191"
            row[
                "boolean_fuzzy_tolerance_topology_name_face_ancestry_count_volume_centroid_shape_generation_identity"
            ] = {
                "boolean_generation": generation,
                "tolerance_generation": generation,
                "topology_generation": generation,
                "ancestry_generation": generation,
                "mass_generation": generation,
                "shape_generation": generation,
                "result_generation": generation,
                "fuzzy_tolerance_m": 1.0e-7,
                "result_fuzzy_tolerance_m": 1.0e-7,
                "surviving_topology_names": ["body", "cut_face", "outer_shell"],
                "result_surviving_topology_names": ["body", "cut_face", "outer_shell"],
                "face_ancestry": [["box:f1", "result:f3"], ["cylinder:f2", "result:f7"]],
                "result_face_ancestry": [["box:f1", "result:f3"], ["cylinder:f2", "result:f7"]],
                "solid_count": 1,
                "result_solid_count": 1,
                "volume_m3": 9.0e-4,
                "result_volume_m3": 9.0e-4,
                "centroid_m": [0.01, 0.0, 0.0],
                "result_centroid_m": [0.01, 0.0, 0.0],
                "boolean_shape_sha256": suffix * 64,
                "result_boolean_shape_sha256": suffix * 64,
            }
            generation = "sweep-frame-closure-191"
            sweep_suffix = "3" if index == 0 else "4"
            row[
                "sweep_frenet_frame_twist_transition_profile_orientation_self_intersection_volume_owner_shape_generation_identity"
            ] = {
                "sweep_generation": generation,
                "frame_generation": generation,
                "twist_generation": generation,
                "transition_generation": generation,
                "orientation_generation": generation,
                "intersection_generation": generation,
                "mass_generation": generation,
                "owner_generation": generation,
                "result_generation": generation,
                "frame_convention": "corrected_frenet",
                "result_frame_convention": "corrected_frenet",
                "twist_parameters_rad": [0.0, 0.2, 0.4],
                "result_twist_parameters_rad": [0.0, 0.2, 0.4],
                "transition_mode": "right_corner",
                "result_transition_mode": "right_corner",
                "profile_orientation_signs": [1, 1, 1],
                "result_profile_orientation_signs": [1, 1, 1],
                "self_intersection": False,
                "result_self_intersection": False,
                "sweep_volume_m3": 1.2e-3,
                "result_sweep_volume_m3": 1.2e-3,
                "shape_owner": "sweep/body1",
                "result_shape_owner": "sweep/body1",
                "sweep_shape_sha256": sweep_suffix * 64,
                "result_sweep_shape_sha256": sweep_suffix * 64,
            }
    return reference, measured


def _source_v32():
    row = _source_v31()
    identity = row["replay_identity"]
    generation = "step-assembly-closure-191"
    identity[
        "step_assembly_instance_transform_unit_color_material_uuid_component_volume_file_generation_identity"
    ] = {
        "step_generation": generation,
        "instance_generation": generation,
        "transform_generation": generation,
        "unit_generation": generation,
        "color_generation": generation,
        "material_generation": generation,
        "uuid_generation": generation,
        "component_generation": generation,
        "volume_generation": generation,
        "file_generation": generation,
        "result_generation": generation,
        "instance_ids": ["housing-1", "shaft-1"],
        "decoded_instance_ids": ["housing-1", "shaft-1"],
        "instance_transform_sha256": [["housing-1", "5" * 64], ["shaft-1", "6" * 64]],
        "decoded_instance_transform_sha256": [["housing-1", "5" * 64], ["shaft-1", "6" * 64]],
        "length_unit": "m",
        "decoded_length_unit": "m",
        "instance_colors_rgb": [["housing-1", 0.8, 0.8, 0.8], ["shaft-1", 0.3, 0.3, 0.3]],
        "decoded_instance_colors_rgb": [["housing-1", 0.8, 0.8, 0.8], ["shaft-1", 0.3, 0.3, 0.3]],
        "material_labels": [["housing-1", "aluminum"], ["shaft-1", "steel"]],
        "decoded_material_labels": [["housing-1", "aluminum"], ["shaft-1", "steel"]],
        "product_uuids": [["housing-1", "uuid-housing-v4"], ["shaft-1", "uuid-shaft-v6"]],
        "decoded_product_uuids": [["housing-1", "uuid-housing-v4"], ["shaft-1", "uuid-shaft-v6"]],
        "component_shape_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64]],
        "decoded_component_shape_sha256": [["housing-1", "7" * 64], ["shaft-1", "8" * 64]],
        "total_volume_m3": 1.2e-3,
        "decoded_total_volume_m3": 1.2e-3,
        "step_file_sha256": "9" * 64,
        "decoded_step_file_sha256": "9" * 64,
    }
    generation = "stl-repair-closure-191"
    identity[
        "stl_repair_merge_tolerance_normal_duplicate_boundary_watertight_volume_unit_file_generation_identity"
    ] = {
        "repair_generation": generation,
        "tolerance_generation": generation,
        "normal_generation": generation,
        "duplicate_generation": generation,
        "boundary_generation": generation,
        "watertight_generation": generation,
        "volume_generation": generation,
        "unit_generation": generation,
        "file_generation": generation,
        "result_generation": generation,
        "merge_tolerance_m": 1.0e-6,
        "decoded_merge_tolerance_m": 1.0e-6,
        "facet_normal_sha256": "a" * 64,
        "decoded_facet_normal_sha256": "a" * 64,
        "duplicate_vertex_count": 0,
        "decoded_duplicate_vertex_count": 0,
        "boundary_edge_count": 0,
        "decoded_boundary_edge_count": 0,
        "watertight_component_count": 1,
        "decoded_watertight_component_count": 1,
        "repaired_volume_m3": 1.5e-3,
        "decoded_repaired_volume_m3": 1.5e-3,
        "length_unit": "m",
        "decoded_length_unit": "m",
        "source_stl_sha256": "b" * 64,
        "decoded_source_stl_sha256": "b" * 64,
        "repaired_stl_sha256": "c" * 64,
        "decoded_repaired_stl_sha256": "c" * 64,
    }
    return row


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


def _public_v34():
    reference, measured = _public_v33()
    identity4 = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    rotations = [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ]
    inertia_local = [
        [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]],
        [[0.4, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.6]],
        [[0.7, 0.0, 0.0], [0.0, 0.8, 0.0], [0.0, 0.0, 0.9]],
    ]
    inertia_global = [
        [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]],
        [[0.5, 0.0, 0.0], [0.0, 0.4, 0.0], [0.0, 0.0, 0.6]],
        [[0.7, 0.0, 0.0], [0.0, 0.8, 0.0], [0.0, 0.0, 0.9]],
    ]
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            suffix = "1" if index == 0 else "2"
            generation = "assembly-mate-211"
            row[
                "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
            ] = {
                "assembly_generation": generation,
                **{
                    key: generation
                    for key in (
                        "mate_generation",
                        "cycle_generation",
                        "frame_generation",
                        "mass_generation",
                        "center_generation",
                        "inertia_generation",
                        "owner_generation",
                        "shape_generation",
                        "result_generation",
                    )
                },
                "part_ids": ["base", "arm", "tool"],
                "result_part_ids": ["base", "arm", "tool"],
                "mate_edges": [["base", "arm"], ["arm", "tool"], ["tool", "base"]],
                "result_mate_edges": [["base", "arm"], ["arm", "tool"], ["tool", "base"]],
                "mate_cycle_transform": identity4,
                "result_mate_cycle_transform": identity4,
                "frame_determinants": [1.0, 1.0, 1.0],
                "result_frame_determinants": [1.0, 1.0, 1.0],
                "part_masses_kg": [2.0, 3.0, 1.0],
                "result_part_masses_kg": [2.0, 3.0, 1.0],
                "part_centers_m": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                "result_part_centers_m": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                "assembly_mass_kg": 6.0,
                "result_assembly_mass_kg": 6.0,
                "assembly_center_of_mass_m": [2.0 / 3.0, 1.0 / 6.0, 0.0],
                "result_assembly_center_of_mass_m": [2.0 / 3.0, 1.0 / 6.0, 0.0],
                "part_rotation_matrices": rotations,
                "result_part_rotation_matrices": rotations,
                "part_inertia_local_kg_m2": inertia_local,
                "result_part_inertia_local_kg_m2": inertia_local,
                "part_inertia_global_kg_m2": inertia_global,
                "result_part_inertia_global_kg_m2": inertia_global,
                "assembly_owner": "assembly/root",
                "result_assembly_owner": "assembly/root",
                "assembly_shape_sha256": suffix * 64,
                "accepted_assembly_shape_sha256": suffix * 64,
            }
            generation = "shell-fillet-211"
            shell_suffix = "3" if index == 0 else "4"
            row[
                "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
            ] = {
                "shell_fillet_generation": generation,
                **{
                    key: generation
                    for key in (
                        "topology_generation",
                        "thickness_generation",
                        "volume_generation",
                        "area_generation",
                        "inertia_generation",
                        "convergence_generation",
                        "brep_generation",
                        "result_generation",
                    )
                },
                "vertex_count": 16,
                "result_vertex_count": 16,
                "edge_count": 24,
                "result_edge_count": 24,
                "face_count": 10,
                "result_face_count": 10,
                "euler_characteristic": 2,
                "result_euler_characteristic": 2,
                "edge_face_incidence_counts": [2] * 24,
                "result_edge_face_incidence_counts": [2] * 24,
                "nominal_wall_thickness_m": 0.002,
                "result_nominal_wall_thickness_m": 0.002,
                "wall_thickness_samples_m": [0.002, 0.0020001, 0.0019999],
                "result_wall_thickness_samples_m": [0.002, 0.0020001, 0.0019999],
                "original_volume_m3": 1.0,
                "result_original_volume_m3": 1.0,
                "removed_volume_m3": 0.8,
                "result_removed_volume_m3": 0.8,
                "shell_volume_m3": 0.2,
                "result_shell_volume_m3": 0.2,
                "surface_area_m2": 5.0,
                "result_surface_area_m2": 5.0,
                "inertia_tensor_kg_m2": [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.4]],
                "result_inertia_tensor_kg_m2": [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.4]],
                "convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6],
                "result_convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6],
                "convergence_volumes_m3": [0.2001, 0.200001, 0.2],
                "result_convergence_volumes_m3": [0.2001, 0.200001, 0.2],
                "shell_brep_sha256": shell_suffix * 64,
                "accepted_shell_brep_sha256": shell_suffix * 64,
            }
    return reference, measured


def _source_v34():
    row = _source_v33()
    identity = row["replay_identity"]
    generation = "step-semantic-roundtrip-211"
    identity[
        "step_unit_product_entity_color_assembly_transform_shape_validity_owner_file_result_generation_identity"
    ] = {
        "step_generation": generation,
        **{
            key: generation
            for key in (
                "unit_generation",
                "entity_generation",
                "color_generation",
                "transform_generation",
                "shape_generation",
                "validity_generation",
                "owner_generation",
                "file_generation",
                "result_generation",
            )
        },
        "length_unit": "m",
        "decoded_length_unit": "m",
        "product_entities": [[1, "base"], [2, "arm"]],
        "decoded_product_entities": [[1, "base"], [2, "arm"]],
        "entity_colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]],
        "decoded_entity_colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]],
        "assembly_transform_sha256": [[1, "5" * 64], [2, "6" * 64]],
        "decoded_assembly_transform_sha256": [[1, "5" * 64], [2, "6" * 64]],
        "shape_count": 2,
        "decoded_shape_count": 2,
        "solid_validity": [[1, True], [2, True]],
        "decoded_solid_validity": [[1, True], [2, True]],
        "source_owner": "assembly/root",
        "decoded_source_owner": "assembly/root",
        "step_file_sha256": "7" * 64,
        "decoded_step_file_sha256": "7" * 64,
    }
    generation = "brep-periodic-seam-211"
    identity[
        "brep_periodic_face_seam_edge_orientation_pcurve_tolerance_manifold_serializer_shape_result_generation_identity"
    ] = {
        "periodic_brep_generation": generation,
        **{
            key: generation
            for key in (
                "face_generation",
                "seam_generation",
                "orientation_generation",
                "pcurve_generation",
                "tolerance_generation",
                "manifold_generation",
                "serializer_generation",
                "shape_generation",
                "result_generation",
            )
        },
        "periodic_face_ids": [1, 2],
        "decoded_periodic_face_ids": [1, 2],
        "seam_edge_multiplicity": [[1, 11, 2], [2, 12, 2]],
        "decoded_seam_edge_multiplicity": [[1, 11, 2], [2, 12, 2]],
        "face_orientation_signs": [[1, 1], [2, -1]],
        "decoded_face_orientation_signs": [[1, 1], [2, -1]],
        "edge_pcurve_max_deviation_m": [[11, 1.0e-8], [12, 2.0e-8]],
        "decoded_edge_pcurve_max_deviation_m": [[11, 1.0e-8], [12, 2.0e-8]],
        "vertex_tolerances_m": [[101, 1.0e-7], [102, 2.0e-7]],
        "decoded_vertex_tolerances_m": [[101, 1.0e-7], [102, 2.0e-7]],
        "edge_face_incidence_counts": [[11, 2], [12, 2]],
        "decoded_edge_face_incidence_counts": [[11, 2], [12, 2]],
        "serializer_version": "occt-brep-v3",
        "decoded_serializer_version": "occt-brep-v3",
        "periodic_brep_shape_sha256": "8" * 64,
        "decoded_periodic_brep_shape_sha256": "8" * 64,
    }
    return row


def _public_v35():
    reference, measured = _public_v34()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            suffix = "1" if index == 0 else "2"
            generation = "mirrored-pattern-221"
            row["mirrored_pattern_occurrence_handedness_transform_suppression_volume_mass_center_inertia_owner_shape_result_generation_identity"] = {
                "pattern_generation": generation,
                **{key: generation for key in ("occurrence_generation", "handedness_generation", "transform_generation", "suppression_generation", "volume_generation", "mass_generation", "inertia_generation", "owner_generation", "shape_generation", "result_generation")},
                "occurrence_ids": ["part:0", "part:mirror_x"], "result_occurrence_ids": ["part:0", "part:mirror_x"],
                "transform_determinants": [1.0, -1.0], "result_transform_determinants": [1.0, -1.0],
                "handedness": ["right", "left"], "result_handedness": ["right", "left"],
                "suppressed": [False, False], "result_suppressed": [False, False],
                "occurrence_volumes_m3": [1.0, 1.0], "result_occurrence_volumes_m3": [1.0, 1.0],
                "occurrence_masses_kg": [2.0, 2.0], "result_occurrence_masses_kg": [2.0, 2.0],
                "occurrence_centers_m": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], "result_occurrence_centers_m": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
                "occurrence_inertia_principal_kg_m2": [[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]], "result_occurrence_inertia_principal_kg_m2": [[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]],
                "assembly_volume_m3": 2.0, "result_assembly_volume_m3": 2.0,
                "assembly_mass_kg": 4.0, "result_assembly_mass_kg": 4.0,
                "assembly_center_m": [0.0, 0.0, 0.0], "result_assembly_center_m": [0.0, 0.0, 0.0],
                "assembly_owner": "pattern/root", "result_assembly_owner": "pattern/root",
                "mirrored_shape_sha256": suffix * 64, "accepted_mirrored_shape_sha256": suffix * 64,
            }
            generation = "offset-thicken-221"
            digest = ("3" if index == 0 else "4") * 64
            row["offset_thicken_curvature_sign_selfintersection_repair_thickness_volume_topology_convergence_brep_result_generation_identity"] = {
                "thicken_generation": generation,
                **{key: generation for key in ("curvature_generation", "offset_generation", "intersection_generation", "repair_generation", "thickness_generation", "volume_generation", "topology_generation", "convergence_generation", "brep_generation", "result_generation")},
                "minimum_curvature_radius_m": 0.05, "result_minimum_curvature_radius_m": 0.05,
                "offset_m": -0.01, "result_offset_m": -0.01,
                "self_intersection_count": 0, "result_self_intersection_count": 0,
                "repair_mode": "none_required", "result_repair_mode": "none_required",
                "wall_thickness_samples_m": [0.01, 0.010001, 0.009999], "result_wall_thickness_samples_m": [0.01, 0.010001, 0.009999],
                "original_volume_m3": 1.0, "result_original_volume_m3": 1.0,
                "removed_volume_m3": 0.7, "result_removed_volume_m3": 0.7,
                "thickened_volume_m3": 0.3, "result_thickened_volume_m3": 0.3,
                "solid_count": 1, "result_solid_count": 1, "shell_count": 1, "result_shell_count": 1,
                "convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6], "result_convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6],
                "convergence_volumes_m3": [0.3001, 0.300001, 0.3], "result_convergence_volumes_m3": [0.3001, 0.300001, 0.3],
                "thicken_brep_sha256": digest, "accepted_thicken_brep_sha256": digest,
            }
    return reference, measured


def _source_v35():
    row = _source_v34()
    identity = row["replay_identity"]
    generation = "step-ap242-pmi-221"
    identity["step_ap242_pmi_unit_name_color_occurrence_transform_validity_owner_file_result_generation_identity"] = {
        "ap242_generation": generation,
        **{key: generation for key in ("pmi_generation", "unit_generation", "name_generation", "color_generation", "occurrence_generation", "transform_generation", "validity_generation", "owner_generation", "file_generation", "result_generation")},
        "schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING", "decoded_schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING",
        "length_unit": "m", "decoded_length_unit": "m",
        "pmi_annotations": [[1, "linear_dimension", "mm", 25.0], [2, "diameter", "mm", 10.0]], "decoded_pmi_annotations": [[1, "linear_dimension", "mm", 25.0], [2, "diameter", "mm", 10.0]],
        "product_names": [[1, "base"], [2, "arm"]], "decoded_product_names": [[1, "base"], [2, "arm"]],
        "colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]], "decoded_colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]],
        "occurrence_paths": [[1, "root/base"], [2, "root/arm"]], "decoded_occurrence_paths": [[1, "root/base"], [2, "root/arm"]],
        "transform_sha256": [[1, "5" * 64], [2, "6" * 64]], "decoded_transform_sha256": [[1, "5" * 64], [2, "6" * 64]],
        "solid_validity": [[1, True], [2, True]], "decoded_solid_validity": [[1, True], [2, True]],
        "source_owner": "assembly/root", "decoded_source_owner": "assembly/root",
        "ap242_file_sha256": "7" * 64, "decoded_ap242_file_sha256": "7" * 64,
    }
    generation = "dxf-profile-221"
    identity["dxf_profile_unit_plane_layer_arc_bulge_loop_winding_extrusion_topology_owner_file_result_generation_identity"] = {
        "dxf_generation": generation,
        **{key: generation for key in ("unit_generation", "plane_generation", "layer_generation", "arc_generation", "loop_generation", "winding_generation", "extrusion_generation", "topology_generation", "owner_generation", "file_generation", "result_generation")},
        "length_unit": "mm", "decoded_length_unit": "mm", "sketch_plane": "XY", "decoded_sketch_plane": "XY",
        "layer_names": ["OUTER", "HOLES"], "decoded_layer_names": ["OUTER", "HOLES"],
        "arc_bulges": [[11, 0.41421356237], [12, -0.41421356237]], "decoded_arc_bulges": [[11, 0.41421356237], [12, -0.41421356237]],
        "loop_winding_signs": [["outer", 1], ["hole", -1]], "decoded_loop_winding_signs": [["outer", 1], ["hole", -1]],
        "closed_loop_count": 2, "decoded_closed_loop_count": 2,
        "extrusion_height_m": 0.01, "decoded_extrusion_height_m": 0.01,
        "profile_area_m2": 0.1, "decoded_profile_area_m2": 0.1,
        "solid_count": 1, "decoded_solid_count": 1,
        "extruded_volume_m3": 0.001, "decoded_extruded_volume_m3": 0.001,
        "profile_owner": "dxf/profile1", "decoded_profile_owner": "dxf/profile1",
        "dxf_file_sha256": "8" * 64, "decoded_dxf_file_sha256": "8" * 64,
    }
    return row


def _public_v36():
    reference, measured = _public_v35()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "fillet-chamfer-contract-231"
            row[
                "fillet_chamfer_edge_selection_radius_distance_tolerance_euler_volume_area_owner_brep_result_generation_identity"
            ] = {
                "feature_generation": generation,
                **{key: generation for key in ("selection_generation", "size_generation", "tolerance_generation", "topology_generation", "volume_generation", "area_generation", "owner_generation", "brep_generation", "result_generation")},
                "feature_kind": "fillet",
                "result_feature_kind": "fillet",
                "selected_edge_ids": [11, 12, 13, 14],
                "result_selected_edge_ids": [11, 12, 13, 14],
                "radius_m": 0.01,
                "result_radius_m": 0.01,
                "distance_m": 0.0,
                "result_distance_m": 0.0,
                "modeling_tolerance_m": 1.0e-6,
                "result_modeling_tolerance_m": 1.0e-6,
                "topology_before_v_e_f": [8, 12, 6],
                "result_topology_before_v_e_f": [8, 12, 6],
                "topology_after_v_e_f": [16, 24, 10],
                "result_topology_after_v_e_f": [16, 24, 10],
                "volume_before_m3": 1.0,
                "result_volume_before_m3": 1.0,
                "volume_after_m3": 0.99,
                "result_volume_after_m3": 0.99,
                "surface_area_before_m2": 6.0,
                "result_surface_area_before_m2": 6.0,
                "surface_area_after_m2": 5.95,
                "result_surface_area_after_m2": 5.95,
                "shape_owner": "part:fillet31",
                "result_shape_owner": "part:fillet31",
                "feature_brep_sha256": suffix * 64,
                "accepted_feature_brep_sha256": suffix * 64,
            }
            generation = "loft-contract-231"
            row[
                "loft_section_order_orientation_guide_continuity_volume_centroid_owner_brep_result_generation_identity"
            ] = {
                "loft_generation": generation,
                **{key: generation for key in ("section_generation", "orientation_generation", "guide_generation", "continuity_generation", "volume_generation", "centroid_generation", "owner_generation", "brep_generation", "result_generation")},
                "section_ids": ["wire:z0", "wire:z1", "wire:z2"],
                "result_section_ids": ["wire:z0", "wire:z1", "wire:z2"],
                "section_parameters": [0.0, 0.5, 1.0],
                "result_section_parameters": [0.0, 0.5, 1.0],
                "section_orientation_signs": [1, 1, 1],
                "result_section_orientation_signs": [1, 1, 1],
                "guide_correspondence": [[0, 0], [1, 1], [2, 2], [3, 3]],
                "result_guide_correspondence": [[0, 0], [1, 1], [2, 2], [3, 3]],
                "continuity": "G1",
                "result_continuity": "G1",
                "loft_volume_m3": 0.75,
                "result_loft_volume_m3": 0.75,
                "loft_centroid_m": [0.0, 0.0, 0.5],
                "result_loft_centroid_m": [0.0, 0.0, 0.5],
                "loft_owner": "part:loft31",
                "result_loft_owner": "part:loft31",
                "loft_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_loft_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v36():
    row = _source_v35()
    identity = row["replay_identity"]
    generation = "step-assembly-contract-231"
    identity[
        "step_assembly_hierarchy_occurrence_transform_repeated_part_name_color_unit_owner_file_result_generation_identity"
    ] = {
        "assembly_generation": generation,
        **{key: generation for key in ("hierarchy_generation", "occurrence_generation", "transform_generation", "part_generation", "name_generation", "color_generation", "unit_generation", "owner_generation", "file_generation", "result_generation")},
        "occurrence_paths": ["root/base", "root/arm:0", "root/arm:1"],
        "decoded_occurrence_paths": ["root/base", "root/arm:0", "root/arm:1"],
        "part_ids": ["base", "arm", "arm"],
        "decoded_part_ids": ["base", "arm", "arm"],
        "transform_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "decoded_transform_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "product_names": ["base", "arm", "arm"],
        "decoded_product_names": ["base", "arm", "arm"],
        "colors_rgb": [[0.2, 0.3, 0.4], [0.8, 0.1, 0.1], [0.8, 0.1, 0.1]],
        "decoded_colors_rgb": [[0.2, 0.3, 0.4], [0.8, 0.1, 0.1], [0.8, 0.1, 0.1]],
        "length_unit": "m",
        "decoded_length_unit": "m",
        "assembly_owner": "assembly:root31",
        "decoded_assembly_owner": "assembly:root31",
        "step_file_sha256": "8" * 64,
        "decoded_step_file_sha256": "8" * 64,
    }
    generation = "boolean-history-contract-231"
    identity[
        "boolean_tolerance_healing_sliver_nonmanifold_operation_history_input_output_owner_brep_result_generation_identity"
    ] = {
        "boolean_generation": generation,
        **{key: generation for key in ("tolerance_generation", "healing_generation", "sliver_generation", "manifold_generation", "history_generation", "input_generation", "output_generation", "owner_generation", "brep_generation", "result_generation")},
        "operation": "fuse",
        "decoded_operation": "fuse",
        "fuzzy_tolerance_m": 1.0e-6,
        "decoded_fuzzy_tolerance_m": 1.0e-6,
        "healing_actions": ["same_domain", "sew"],
        "decoded_healing_actions": ["same_domain", "sew"],
        "sliver_face_count": 0,
        "decoded_sliver_face_count": 0,
        "nonmanifold_edge_count": 0,
        "decoded_nonmanifold_edge_count": 0,
        "input_shape_ids": ["solid:a", "solid:b"],
        "decoded_input_shape_ids": ["solid:a", "solid:b"],
        "output_shape_id": "solid:fused",
        "decoded_output_shape_id": "solid:fused",
        "operation_history": [["solid:a", "solid:fused"], ["solid:b", "solid:fused"]],
        "decoded_operation_history": [["solid:a", "solid:fused"], ["solid:b", "solid:fused"]],
        "input_owner": "boolean:inputs31",
        "decoded_input_owner": "boolean:inputs31",
        "output_owner": "boolean:result31",
        "decoded_output_owner": "boolean:result31",
        "boolean_brep_sha256": "9" * 64,
        "decoded_boolean_brep_sha256": "9" * 64,
    }
    return row


def _public_v37():
    reference, measured = _public_v36()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1); generation = "mass-placement-contract-241"
            row["mass_properties_centroid_inertia_principal_axes_placement_density_shape_brep_generation_identity"] = {
                "mass_generation": generation, **{key: generation for key in ("density_generation", "centroid_generation", "inertia_generation", "principal_generation", "placement_generation", "owner_generation", "brep_generation", "result_generation")},
                "density_kg_m3": 7800.0, "result_density_kg_m3": 7800.0, "volume_m3": 1.0e-3, "result_volume_m3": 1.0e-3, "mass_kg": 7.8, "result_mass_kg": 7.8,
                "centroid_world_m": [1.0, 2.0, 3.0], "result_centroid_world_m": [1.0, 2.0, 3.0],
                "inertia_world_kg_m2": [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]], "result_inertia_world_kg_m2": [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
                "principal_moments_kg_m2": [1.0, 2.0, 3.0], "result_principal_moments_kg_m2": [1.0, 2.0, 3.0],
                "principal_axes_world": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "result_principal_axes_world": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "placement_transform": [[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 2.0], [0.0, 0.0, 1.0, 3.0], [0.0, 0.0, 0.0, 1.0]], "result_placement_transform": [[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 2.0], [0.0, 0.0, 1.0, 3.0], [0.0, 0.0, 0.0, 1.0]],
                "shape_owner": "part:placed-solid-241", "result_shape_owner": "part:placed-solid-241", "shape_brep_sha256": suffix * 64, "accepted_shape_brep_sha256": suffix * 64,
            }
            generation = "shell-offset-contract-241"
            row["shell_offset_thickness_normal_side_removed_face_topology_volume_input_owner_brep_generation_identity"] = {
                "shell_generation": generation, **{key: generation for key in ("thickness_generation", "normal_generation", "side_generation", "removed_generation", "topology_generation", "volume_generation", "owner_generation", "brep_generation", "result_generation")},
                "thickness_m": 0.002, "result_thickness_m": 0.002, "offset_side": "inside", "result_offset_side": "inside",
                "face_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]], "result_face_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
                "removed_face_ids": [6], "result_removed_face_ids": [6], "wall_topology_v_e_f": [16, 24, 10], "result_wall_topology_v_e_f": [16, 24, 10],
                "analytical_shell_volume_m3": 0.012, "result_analytical_shell_volume_m3": 0.012, "input_owner": "part:shell-input-241", "result_input_owner": "part:shell-input-241",
                "shell_brep_sha256": ("3" if index == 0 else "4") * 64, "accepted_shell_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v37():
    row = _source_v36(); identity = row["replay_identity"]; generation = "sketch-solve-contract-241"
    identity["sketch_constraint_dof_solver_reference_unit_owner_source_result_generation_identity"] = {
        "sketch_generation": generation, **{key: generation for key in ("constraint_generation", "dof_generation", "solver_generation", "reference_generation", "unit_generation", "owner_generation", "source_generation", "result_generation")},
        "constraint_ids": ["horizontal:1", "vertical:2", "distance:3"], "solved_constraint_ids": ["horizontal:1", "vertical:2", "distance:3"], "remaining_dof": 0, "solved_remaining_dof": 0,
        "solver_status": "fully_constrained", "solved_solver_status": "fully_constrained", "reference_geometry_ids": ["axis:x", "origin"], "solved_reference_geometry_ids": ["axis:x", "origin"],
        "length_unit": "m", "solved_length_unit": "m", "sketch_owner": "sketch:profile-241", "solved_sketch_owner": "sketch:profile-241", "sketch_source_sha256": "5" * 64, "solved_sketch_source_sha256": "5" * 64, "sketch_result_sha256": "6" * 64, "accepted_sketch_result_sha256": "6" * 64,
    }
    generation = "toponame-contract-241"
    identity["topological_naming_edge_face_history_ocp_selector_shape_feature_source_brep_generation_identity"] = {
        "toponame_generation": generation, **{key: generation for key in ("edge_generation", "face_generation", "history_generation", "ocp_generation", "selector_generation", "shape_generation", "owner_generation", "source_generation", "brep_generation", "result_generation")},
        "edge_names": ["edge:fillet:0", "edge:fillet:1"], "replayed_edge_names": ["edge:fillet:0", "edge:fillet:1"], "face_names": ["face:top", "face:side"], "replayed_face_names": ["face:top", "face:side"],
        "operation_history": ["box", "fillet", "select:face:top"], "replayed_operation_history": ["box", "fillet", "select:face:top"], "ocp_version": "7.8.1", "replayed_ocp_version": "7.8.1", "selector_result": ["face:top"], "replayed_selector_result": ["face:top"],
        "shape_generation_id": 41, "replayed_shape_generation_id": 41, "feature_owner": "feature:fillet-241", "replayed_feature_owner": "feature:fillet-241", "feature_source_sha256": "7" * 64, "replayed_feature_source_sha256": "7" * 64, "feature_brep_sha256": "8" * 64, "accepted_feature_brep_sha256": "8" * 64,
    }
    return row


def _public_v38():
    reference, measured = _public_v37()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "revolve-contract-258"
            pappus_volume = 1.0e-2 * 2.0 * math.pi * 5.0e-2
            row[
                "revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_shape_brep_generation_identity"
            ] = {
                "revolve_generation": generation,
                **{
                    key: generation
                    for key in (
                        "axis_generation",
                        "angle_generation",
                        "profile_generation",
                        "orientation_generation",
                        "volume_generation",
                        "centroid_generation",
                        "topology_generation",
                        "owner_generation",
                        "brep_generation",
                        "result_generation",
                    )
                },
                "axis_origin_m": [0.0, 0.0, 0.0],
                "result_axis_origin_m": [0.0, 0.0, 0.0],
                "axis_direction": [0.0, 0.0, 1.0],
                "result_axis_direction": [0.0, 0.0, 1.0],
                "sweep_angle_deg": 360.0,
                "result_sweep_angle_deg": 360.0,
                "profile_axis_crossing": False,
                "result_profile_axis_crossing": False,
                "profile_loop_orientation": "counterclockwise_outward",
                "result_profile_loop_orientation": "counterclockwise_outward",
                "profile_area_m2": 1.0e-2,
                "result_profile_area_m2": 1.0e-2,
                "profile_centroid_radius_m": 5.0e-2,
                "result_profile_centroid_radius_m": 5.0e-2,
                "analytic_volume_m3": pappus_volume,
                "result_analytic_volume_m3": pappus_volume,
                "volume_tolerance_m3": 1.0e-12,
                "result_volume_tolerance_m3": 1.0e-12,
                "centroid_world_m": [0.0, 0.0, 0.1],
                "result_centroid_world_m": [0.0, 0.0, 0.1],
                "solid_count": 1,
                "result_solid_count": 1,
                "shell_count": 1,
                "result_shell_count": 1,
                "boundary_genus": 1,
                "result_boundary_genus": 1,
                "boundary_euler_characteristic": 0,
                "result_boundary_euler_characteristic": 0,
                "shape_owner": "part:revolve-258",
                "result_shape_owner": "part:revolve-258",
                "revolve_brep_sha256": suffix * 64,
                "accepted_revolve_brep_sha256": suffix * 64,
            }

            generation = "involute-gear-contract-258"
            module_m = 2.0e-3
            teeth = 20
            pressure_angle_deg = 20.0
            pitch_diameter_m = module_m * teeth
            base_diameter_m = pitch_diameter_m * math.cos(math.radians(pressure_angle_deg))
            addendum_diameter_m = pitch_diameter_m + 2.0 * module_m
            row[
                "involute_gear_module_teeth_pressure_backlash_pitch_base_addendum_periodicity_volume_shape_brep_generation_identity"
            ] = {
                "gear_generation": generation,
                **{
                    key: generation
                    for key in (
                        "module_generation",
                        "tooth_generation",
                        "pressure_generation",
                        "backlash_generation",
                        "diameter_generation",
                        "periodicity_generation",
                        "volume_generation",
                        "owner_generation",
                        "brep_generation",
                        "result_generation",
                    )
                },
                "module_m": module_m,
                "result_module_m": module_m,
                "tooth_count": teeth,
                "result_tooth_count": teeth,
                "pressure_angle_deg": pressure_angle_deg,
                "result_pressure_angle_deg": pressure_angle_deg,
                "backlash_m": 1.0e-4,
                "result_backlash_m": 1.0e-4,
                "pitch_diameter_m": pitch_diameter_m,
                "result_pitch_diameter_m": pitch_diameter_m,
                "base_diameter_m": base_diameter_m,
                "result_base_diameter_m": base_diameter_m,
                "addendum_diameter_m": addendum_diameter_m,
                "result_addendum_diameter_m": addendum_diameter_m,
                "tooth_period_angle_deg": 360.0 / teeth,
                "result_tooth_period_angle_deg": 360.0 / teeth,
                "gear_volume_m3": 8.0e-6,
                "result_gear_volume_m3": 8.0e-6,
                "shape_owner": "part:involute-gear-258",
                "result_shape_owner": "part:involute-gear-258",
                "gear_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_gear_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v38():
    row = _source_v37()
    identity = row["replay_identity"]
    generation = "brep-roundtrip-contract-258"
    identity[
        "brep_roundtrip_tolerance_ocp_version_subshape_count_bounds_volume_shape_owner_source_restored_generation_identity"
    ] = {
        "brep_generation": generation,
        **{
            key: generation
            for key in (
                "tolerance_generation",
                "ocp_generation",
                "topology_generation",
                "bounds_generation",
                "volume_generation",
                "owner_generation",
                "source_generation",
                "restored_generation",
                "result_generation",
            )
        },
        "serialization_tolerance_m": 1.0e-7,
        "restored_serialization_tolerance_m": 1.0e-7,
        "ocp_version": "7.8.1",
        "restored_ocp_version": "7.8.1",
        "subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12, "vertex": 8},
        "restored_subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12, "vertex": 8},
        "bounding_box_min_m": [-0.5, -0.25, 0.0],
        "restored_bounding_box_min_m": [-0.5, -0.25, 0.0],
        "bounding_box_max_m": [0.5, 0.25, 0.1],
        "restored_bounding_box_max_m": [0.5, 0.25, 0.1],
        "volume_m3": 5.0e-2,
        "restored_volume_m3": 5.0e-2,
        "shape_owner": "part:brep-roundtrip-258",
        "restored_shape_owner": "part:brep-roundtrip-258",
        "source_brep_sha256": "5" * 64,
        "restored_source_brep_sha256": "5" * 64,
        "restored_brep_sha256": "6" * 64,
        "accepted_restored_brep_sha256": "6" * 64,
    }

    generation = "svg-extrusion-contract-258"
    identity[
        "svg_path_fillrule_curve_transform_unit_wire_face_extrusion_source_digest_generation_identity"
    ] = {
        "svg_generation": generation,
        **{
            key: generation
            for key in (
                "path_generation",
                "fill_generation",
                "transform_generation",
                "unit_generation",
                "wire_generation",
                "face_generation",
                "volume_generation",
                "owner_generation",
                "digest_generation",
                "result_generation",
            )
        },
        "path_commands": ["M", "C", "C", "Z"],
        "replayed_path_commands": ["M", "C", "C", "Z"],
        "fill_rule": "nonzero",
        "replayed_fill_rule": "nonzero",
        "curve_transform": [[1.0, 0.0, 10.0], [0.0, 1.0, 20.0], [0.0, 0.0, 1.0]],
        "replayed_curve_transform": [[1.0, 0.0, 10.0], [0.0, 1.0, 20.0], [0.0, 0.0, 1.0]],
        "document_unit": "mm",
        "replayed_document_unit": "mm",
        "document_to_meter_scale": 1.0e-3,
        "replayed_document_to_meter_scale": 1.0e-3,
        "wire_closed": True,
        "replayed_wire_closed": True,
        "face_orientation": "counterclockwise_positive",
        "replayed_face_orientation": "counterclockwise_positive",
        "profile_area_m2": 2.0e-3,
        "replayed_profile_area_m2": 2.0e-3,
        "extrusion_distance_m": 1.0e-2,
        "replayed_extrusion_distance_m": 1.0e-2,
        "extrusion_volume_m3": 2.0e-5,
        "replayed_extrusion_volume_m3": 2.0e-5,
        "source_owner": "svg:path-logo-258",
        "replayed_source_owner": "svg:path-logo-258",
        "svg_source_sha256": "7" * 64,
        "replayed_svg_source_sha256": "7" * 64,
        "svg_result_sha256": "8" * 64,
        "accepted_svg_result_sha256": "8" * 64,
    }
    return row


_LOFT_V39 = "loft_multiwire_hole_correspondence_section_orientation_selfintersection_volume_centroid_euler_owner_brep_generation_identity"

_OFFSET = "offset_shell_signed_thickness_join_repair_area_volume_mass_centroid_inertia_owner_brep_generation_identity"

_HEAL = "occt_heal_tolerance_sew_orientation_degenerate_stablename_roundtrip_owner_digest_generation_identity"

_ASSEMBLY_V39 = "assembly_hierarchy_location_joint_axis_collision_quantity_bom_owner_result_generation_identity"


def _generation(generation: str, *keys: str) -> dict[str, str]:
    return {key: generation for key in keys}


def _public_v39():
    reference, measured = _public_v38()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "multiwire-loft-271"
            row[_LOFT_V39] = {
                "loft_generation": generation,
                **_generation(generation, "wire_generation", "correspondence_generation", "orientation_generation", "hole_generation", "intersection_generation", "mass_generation", "topology_generation", "owner_generation", "brep_generation", "result_generation"),
                "outer_wire_ids": [11, 21, 31],
                "result_outer_wire_ids": [11, 21, 31],
                "inner_wire_ids": [12, 22, 32],
                "result_inner_wire_ids": [12, 22, 32],
                "section_order": [0, 1, 2],
                "result_section_order": [0, 1, 2],
                "outer_orientation": ["ccw", "ccw", "ccw"],
                "result_outer_orientation": ["ccw", "ccw", "ccw"],
                "inner_orientation": ["cw", "cw", "cw"],
                "result_inner_orientation": ["cw", "cw", "cw"],
                "hole_continuity": [[12, 22], [22, 32]],
                "result_hole_continuity": [[12, 22], [22, 32]],
                "self_intersection_free": True,
                "result_self_intersection_free": True,
                "volume_m3": 2.0e-3,
                "result_volume_m3": 2.0e-3,
                "centroid_m": [0.0, 0.0, 5.0e-2],
                "result_centroid_m": [0.0, 0.0, 5.0e-2],
                "solid_count": 1,
                "result_solid_count": 1,
                "boundary_euler_characteristic": 0,
                "result_boundary_euler_characteristic": 0,
                "shape_owner": "part:multiwire-loft-271",
                "result_shape_owner": "part:multiwire-loft-271",
                "loft_brep_sha256": suffix * 64,
                "accepted_loft_brep_sha256": suffix * 64,
            }
            generation = "offset-shell-271"
            row[_OFFSET] = {
                "offset_generation": generation,
                **_generation(generation, "thickness_generation", "join_generation", "repair_generation", "area_generation", "volume_generation", "mass_generation", "owner_generation", "brep_generation", "result_generation"),
                "signed_thickness_m": -2.0e-3,
                "result_signed_thickness_m": -2.0e-3,
                "offset_direction": "inward",
                "result_offset_direction": "inward",
                "join_mode": "arc",
                "result_join_mode": "arc",
                "self_intersection_detected": True,
                "result_self_intersection_detected": True,
                "self_intersection_repaired": True,
                "result_self_intersection_repaired": True,
                "outer_area_m2": 1.0,
                "result_outer_area_m2": 1.0,
                "inner_area_m2": 0.8,
                "result_inner_area_m2": 0.8,
                "enclosed_volume_m3": 1.6e-3,
                "result_enclosed_volume_m3": 1.6e-3,
                "density_kg_per_m3": 1000.0,
                "result_density_kg_per_m3": 1000.0,
                "mass_kg": 1.6,
                "result_mass_kg": 1.6,
                "centroid_m": [0.0, 0.0, 0.0],
                "result_centroid_m": [0.0, 0.0, 0.0],
                "principal_inertia_kg_m2": [0.01, 0.01, 0.02],
                "result_principal_inertia_kg_m2": [0.01, 0.01, 0.02],
                "shape_owner": "part:offset-shell-271",
                "result_shape_owner": "part:offset-shell-271",
                "offset_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_offset_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v39():
    row = _source_v38()
    identity = row["replay_identity"]
    generation = "occt-heal-271"
    identity[_HEAL] = {
        "heal_generation": generation,
        **_generation(generation, "tolerance_generation", "sew_generation", "orientation_generation", "degenerate_generation", "name_generation", "roundtrip_generation", "owner_generation", "result_generation"),
        "healing_tolerance_m": 1.0e-7,
        "replayed_healing_tolerance_m": 1.0e-7,
        "sewn_shell_count": 1,
        "replayed_sewn_shell_count": 1,
        "solid_count": 1,
        "replayed_solid_count": 1,
        "face_orientations": [1, 1, 1, 1, 1, 1],
        "replayed_face_orientations": [1, 1, 1, 1, 1, 1],
        "removed_degenerate_edge_ids": [91, 92],
        "replayed_removed_degenerate_edge_ids": [91, 92],
        "stable_subshape_names": {"face:top": 5, "face:bottom": 6},
        "replayed_stable_subshape_names": {"face:top": 5, "face:bottom": 6},
        "roundtrip_subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12},
        "replayed_roundtrip_subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12},
        "shape_owner": "headless:occt-heal-271",
        "replayed_shape_owner": "headless:occt-heal-271",
        "healed_brep_sha256": "5" * 64,
        "replayed_healed_brep_sha256": "5" * 64,
        "heal_result_sha256": "6" * 64,
        "accepted_heal_result_sha256": "6" * 64,
    }
    generation = "assembly-hierarchy-271"
    identity[_ASSEMBLY_V39] = {
        "assembly_generation": generation,
        **_generation(generation, "hierarchy_generation", "location_generation", "joint_generation", "collision_generation", "quantity_generation", "bom_generation", "owner_generation", "result_generation"),
        "hierarchy": {"root": ["frame", "shaft"], "shaft": ["rotor"]},
        "replayed_hierarchy": {"root": ["frame", "shaft"], "shaft": ["rotor"]},
        "local_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.2]},
        "replayed_local_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.2]},
        "global_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.3]},
        "replayed_global_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.3]},
        "joint_axis": [0.0, 0.0, 1.0],
        "replayed_joint_axis": [0.0, 0.0, 1.0],
        "collision_pairs": [["frame", "rotor"]],
        "replayed_collision_pairs": [["frame", "rotor"]],
        "component_quantities": {"frame": 1, "shaft": 1, "rotor": 1},
        "replayed_component_quantities": {"frame": 1, "shaft": 1, "rotor": 1},
        "bom_identity": {"frame": "FRAME-001", "shaft": "SHAFT-001", "rotor": "ROTOR-001"},
        "replayed_bom_identity": {"frame": "FRAME-001", "shaft": "SHAFT-001", "rotor": "ROTOR-001"},
        "assembly_owner": "headless:assembly-271",
        "replayed_assembly_owner": "headless:assembly-271",
        "assembly_result_sha256": "7" * 64,
        "accepted_assembly_result_sha256": "7" * 64,
    }
    return row


_DRAFT = (
    "draft_neutral_plane_pull_direction_signed_angle_face_tangent_topology_"
    "volume_centroid_owner_brep_generation_identity"
)

_THREAD = (
    "thread_pitch_handedness_flank_profile_diameter_runout_selfintersection_"
    "volume_owner_brep_generation_identity"
)

_OBB = (
    "oriented_bounding_box_principal_inertia_axis_com_frame_transform_unit_"
    "owner_result_generation_identity"
)

_TESSELLATION = (
    "tessellation_linear_angular_deflection_triangle_vertex_orientation_"
    "watertight_stl_brep_owner_result_generation_identity"
)


def _public_v40():
    reference, measured = _public_v39()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "draft-solid-311"
            row[_DRAFT] = {
                "draft_generation": generation,
                **_generation(
                    generation,
                    "plane_generation",
                    "pull_generation",
                    "angle_generation",
                    "face_generation",
                    "tangent_generation",
                    "topology_generation",
                    "mass_generation",
                    "owner_generation",
                    "brep_generation",
                    "result_generation",
                ),
                "neutral_plane_origin_m": [0.0, 0.0, 0.0],
                "result_neutral_plane_origin_m": [0.0, 0.0, 0.0],
                "neutral_plane_normal": [0.0, 0.0, 1.0],
                "result_neutral_plane_normal": [0.0, 0.0, 1.0],
                "pull_direction": [0.0, 0.0, 1.0],
                "result_pull_direction": [0.0, 0.0, 1.0],
                "signed_draft_angle_rad": 0.05235987755982988,
                "result_signed_draft_angle_rad": 0.05235987755982988,
                "selected_face_ids": [1, 2, 3, 4],
                "result_selected_face_ids": [1, 2, 3, 4],
                "tangent_continuity": True,
                "result_tangent_continuity": True,
                "topology_signature": {
                    "solid": 1,
                    "shell": 1,
                    "face": 10,
                    "edge": 24,
                    "vertex": 16,
                },
                "result_topology_signature": {
                    "solid": 1,
                    "shell": 1,
                    "face": 10,
                    "edge": 24,
                    "vertex": 16,
                },
                "volume_m3": 1.9e-3,
                "result_volume_m3": 1.9e-3,
                "centroid_m": [0.0, 0.0, 5.0e-2],
                "result_centroid_m": [0.0, 0.0, 5.0e-2],
                "shape_owner": "part:draft-solid-311",
                "result_shape_owner": "part:draft-solid-311",
                "draft_brep_sha256": suffix * 64,
                "accepted_draft_brep_sha256": suffix * 64,
            }

            generation = "modeled-thread-311"
            row[_THREAD] = {
                "thread_generation": generation,
                **_generation(
                    generation,
                    "pitch_generation",
                    "handedness_generation",
                    "profile_generation",
                    "diameter_generation",
                    "runout_generation",
                    "intersection_generation",
                    "mass_generation",
                    "owner_generation",
                    "brep_generation",
                    "result_generation",
                ),
                "pitch_m": 2.0e-3,
                "result_pitch_m": 2.0e-3,
                "handedness": "right",
                "result_handedness": "right",
                "flank_angle_rad": 1.0471975511965976,
                "result_flank_angle_rad": 1.0471975511965976,
                "profile_type": "iso_v",
                "result_profile_type": "iso_v",
                "major_diameter_m": 2.0e-2,
                "result_major_diameter_m": 2.0e-2,
                "minor_diameter_m": 1.70e-2,
                "result_minor_diameter_m": 1.70e-2,
                "thread_length_m": 2.0e-2,
                "result_thread_length_m": 2.0e-2,
                "turn_count": 10.0,
                "result_turn_count": 10.0,
                "runout_length_m": 4.0e-3,
                "result_runout_length_m": 4.0e-3,
                "self_intersection_free": True,
                "result_self_intersection_free": True,
                "volume_m3": 5.0e-6,
                "result_volume_m3": 5.0e-6,
                "shape_owner": "part:modeled-thread-311",
                "result_shape_owner": "part:modeled-thread-311",
                "thread_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_thread_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v40():
    row = _source_v39()
    identity = row["replay_identity"]
    generation = "obb-inertia-311"
    identity[_OBB] = {
        "obb_generation": generation,
        **_generation(
            generation,
            "box_generation",
            "inertia_generation",
            "axis_generation",
            "com_generation",
            "transform_generation",
            "unit_generation",
            "owner_generation",
            "result_generation",
        ),
        "obb_center_m": [1.0, 2.0, 3.0],
        "replayed_obb_center_m": [1.0, 2.0, 3.0],
        "obb_half_extents_m": [0.5, 0.25, 0.125],
        "replayed_obb_half_extents_m": [0.5, 0.25, 0.125],
        "obb_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "replayed_obb_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "principal_moments_kg_m2": [0.01, 0.02, 0.025],
        "replayed_principal_moments_kg_m2": [0.01, 0.02, 0.025],
        "principal_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "replayed_principal_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "center_of_mass_local_m": [0.0, 0.0, 0.0],
        "replayed_center_of_mass_local_m": [0.0, 0.0, 0.0],
        "center_of_mass_world_m": [1.0, 2.0, 3.0],
        "replayed_center_of_mass_world_m": [1.0, 2.0, 3.0],
        "local_to_world_transform": [
            [0.0, -1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 1.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "replayed_local_to_world_transform": [
            [0.0, -1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 1.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "length_unit": "m",
        "replayed_length_unit": "m",
        "inertia_unit": "kg*m^2",
        "replayed_inertia_unit": "kg*m^2",
        "shape_owner": "headless:obb-inertia-311",
        "replayed_shape_owner": "headless:obb-inertia-311",
        "obb_result_sha256": "5" * 64,
        "accepted_obb_result_sha256": "5" * 64,
    }

    generation = "tessellation-311"
    vertices = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
    identity[_TESSELLATION] = {
        "tessellation_generation": generation,
        **_generation(
            generation,
            "linear_generation",
            "angular_generation",
            "triangle_generation",
            "vertex_generation",
            "orientation_generation",
            "watertight_generation",
            "stl_generation",
            "brep_generation",
            "owner_generation",
            "result_generation",
        ),
        "linear_deflection_m": 1.0e-4,
        "replayed_linear_deflection_m": 1.0e-4,
        "angular_deflection_rad": 0.2,
        "replayed_angular_deflection_rad": 0.2,
        "triangle_count": 12,
        "replayed_triangle_count": 12,
        "vertex_coordinates_m": vertices,
        "replayed_vertex_coordinates_m": vertices,
        "outward_orientation": True,
        "replayed_outward_orientation": True,
        "watertight": True,
        "replayed_watertight": True,
        "nonmanifold_edge_count": 0,
        "replayed_nonmanifold_edge_count": 0,
        "stl_owner": "headless:tessellation-311",
        "replayed_stl_owner": "headless:tessellation-311",
        "source_brep_sha256": "6" * 64,
        "replayed_source_brep_sha256": "6" * 64,
        "stl_sha256": "7" * 64,
        "replayed_stl_sha256": "7" * 64,
        "tessellation_result_sha256": "8" * 64,
        "accepted_tessellation_result_sha256": "8" * 64,
    }
    return row


_LOFT_V41 = (
    "loft_section_order_parameterization_seam_twist_closure_topology_volume_"
    "centroid_owner_brep_generation_identity"
)

_SHELL = (
    "shell_offset_removedface_direction_thickness_join_selfintersection_"
    "validity_mass_owner_brep_generation_identity"
)

_ASSEMBLY_V41 = (
    "assembly_interference_part_transform_unit_contact_overlap_clearance_"
    "owner_result_generation_identity"
)

_STEP_V41 = (
    "step_name_color_layer_hierarchy_subshape_label_source_import_owner_"
    "result_generation_identity"
)


def _public_v41():
    reference, measured = _public_v40()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "loft-solid-724"
            row[_LOFT_V41] = {
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
    identity[_ASSEMBLY_V41] = {
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
    identity[_STEP_V41] = {
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


_HELIX = (
    "helicalsweep_pitch_turns_profile_frame_frenet_torsion_validity_volume_"
    "centroid_owner_brep_generation_identity"
)

_BOOLEAN_V42 = (
    "boolean_fuzzy_tolerance_sliver_topology_validity_volume_surfacearea_"
    "owner_brep_generation_identity"
)

_SELECTOR = (
    "selector_cache_topology_renumber_geometry_predicate_feature_parent_"
    "owner_result_generation_identity"
)

_STEP_V42 = (
    "step_assembly_unit_color_name_transform_shape_hierarchy_export_owner_"
    "file_generation_identity"
)


def _public_v42():
    reference, measured = _public_v41()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "helical-sweep-725"
            radius, pitch, turns, area = 0.01, 0.02, 5.0, 1.0e-4
            rise = pitch * turns
            path = turns * math.sqrt((2.0 * math.pi * radius) ** 2 + pitch**2)
            b = pitch / (2.0 * math.pi)
            torsion = b / (radius**2 + b**2)
            row[_HELIX] = {
                "helical_generation": generation,
                **_generation(generation, "pitch_generation", "turn_generation", "profile_generation", "frame_generation", "torsion_generation", "validity_generation", "mass_generation", "owner_generation", "brep_generation", "result_generation"),
                "helix_radius_m": radius, "result_helix_radius_m": radius,
                "pitch_m": pitch, "result_pitch_m": pitch,
                "turns": turns, "result_turns": turns,
                "axial_rise_m": rise, "result_axial_rise_m": rise,
                "profile_area_m2": area, "result_profile_area_m2": area,
                "profile_frame": "frenet", "result_profile_frame": "frenet",
                "frenet_transport": True, "result_frenet_transport": True,
                "path_length_m": path, "result_path_length_m": path,
                "torsion_per_m": torsion, "result_torsion_per_m": torsion,
                "solid_valid": True, "result_solid_valid": True,
                "volume_m3": area * path, "result_volume_m3": area * path,
                "centroid_m": [0.0, 0.0, rise / 2.0], "result_centroid_m": [0.0, 0.0, rise / 2.0],
                "shape_owner": "part:helical-sweep-725", "result_shape_owner": "part:helical-sweep-725",
                "helical_brep_sha256": suffix * 64, "accepted_helical_brep_sha256": suffix * 64,
            }
            generation = "fuzzy-boolean-725"
            topology = {"solid": 1, "shell": 1, "face": 10, "edge": 24, "vertex": 16}
            row[_BOOLEAN_V42] = {
                "boolean_generation": generation,
                **_generation(generation, "tolerance_generation", "sliver_generation", "topology_generation", "validity_generation", "mass_generation", "owner_generation", "brep_generation", "result_generation"),
                "operation": "cut", "result_operation": "cut",
                "fuzzy_tolerance_m": 1.0e-7, "result_fuzzy_tolerance_m": 1.0e-7,
                "minimum_feature_size_m": 1.0e-4, "result_minimum_feature_size_m": 1.0e-4,
                "sliver_face_count_before": 2, "result_sliver_face_count_before": 2,
                "sliver_face_count_after": 0, "result_sliver_face_count_after": 0,
                "topology_signature": topology, "result_topology_signature": topology,
                "solid_valid": True, "result_solid_valid": True,
                "volume_m3": 9.9e-4, "result_volume_m3": 9.9e-4,
                "surface_area_m2": 6.1e-2, "result_surface_area_m2": 6.1e-2,
                "shape_owner": "part:fuzzy-boolean-725", "result_shape_owner": "part:fuzzy-boolean-725",
                "boolean_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_boolean_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v42():
    row = _source_v41()
    identity = row["replay_identity"]
    generation = "selector-cache-725"
    renumber = {"face:11": "face:21", "face:12": "face:22"}
    identity[_SELECTOR] = {
        "selector_generation": generation,
        **_generation(generation, "topology_generation", "renumber_generation", "geometry_generation", "predicate_generation", "feature_generation", "owner_generation", "result_generation"),
        "topology_renumber_map": renumber, "replayed_topology_renumber_map": renumber,
        "geometry_generation_id": 725, "replayed_geometry_generation_id": 725,
        "selector_predicate": "Axis.Z and Area>1e-4", "replayed_selector_predicate": "Axis.Z and Area>1e-4",
        "selected_feature_ids": ["face:21", "face:22"], "replayed_selected_feature_ids": ["face:21", "face:22"],
        "parent_shape_owner": "headless:selector-parent-725", "replayed_parent_shape_owner": "headless:selector-parent-725",
        "selector_result_sha256": "5" * 64, "accepted_selector_result_sha256": "5" * 64,
    }
    generation = "step-assembly-725"
    transforms = {
        "part:1": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
        "part:2": [[1, 0, 0, 100], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    }
    names = {"part:1": "base", "part:2": "slider"}
    colors = {"part:1": [0.2, 0.4, 0.8], "part:2": [0.8, 0.4, 0.2]}
    shape_ids = {"part:1": "shape:101", "part:2": "shape:102"}
    hierarchy = {"assembly:1": ["part:1", "part:2"]}
    identity[_STEP_V42] = {
        "assembly_generation": generation,
        **_generation(generation, "unit_generation", "color_generation", "name_generation", "transform_generation", "shape_generation", "hierarchy_generation", "owner_generation", "file_generation", "result_generation"),
        "length_unit": "mm", "replayed_length_unit": "mm",
        "unit_scale_to_m": 1.0e-3, "replayed_unit_scale_to_m": 1.0e-3,
        "part_names": names, "replayed_part_names": names,
        "part_colors_rgb": colors, "replayed_part_colors_rgb": colors,
        "part_transforms_in_source_units": transforms, "replayed_part_transforms_in_source_units": transforms,
        "part_shape_ids": shape_ids, "replayed_part_shape_ids": shape_ids,
        "assembly_hierarchy": hierarchy, "replayed_assembly_hierarchy": hierarchy,
        "export_owner": "headless:step-assembly-725", "replayed_export_owner": "headless:step-assembly-725",
        "step_file_sha256": "6" * 64, "replayed_step_file_sha256": "6" * 64,
        "assembly_result_sha256": "7" * 64, "accepted_assembly_result_sha256": "7" * 64,
    }
    return row


_GEAR = "involute_gear_module_teeth_pressureangle_backlash_volume_inertia_brep_generation_identity"

_SHEET = "sheetmetal_bend_radius_kfactor_thickness_neutralaxis_volume_brep_generation_identity"

_LOCATION = "location_composition_local_global_rotation_order_subshape_owner_result_generation_identity"

_REBUILD = "parametric_rebuild_dependency_cache_property_invalidation_export_owner_result_generation_identity"


def _public_v43():
    reference, measured = _public_v42()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "involute-gear-726"
            row[_GEAR] = {
                "gear_generation": generation,
                **_generation(generation, "module_generation", "teeth_generation", "pressure_angle_generation", "backlash_generation", "volume_generation", "inertia_generation", "owner_generation", "brep_generation", "result_generation"),
                "module_m": 0.002, "result_module_m": 0.002,
                "tooth_count": 24, "result_tooth_count": 24,
                "pressure_angle_deg": 20.0, "result_pressure_angle_deg": 20.0,
                "backlash_m": 2.0e-5, "result_backlash_m": 2.0e-5,
                "pitch_diameter_m": 0.048, "result_pitch_diameter_m": 0.048,
                "base_diameter_m": 0.048 * __import__("math").cos(__import__("math").radians(20.0)),
                "result_base_diameter_m": 0.048 * __import__("math").cos(__import__("math").radians(20.0)),
                "gear_volume_m3": 2.4e-5, "result_gear_volume_m3": 2.4e-5,
                "inertia_tensor_kg_m2": [[1.0e-8, 0.0, 0.0], [0.0, 2.0e-8, 0.0], [0.0, 0.0, 3.0e-8]],
                "result_inertia_tensor_kg_m2": [[1.0e-8, 0.0, 0.0], [0.0, 2.0e-8, 0.0], [0.0, 0.0, 3.0e-8]],
                "shape_owner": "part:involute-gear-726", "result_shape_owner": "part:involute-gear-726",
                "gear_brep_sha256": suffix * 64, "accepted_gear_brep_sha256": suffix * 64,
            }
            generation = "sheetmetal-bend-726"
            row[_SHEET] = {
                "sheet_generation": generation,
                **_generation(generation, "bend_radius_generation", "kfactor_generation", "thickness_generation", "neutral_axis_generation", "volume_generation", "owner_generation", "brep_generation", "result_generation"),
                "bend_radius_m": 0.003, "result_bend_radius_m": 0.003,
                "bend_angle_deg": 90.0, "result_bend_angle_deg": 90.0,
                "k_factor": 0.42, "result_k_factor": 0.42,
                "thickness_m": 0.001, "result_thickness_m": 0.001,
                "neutral_axis_radius_m": 0.00342, "result_neutral_axis_radius_m": 0.00342,
                "volume_m3": 1.25e-5, "result_volume_m3": 1.25e-5,
                "surface_area_m2": 0.025, "result_surface_area_m2": 0.025,
                "shape_owner": "part:sheetmetal-726", "result_shape_owner": "part:sheetmetal-726",
                "sheet_brep_sha256": suffix * 64, "accepted_sheet_brep_sha256": suffix * 64,
            }
    return reference, measured


def _source_v43():
    row = _source_v42()
    identity = row["replay_identity"]
    generation = "location-compose-726"
    identity[_LOCATION] = {
        "location_generation": generation,
        **_generation(generation, "local_generation", "global_generation", "rotation_generation", "order_generation", "subshape_generation", "owner_generation", "result_generation"),
        "local_location": {"translation_m": [0.01, 0.0, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "global_location": {"translation_m": [0.01, 0.02, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "replayed_local_location": {"translation_m": [0.01, 0.0, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "replayed_global_location": {"translation_m": [0.01, 0.02, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "rotation_order": ["Z", "Y", "X"], "replayed_rotation_order": ["Z", "Y", "X"],
        "translation_frame": "parent", "replayed_translation_frame": "parent",
        "subshape_id": "face:31", "replayed_subshape_id": "face:31",
        "shape_owner": "headless:location-compose-726", "replayed_shape_owner": "headless:location-compose-726",
        "location_result_sha256": "8" * 64, "accepted_location_result_sha256": "8" * 64,
    }
    generation = "parametric-rebuild-726"
    identity[_REBUILD] = {
        "rebuild_generation": generation,
        **_generation(generation, "dependency_generation", "cache_generation", "property_generation", "invalidation_generation", "topology_generation", "export_generation", "owner_generation", "result_generation"),
        "dependency_values": {"length_m": 0.1, "radius_m": 0.02},
        "replayed_dependency_values": {"length_m": 0.1, "radius_m": 0.02},
        "dependency_identity_sha256": "04d499f658a26ae1536b0756e01599404c8e2694d14fe7149437857a10ad60d0",
        "replayed_dependency_identity_sha256": "04d499f658a26ae1536b0756e01599404c8e2694d14fe7149437857a10ad60d0",
        "cache_key": "length=0.1;radius=0.02", "replayed_cache_key": "length=0.1;radius=0.02",
        "invalidated_properties": ["volume", "surface_area", "center_of_mass"],
        "replayed_invalidated_properties": ["volume", "surface_area", "center_of_mass"],
        "topology_signature": {"solid": 1, "face": 8}, "replayed_topology_signature": {"solid": 1, "face": 8},
        "export_owner": "headless:parametric-rebuild-726", "replayed_export_owner": "headless:parametric-rebuild-726",
        "rebuild_result_sha256": "9" * 64, "accepted_rebuild_result_sha256": "9" * 64,
    }
    return row


_BOOLEAN_V44 = "boolean_shell_fillet_massproperties_centerofmass_volume_brep_owner_generation_identity"

_LOFT_V44 = "loft_sweep_section_orientation_tangent_area_volume_inertia_brep_generation_identity"

_SKETCH = "sketch_constraint_solver_order_plane_frame_parameter_cache_shape_owner_generation_identity"

_STEP_V44 = "step_export_units_tessellation_tolerance_facecount_brep_digest_owner_generation_identity"


def _public_v44():
    reference, measured = _public_v43()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = str(index + 1) * 64
            generation = "boolean-shell-fillet-v44-731"
            row[_BOOLEAN_V44] = {
                "boolean_generation": generation,
                **{name: generation for name in ("history_generation", "shell_generation", "fillet_generation", "mass_generation", "center_generation", "volume_generation", "owner_generation", "brep_generation", "result_generation")},
                "operation": "cut", "result_operation": "cut",
                "shell_thickness_m": 0.001, "result_shell_thickness_m": 0.001,
                "fillet_radius_m": 0.0005, "result_fillet_radius_m": 0.0005,
                "center_of_mass_m": [0.01, 0.02, 0.03], "result_center_of_mass_m": [0.01, 0.02, 0.03],
                "volume_m3": 9.9e-4, "result_volume_m3": 9.9e-4,
                "surface_area_m2": 6.1e-2, "result_surface_area_m2": 6.1e-2,
                "topology_signature": {"solid": 1, "shell": 1, "face": 10}, "result_topology_signature": {"solid": 1, "shell": 1, "face": 10},
                "shape_owner": "part:boolean-shell-fillet-731", "result_shape_owner": "part:boolean-shell-fillet-731",
                "boolean_brep_sha256": digest, "accepted_boolean_brep_sha256": digest,
            }
            generation = "loft-sweep-v44-731"
            row[_LOFT_V44] = {
                "loft_generation": generation,
                **{name: generation for name in ("section_generation", "orientation_generation", "tangent_generation", "area_generation", "volume_generation", "inertia_generation", "owner_generation", "brep_generation", "result_generation")},
                "section_count": 3, "result_section_count": 3,
                "section_orientation": "consistent_ccw", "result_section_orientation": "consistent_ccw",
                "tangent_continuity": "C1", "result_tangent_continuity": "C1",
                "section_area_m2": [1e-4, 1.2e-4, 1e-4], "result_section_area_m2": [1e-4, 1.2e-4, 1e-4],
                "volume_m3": 2.5e-5, "result_volume_m3": 2.5e-5,
                "inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]],
                "result_inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]],
                "shape_owner": "part:loft-sweep-731", "result_shape_owner": "part:loft-sweep-731",
                "loft_brep_sha256": digest, "accepted_loft_brep_sha256": digest,
            }
    return reference, measured


def _source_v44():
    row = _source_v43()
    identity = row["replay_identity"]
    generation = "sketch-replay-v44-731"
    identity[_SKETCH] = {
        "sketch_generation": generation,
        **{name: generation for name in ("constraint_generation", "solver_order_generation", "plane_frame_generation", "parameter_cache_generation", "shape_generation", "owner_generation", "result_generation")},
        "constraint_order": ["coincident:1", "distance:2", "horizontal:3"], "replayed_constraint_order": ["coincident:1", "distance:2", "horizontal:3"],
        "plane_frame": "XY", "replayed_plane_frame": "XY", "parameter_cache_key": "width=0.1;radius=0.02", "replayed_parameter_cache_key": "width=0.1;radius=0.02",
        "solver_status": "solved", "replayed_solver_status": "solved", "shape_generation_id": 731, "replayed_shape_generation_id": 731,
        "shape_owner": "headless:sketch-731", "replayed_shape_owner": "headless:sketch-731", "sketch_result_sha256": "a" * 64, "accepted_sketch_result_sha256": "a" * 64,
    }
    generation = "step-export-v44-731"
    identity[_STEP_V44] = {
        "step_generation": generation,
        **{name: generation for name in ("unit_generation", "tessellation_generation", "tolerance_generation", "facecount_generation", "brep_generation", "digest_generation", "owner_generation", "result_generation")},
        "length_unit": "mm", "replayed_length_unit": "mm", "tessellation_tolerance_m": 1e-5, "replayed_tessellation_tolerance_m": 1e-5,
        "face_count": 42, "replayed_face_count": 42, "topology_signature": {"solid": 1, "shell": 1, "face": 42}, "replayed_topology_signature": {"solid": 1, "shell": 1, "face": 42},
        "brep_generation_id": 731, "replayed_brep_generation_id": 731, "export_owner": "headless:step-export-731", "replayed_export_owner": "headless:step-export-731",
        "step_digest_sha256": "b" * 64, "replayed_step_digest_sha256": "b" * 64, "accepted_step_result_sha256": "c" * 64,
    }
    return row


def _public_v45():
    reference, measured = _public_v44()
    for rows in (reference, *measured.values()):
        for row in rows:
            row["boolean_fillet_shell_massproperties_volume_area_brep_owner_identity"] = {
                "generation": "boolean-shell-fillet-v45-812", "operation": "cut", "result_operation": "cut",
                "shell_thickness_m": 0.001, "result_shell_thickness_m": 0.001, "fillet_radius_m": 0.0005, "result_fillet_radius_m": 0.0005,
                "center_of_mass_m": [0.01, 0.02, 0.03], "result_center_of_mass_m": [0.01, 0.02, 0.03], "volume_m3": 9.9e-4, "result_volume_m3": 9.9e-4, "surface_area_m2": 6.1e-2, "result_surface_area_m2": 6.1e-2,
                "topology_signature": {"solid": 1, "shell": 1, "face": 10}, "result_topology_signature": {"solid": 1, "shell": 1, "face": 10}, "shape_owner": "part:boolean-shell-fillet-v45-812", "result_shape_owner": "part:boolean-shell-fillet-v45-812", "boolean_brep_sha256": "1" * 64, "accepted_boolean_brep_sha256": "1" * 64,
            }
            row["loft_sweep_section_frame_tangent_continuity_inertia_export_digest_identity"] = {
                "generation": "loft-sweep-v45-812", "section_count": 3, "result_section_count": 3, "section_orientation": "consistent_ccw", "result_section_orientation": "consistent_ccw", "tangent_continuity": "C1", "result_tangent_continuity": "C1", "section_area_m2": [1e-4, 1.2e-4, 1e-4], "result_section_area_m2": [1e-4, 1.2e-4, 1e-4], "volume_m3": 2.5e-5, "result_volume_m3": 2.5e-5, "inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]], "result_inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]], "shape_owner": "part:loft-sweep-v45-812", "result_shape_owner": "part:loft-sweep-v45-812", "loft_brep_sha256": "2" * 64, "accepted_loft_brep_sha256": "2" * 64,
            }
    return reference, measured


def _source_v45():
    row = _source_v44()
    replay = row["replay_identity"]
    replay["sketch_constraint_order_plane_frame_solver_cache_shape_generation_owner_identity"] = {"generation": "sketch-v45-812", "constraint_order": ["coincident:1", "distance:2"], "replayed_constraint_order": ["coincident:1", "distance:2"], "plane_frame": "XY", "replayed_plane_frame": "XY", "parameter_cache_key": "width=0.1", "replayed_parameter_cache_key": "width=0.1", "solver_status": "solved", "replayed_solver_status": "solved", "shape_generation_id": 812, "replayed_shape_generation_id": 812, "shape_owner": "headless:sketch-v45-812", "replayed_shape_owner": "headless:sketch-v45-812", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64}
    replay["step_units_tessellation_tolerance_face_topology_brep_export_owner_identity"] = {"generation": "step-v45-812", "length_unit": "mm", "replayed_length_unit": "mm", "tessellation_tolerance_m": 1e-5, "replayed_tessellation_tolerance_m": 1e-5, "face_topology": {"solid": 1, "shell": 1, "face": 42}, "replayed_face_topology": {"solid": 1, "shell": 1, "face": 42}, "brep_generation_id": 812, "replayed_brep_generation_id": 812, "export_owner": "headless:step-v45-812", "replayed_export_owner": "headless:step-v45-812", "step_digest_sha256": "4" * 64, "replayed_step_digest_sha256": "4" * 64, "accepted_result_sha256": "4" * 64}
    return row


PLACEMENT = "placement_unit_scale_inertia_coordinate_frame_nonfinite_identity"

BOOLEAN_V46 = "boolean_partial_shape_fillet_recovery_topology_identity"

SKETCH_V46 = "sketch_solver_partial_constraint_warning_plane_unit_identity"

STEP_V46 = "step_import_partial_face_tolerance_coordinate_frame_checksum_identity"


def _public_payload() -> dict[str, object]:
    placement = {
        "generation": "placement-unit-inertia-test",
        "placement_generation": "placement-unit-inertia-test",
        "result_placement_generation": "placement-unit-inertia-test",
        "unit_scale_generation": "placement-unit-inertia-test",
        "result_unit_scale_generation": "placement-unit-inertia-test",
        "inertia_frame_generation": "placement-unit-inertia-test",
        "result_inertia_frame_generation": "placement-unit-inertia-test",
        "mass_property_generation": "placement-unit-inertia-test",
        "result_mass_property_generation": "placement-unit-inertia-test",
        "placement_m": [0.01, 0.02, 0.03],
        "result_placement_m": [0.01, 0.02, 0.03],
        "unit_scale_to_si": 0.001,
        "result_unit_scale_to_si": 0.001,
        "inertia_frame": "global_cartesian",
        "result_inertia_frame": "global_cartesian",
        "finite_mass_property_status": "finite",
        "result_finite_mass_property_status": "finite",
        "shape_owner": "part:placement-unit-inertia-test",
        "accepted_shape_owner": "part:placement-unit-inertia-test",
        "result_sha256": "1" * 64,
        "accepted_result_sha256": "1" * 64,
    }
    boolean = {
        "generation": "boolean-recovery-topology-test",
        "boolean_generation": "boolean-recovery-topology-test",
        "result_boolean_generation": "boolean-recovery-topology-test",
        "fillet_generation": "boolean-recovery-topology-test",
        "result_fillet_generation": "boolean-recovery-topology-test",
        "recovery_generation": "boolean-recovery-topology-test",
        "result_recovery_generation": "boolean-recovery-topology-test",
        "topology_generation": "boolean-recovery-topology-test",
        "result_topology_generation": "boolean-recovery-topology-test",
        "boolean_status": "recovered",
        "result_boolean_status": "recovered",
        "topology_signature": {"solid": 1, "shell": 1, "face": 10},
        "result_topology_signature": {"solid": 1, "shell": 1, "face": 10},
        "partial_shape_status": "complete",
        "result_partial_shape_status": "complete",
        "shape_owner": "part:boolean-recovery-topology-test",
        "accepted_shape_owner": "part:boolean-recovery-topology-test",
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    row = {PLACEMENT: placement, BOOLEAN_V46: boolean}
    return {"reference": [row], "measured": {"headless": [deepcopy(row)]}}


def _source_payload() -> dict[str, object]:
    return {
        "replay_identity": {
            SKETCH_V46: {
                "generation": "sketch-partial-constraint-test",
                "result_generation": "sketch-partial-constraint-test",
                "solver_warning": "none",
                "result_solver_warning": "none",
                "constraint_state": "fully_constrained",
                "result_constraint_state": "fully_constrained",
                "plane": "XY",
                "result_plane": "XY",
                "unit_scale_to_si": 1.0,
                "result_unit_scale_to_si": 1.0,
                "result_sha256": "3" * 64,
                "accepted_result_sha256": "3" * 64,
            },
            STEP_V46: {
                "generation": "step-partial-face-test",
                "result_generation": "step-partial-face-test",
                "partial_face_count": 0,
                "result_partial_face_count": 0,
                "tolerance_m": 1.0e-6,
                "result_tolerance_m": 1.0e-6,
                "coordinate_frame": "global_cartesian",
                "result_coordinate_frame": "global_cartesian",
                "checksum_sha256": "4" * 64,
                "result_checksum_sha256": "4" * 64,
                "result_sha256": "5" * 64,
                "accepted_result_sha256": "5" * 64,
            },
        }
    }


def _payloads_v47() -> tuple[dict[str, object], dict[str, object]]:
    compound_generation = "compound-v47"
    assembly_generation = "assembly-v47"
    children = ["body:a", "body:b", "body:c"]
    masses = [1.0, 2.0, 3.0]
    hierarchy = ["root", "root/arm", "root/arm/tool"]
    owners = {name: "assembly:a1" for name in hierarchy}
    row = {
        COMPOUND: {"generation":compound_generation,"child_generation":compound_generation,"mass_property_generation":compound_generation,"result_generation":compound_generation,"child_ids":children,"result_child_ids":children,"child_masses_kg":masses,"result_child_masses_kg":masses,"aggregate_mass_kg":6.0,"result_aggregate_mass_kg":6.0,"owner":"compound:a1","accepted_owner":"compound:a1","result_sha256":"1"*64,"accepted_result_sha256":"1"*64},
        ASSEMBLY_V47: {"generation":assembly_generation,"mate_generation":assembly_generation,"hierarchy_generation":assembly_generation,"transform_generation":assembly_generation,"result_generation":assembly_generation,"mate_hierarchy":hierarchy,"result_mate_hierarchy":hierarchy,"transform_owner_map":owners,"result_transform_owner_map":owners,"local_to_world_transform_sha256":"2"*64,"result_local_to_world_transform_sha256":"2"*64,"owner":"assembly:a1","accepted_owner":"assembly:a1","result_sha256":"3"*64,"accepted_result_sha256":"3"*64},
    }
    public = {"reference":[row],"measured":{"cad":[deepcopy(row)]}}
    roundtrip_generation = "roundtrip-v47"
    external_generation = "external-v47"
    entities = [{"label":"body","uuid":"uuid-body","owner":"body:1"},{"label":"face-top","uuid":"uuid-top","owner":"body:1"}]
    references = ["edge:10@rev7", "vertex:2@rev7"]
    edges = [["sketch:1","edge:10"],["edge:10","geometry:rev7"]]
    source = {"replay_identity":{
        ROUNDTRIP:{"generation":roundtrip_generation,"roundtrip_generation":roundtrip_generation,"result_generation":roundtrip_generation,"entity_identities":entities,"result_entity_identities":entities,"duplicate_label_count":0,"result_duplicate_label_count":0,"duplicate_uuid_count":0,"result_duplicate_uuid_count":0,"result_sha256":"4"*64,"accepted_result_sha256":"4"*64},
        EXTERNAL:{"generation":external_generation,"reference_generation":external_generation,"dependency_generation":external_generation,"result_generation":external_generation,"geometry_revision":"rev7","result_geometry_revision":"rev7","external_references":references,"result_external_references":references,"dependency_edges":edges,"result_dependency_edges":edges,"dependency_cycle_count":0,"result_dependency_cycle_count":0,"result_sha256":"5"*64,"accepted_result_sha256":"5"*64},
    }}
    return public, source


def _payloads_v48() -> tuple[dict[str, object], dict[str, object]]:
    suppression_generation = "feature-config-mass-v48-901"
    loft_generation = "loft-correspondence-v48-901"
    configuration = {"length_mm": 40.0, "hole_enabled": False, "fillet_enabled": True}
    profiles = ["wire:base", "wire:mid", "wire:top"]
    orientation = {profile: 1 for profile in profiles}
    seams = {"wire:base": "vertex:1", "wire:mid": "vertex:5", "wire:top": "vertex:9"}
    correspondence = [["vertex:1", "vertex:5", "vertex:9"], ["vertex:2", "vertex:6", "vertex:10"]]
    row = {
        SUPPRESSION: {
            "generation": suppression_generation,
            "suppression_generation": suppression_generation,
            "configuration_generation": suppression_generation,
            "shape_generation": suppression_generation,
            "cache_generation": suppression_generation,
            "result_generation": suppression_generation,
            "suppressed_features": ["hole:1"],
            "result_suppressed_features": ["hole:1"],
            "configuration": configuration,
            "result_configuration": configuration,
            "shape_sha256": "6" * 64,
            "cache_shape_sha256": "6" * 64,
            "mass_kg": 1.25,
            "cached_mass_kg": 1.25,
            "result_mass_kg": 1.25,
            "owner": "build:feature-config-v48-901",
            "cache_owner": "build:feature-config-v48-901",
            "accepted_owner": "build:feature-config-v48-901",
            "result_sha256": "7" * 64,
            "accepted_result_sha256": "7" * 64,
        },
        LOFT_V48: {
            "generation": loft_generation,
            "profile_generation": loft_generation,
            "orientation_generation": loft_generation,
            "seam_generation": loft_generation,
            "correspondence_generation": loft_generation,
            "diagnostic_generation": loft_generation,
            "result_generation": loft_generation,
            "profile_order": profiles,
            "result_profile_order": profiles,
            "profile_orientation": orientation,
            "result_profile_orientation": orientation,
            "wire_seams": seams,
            "result_wire_seams": seams,
            "profile_correspondence": correspondence,
            "result_profile_correspondence": correspondence,
            "self_intersection_count": 0,
            "result_self_intersection_count": 0,
            "owner": "build:loft-v48-901",
            "accepted_owner": "build:loft-v48-901",
            "result_sha256": "8" * 64,
            "accepted_result_sha256": "8" * 64,
        },
    }
    import_generation = "import-metadata-v48-901"
    selector_generation = "selector-witness-v48-901"
    layers = {"body:1": "mechanical", "body:2": "insulation"}
    colors = {"body:1": [0.8, 0.2, 0.1], "body:2": [0.2, 0.4, 0.9]}
    mapping = {"source:body-1": "body:1", "source:body-2": "body:2"}
    entities = ["face:11", "face:12"]
    witnesses = {"face:11": [0.0, 0.0, 1.0], "face:12": [0.0, 0.0, -1.0]}
    source = {
        "replay_identity": {
            IMPORT: {
                "generation": import_generation,
                "unit_generation": import_generation,
                "metadata_generation": import_generation,
                "subshape_generation": import_generation,
                "result_generation": import_generation,
                "source_revision": "source-rev-v48-901",
                "result_source_revision": "source-rev-v48-901",
                "inferred_unit": "mm",
                "result_inferred_unit": "mm",
                "unit_scale_to_m": 0.001,
                "result_unit_scale_to_m": 0.001,
                "layer_map": layers,
                "result_layer_map": layers,
                "color_map": colors,
                "result_color_map": colors,
                "persistent_subshape_map": mapping,
                "result_persistent_subshape_map": mapping,
                "owner": "import:source-rev-v48-901",
                "accepted_owner": "import:source-rev-v48-901",
                "result_sha256": "9" * 64,
                "accepted_result_sha256": "9" * 64,
            },
            SELECTOR_V48: {
                "generation": selector_generation,
                "query_generation": selector_generation,
                "witness_generation": selector_generation,
                "history_generation": selector_generation,
                "result_generation": selector_generation,
                "query": "faces().filter_by(Axis.Z)",
                "result_query": "faces().filter_by(Axis.Z)",
                "expected_cardinality": 2,
                "result_cardinality": 2,
                "selected_entities": entities,
                "result_selected_entities": entities,
                "witness_points": witnesses,
                "result_witness_points": witnesses,
                "feature_history_owner": "history:feature-v48-901",
                "result_feature_history_owner": "history:feature-v48-901",
                "result_sha256": "a" * 64,
                "accepted_result_sha256": "a" * 64,
            },
        }
    }
    return {"reference": [row], "measured": {"cad": [deepcopy(row)]}}, source


def _payloads_v49() -> tuple[dict[str, object], dict[str, object]]:
    generation = "build123d-v49"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    row = {
        ASSEMBLY_V49: {
            "generation": generation, "density_generation": generation, "material_generation": generation,
            "occurrence_generation": generation, "transform_generation": generation, "suppression_generation": generation,
            "mass_generation": generation, "result_generation": generation,
            "material_by_occurrence": {"occurrence:a": "steel"}, "result_material_by_occurrence": {"occurrence:a": "steel"},
            "density_kg_m3_by_occurrence": {"occurrence:a": 7850.0}, "result_density_kg_m3_by_occurrence": {"occurrence:a": 7850.0},
            "occurrence_transforms": {"occurrence:a": [0.0, 0.0, 0.0, 1.0]}, "result_occurrence_transforms": {"occurrence:a": [0.0, 0.0, 0.0, 1.0]},
            "suppressed_occurrences": [], "result_suppressed_occurrences": [], "assembly_mass_kg": 12.5, "result_assembly_mass_kg": 12.5,
            "assembly_owner": "assembly:v49", "result_assembly_owner": "assembly:v49", **result,
        },
        SKETCH_V49: {
            "generation": generation, "constraint_generation": generation, "dof_generation": generation, "plane_generation": generation,
            "unit_generation": generation, "wire_generation": generation, "result_generation": generation,
            "constraints": ["horizontal:e1", "distance:e1=40mm"], "result_constraints": ["horizontal:e1", "distance:e1=40mm"],
            "remaining_dof": 0, "result_remaining_dof": 0, "work_plane": "Plane.XY", "result_work_plane": "Plane.XY",
            "length_unit": "mm", "result_length_unit": "mm", "profile_wires": ["wire:outer"], "result_profile_wires": ["wire:outer"],
            "sketch_owner": "sketch:v49", "result_sketch_owner": "sketch:v49", **result,
        },
    }
    step = {
        "generation": generation, "schema_generation": generation, "assembly_generation": generation, "metadata_generation": generation,
        "unit_generation": generation, "tolerance_generation": generation, "result_generation": generation,
        "ap_schema": "AP242", "result_ap_schema": "AP242", "assembly_structure": {"assembly:root": ["part:a"]},
        "result_assembly_structure": {"assembly:root": ["part:a"]}, "color_map": {"part:a": [0.8, 0.2, 0.1]},
        "result_color_map": {"part:a": [0.8, 0.2, 0.1]}, "layer_map": {"part:a": "structure"}, "result_layer_map": {"part:a": "structure"},
        "length_unit": "mm", "result_length_unit": "mm", "model_tolerance_m": 1e-7, "result_model_tolerance_m": 1e-7,
        "import_owner": "import:v49", "result_import_owner": "import:v49", **result,
    }
    boolean = {
        "generation": generation, "deletion_generation": generation, "selector_generation": generation, "adjacency_generation": generation,
        "mass_generation": generation, "cache_generation": generation, "result_generation": generation,
        "deleted_subshapes": ["face:old"], "result_deleted_subshapes": ["face:old"],
        "selector_results": {"selector:top": ["face:new"]}, "result_selector_results": {"selector:top": ["face:new"]},
        "adjacency_map": {"face:new": ["edge:new"]}, "result_adjacency_map": {"face:new": ["edge:new"]},
        "mass_kg": 3.75, "cached_mass_kg": 3.75, "result_mass_kg": 3.75,
        "cache_shape_sha256": "e" * 64, "result_cache_shape_sha256": "e" * 64,
        "history_owner": "history:v49", "result_history_owner": "history:v49", **result,
    }
    return {"reference": [row], "measured": {"cad": [deepcopy(row)]}}, {"replay_identity": {STEP_V49: step, BOOLEAN_V49: boolean}}


def _payloads_v50() -> tuple[dict[str, object], dict[str, object]]:
    generation = "build123d-v50"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    frames = {
        "base": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        "shaft": [0.0, 0.0, 0.02, 1.0, 0.0, 0.0, 0.0],
    }
    transforms = {
        "occurrence:base": frames["base"],
        "occurrence:shaft": frames["shaft"],
    }
    selectors = {"fillet": ["edge:11", "edge:12"], "chamfer": ["edge:21"]}
    history = {
        "edge:11": ["edge:31", "face:41"],
        "edge:12": ["edge:32", "face:42"],
        "edge:21": ["edge:33", "face:43"],
    }
    row = {
        MATE_V50: {
            "generation": generation,
            "constraint_generation": generation,
            "dof_generation": generation,
            "frame_generation": generation,
            "occurrence_generation": generation,
            "transform_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "mate_constraints": ["fixed:base", "revolute:shaft-base"],
            "result_mate_constraints": ["fixed:base", "revolute:shaft-base"],
            "remaining_dof": 1,
            "result_remaining_dof": 1,
            "mate_frames": frames,
            "result_mate_frames": frames,
            "occurrence_transforms": transforms,
            "result_occurrence_transforms": transforms,
            "assembly_owner": "assembly:mates-v50",
            "result_assembly_owner": "assembly:mates-v50",
            **result,
        },
        FEATURE: {
            "generation": generation,
            "selector_generation": generation,
            "radius_generation": generation,
            "topology_generation": generation,
            "history_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "edge_selectors": selectors,
            "result_edge_selectors": selectors,
            "fillet_radius_m": 0.002,
            "result_fillet_radius_m": 0.002,
            "chamfer_distance_m": 0.001,
            "result_chamfer_distance_m": 0.001,
            "topology_history": history,
            "result_topology_history": history,
            "shape_owner": "shape:features-v50",
            "result_shape_owner": "shape:features-v50",
            **result,
        },
    }
    orientations = {"shell:outer": 1, "shell:inner": -1}
    healing = {
        "generation": generation,
        "tolerance_generation": generation,
        "healing_generation": generation,
        "sewing_generation": generation,
        "shell_generation": generation,
        "solid_generation": generation,
        "orientation_generation": generation,
        "owner_generation": generation,
        "result_generation": generation,
        "input_tolerance_m": 1e-7,
        "result_input_tolerance_m": 1e-7,
        "healing_applied": True,
        "result_healing_applied": True,
        "sewing_tolerance_m": 5e-7,
        "result_sewing_tolerance_m": 5e-7,
        "shell_count": 2,
        "result_shell_count": 2,
        "solid_count": 1,
        "result_solid_count": 1,
        "shell_orientation_signs": orientations,
        "result_shell_orientation_signs": orientations,
        "shape_owner": "shape:healed-v50",
        "result_shape_owner": "shape:healed-v50",
        **result,
    }
    normals = {"outward": 12480, "inward": 0, "degenerate": 0}
    stl = {
        "generation": generation,
        "linear_generation": generation,
        "angular_generation": generation,
        "triangle_generation": generation,
        "normal_generation": generation,
        "unit_generation": generation,
        "owner_generation": generation,
        "result_generation": generation,
        "linear_deflection_m": 1e-4,
        "result_linear_deflection_m": 1e-4,
        "angular_deflection_rad": 0.1,
        "result_angular_deflection_rad": 0.1,
        "triangle_count": 12480,
        "result_triangle_count": 12480,
        "normal_counts": normals,
        "result_normal_counts": normals,
        "length_unit": "m",
        "result_length_unit": "m",
        "mesh_owner": "mesh:stl-v50",
        "result_mesh_owner": "mesh:stl-v50",
        **result,
    }
    public = {"reference": [row], "measured": {"cad": [deepcopy(row)]}}
    source = {"replay_identity": {HEALING: healing, STL: stl}}
    return public, source


def _records() -> tuple[dict[str, object], dict[str, object]]:
    mass = "mass-v51"
    sweep = "sweep-v51"
    step = "step-v51"
    brep = "brep-v51"
    inertia = [[0.10, 0.0, 0.0], [0.0, 0.20, 0.0], [0.0, 0.0, 0.25]]
    shifted = [[0.10, 0.0, 0.0], [0.0, 0.22, 0.0], [0.0, 0.0, 0.27]]
    history = {"profile:1": ["face:11"], "path:1": ["edge:21"], "result": ["solid:31"]}
    public_row = {
        "mass_properties_frame_inertia_parallel_axis_density_shape_owner_identity": {"generation": mass, **{name: mass for name in ("frame_generation", "inertia_generation", "shift_generation", "density_generation", "owner_generation", "result_generation")}, "coordinate_frame": "global_cartesian", "result_coordinate_frame": "global_cartesian", "mass_kg": 2.0, "result_mass_kg": 2.0, "centroidal_inertia_kg_m2": inertia, "result_centroidal_inertia_kg_m2": inertia, "parallel_axis_shift_m": [0.1, 0.0, 0.0], "result_parallel_axis_shift_m": [0.1, 0.0, 0.0], "shifted_inertia_kg_m2": shifted, "result_shifted_inertia_kg_m2": shifted, "density_kg_m3": 7800.0, "result_density_kg_m3": 7800.0, "shape_owner": "shape:mass-v51", "result_shape_owner": "shape:mass-v51", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64},
        "sweep_profile_path_trihedron_transition_selfintersection_history_owner_identity": {"generation": sweep, **{name: sweep for name in ("profile_generation", "path_generation", "trihedron_generation", "transition_generation", "intersection_generation", "history_generation", "owner_generation", "result_generation")}, "profile_id": "profile:1", "result_profile_id": "profile:1", "path_id": "path:1", "result_path_id": "path:1", "trihedron": "corrected_frenet", "result_trihedron": "corrected_frenet", "transition": "transformed", "result_transition": "transformed", "self_intersections": [], "result_self_intersections": [], "topology_history": history, "result_topology_history": history, "shape_owner": "shape:sweep-v51", "result_shape_owner": "shape:sweep-v51", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
    }
    references = {"document:main": ["part:shaft", "part:housing"]}
    names = {"occurrence:1": "shaft", "occurrence:2": "housing"}
    colors = {"occurrence:1": [0.8, 0.8, 0.8], "occurrence:2": [0.2, 0.4, 0.8]}
    location = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    source = {"replay_identity": {
        "step_external_reference_occurrence_name_schema_unit_color_owner_identity": {"generation": step, **{name: step for name in ("reference_generation", "name_generation", "schema_generation", "unit_generation", "color_generation", "owner_generation", "result_generation")}, "external_references": references, "result_external_references": references, "occurrence_names": names, "result_occurrence_names": names, "step_schema": "AP242", "result_step_schema": "AP242", "length_unit": "m", "result_length_unit": "m", "occurrence_colors": colors, "result_occurrence_colors": colors, "document_owner": "document:step-v51", "result_document_owner": "document:step-v51", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64},
        "brep_occt_version_location_precision_triangulation_cache_owner_identity": {"generation": brep, **{name: brep for name in ("version_generation", "location_generation", "precision_generation", "triangulation_generation", "owner_generation", "result_generation")}, "occt_version": "7.9.0", "result_occt_version": "7.9.0", "shape_location": location, "result_shape_location": location, "model_precision_m": 1e-7, "result_model_precision_m": 1e-7, "triangulation_cache_sha256": "4" * 64, "result_triangulation_cache_sha256": "4" * 64, "shape_owner": "shape:brep-v51", "result_shape_owner": "shape:brep-v51", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64},
    }}
    return {"reference": [public_row], "measured": {"cad": [deepcopy(public_row)]}}, source


def _generation_v52(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": prefix, **{name: prefix for name in names}}


def _payloads_v52():
    selected = ["face:1", "face:2"]; labels = {"face:1": "a", "face:2": "b"}
    selector = {**_generation_v52("sel-v52", ("query_generation", "order_generation", "topology_generation", "label_generation", "owner_generation", "result_generation")), "selector_query": "faces().sort_by(Axis.X)", "result_selector_query": "faces().sort_by(Axis.X)", "selected_topology_order": selected, "result_selected_topology_order": selected, "topology_revision": "topology:v52", "result_topology_revision": "topology:v52", "topology_labels": labels, "result_topology_labels": labels, "shape_owner": "shape:sel-v52", "result_shape_owner": "shape:sel-v52", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64}
    frame = {"origin": [0.0, 0.0, 0.0], "x_dir": [1.0, 0.0, 0.0], "z_dir": [0.0, 0.0, 1.0]}; edges = ["edge:1", "edge:2", "edge:3"]
    workplane = {**_generation_v52("wp-v52", ("coordinate_generation", "edge_generation", "closure_generation", "owner_generation", "result_generation")), "local_frame": frame, "result_local_frame": frame, "pending_edge_order": edges, "result_pending_edge_order": edges, "wire_closed": True, "result_wire_closed": True, "builder_owner": "builder:wp-v52", "result_builder_owner": "builder:wp-v52", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64}
    public = {"reference": [{SELECTOR_V52: selector, WORKPLANE: workplane}], "measured": {}}
    location = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    brep = {**_generation_v52("brep-v52", ("version_generation", "location_generation", "tshape_generation", "serialization_generation", "owner_generation", "result_generation")), "occt_version": "7.9.0", "replayed_occt_version": "7.9.0", "shape_location": location, "replayed_shape_location": location, "tshape_sha256": "c" * 64, "replayed_tshape_sha256": "c" * 64, "serialization_sha256": "d" * 64, "replayed_serialization_sha256": "d" * 64, "shape_owner": "shape:brep-v52", "replayed_shape_owner": "shape:brep-v52", "result_sha256": "e" * 64, "accepted_result_sha256": "e" * 64}
    materials = {"part:a": "material:steel"}; instances = {"instance:a": location}
    gltf = {**_generation_v52("gltf-v52", ("axis_generation", "scale_generation", "material_generation", "instance_generation", "owner_generation", "result_generation")), "axis_convention": "Y_up_right_handed", "replayed_axis_convention": "Y_up_right_handed", "length_scale_to_m": 1.0, "replayed_length_scale_to_m": 1.0, "part_materials": materials, "replayed_part_materials": materials, "instance_transforms": instances, "replayed_instance_transforms": instances, "scene_owner": "scene:gltf-v52", "replayed_scene_owner": "scene:gltf-v52", "result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    return public, {"replay_identity": {BREP_V52: brep, GLTF: gltf}}


def _generations_v53(generation: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{name: generation for name in names}}


def _payloads_v53():
    frame = {"origin": [0.0, 0.0, 0.0], "x_dir": [1.0, 0.0, 0.0], "z_dir": [0.0, 0.0, 1.0]}
    pair = ["component:a", "component:b"]
    mate = {**_generations_v53("mate-v53", ("frame_generation", "handedness_generation", "axis_generation", "offset_generation", "component_generation", "owner_generation", "result_generation")), "mate_frame": frame, "result_mate_frame": frame, "handedness": "right", "result_handedness": "right", "mate_axis": [0.0, 0.0, 1.0], "result_mate_axis": [0.0, 0.0, 1.0], "offset_m": -0.001, "result_offset_m": -0.001, "component_pair": pair, "result_component_pair": pair, "assembly_owner": "assembly:v53", "result_assembly_owner": "assembly:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    colors = {"component:a": [0.1, 0.2, 0.3], "component:b": [0.7, 0.8, 0.9]}; hierarchy = {"assembly:root": pair}
    step = {**_generations_v53("step-v53", ("unit_generation", "tolerance_generation", "color_generation", "hierarchy_generation", "owner_generation", "result_generation")), "length_unit": "mm", "result_length_unit": "mm", "linear_tolerance_m": 1.0e-6, "result_linear_tolerance_m": 1.0e-6, "component_colors": colors, "result_component_colors": colors, "component_hierarchy": hierarchy, "result_component_hierarchy": hierarchy, "document_owner": "document:v53", "result_document_owner": "document:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    ancestry = {"face:r1": ["face:a1", "face:b1"]}
    history = {**_generations_v53("history-v53", ("boolean_generation", "ancestry_generation", "fillet_generation", "chamfer_generation", "owner_generation", "result_generation")), "boolean_operation": "cut", "replayed_boolean_operation": "cut", "face_ancestry": ancestry, "replayed_face_ancestry": ancestry, "fillet_edges": ["edge:1"], "replayed_fillet_edges": ["edge:1"], "chamfer_edges": ["edge:2"], "replayed_chamfer_edges": ["edge:2"], "shape_owner": "shape:v53", "replayed_shape_owner": "shape:v53", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64}
    tensor = [[0.02, 0.001, 0.0], [0.001, 0.03, 0.0], [0.0, 0.0, 0.04]]
    mass = {**_generations_v53("mass-v53", ("tensor_generation", "frame_generation", "centroid_generation", "density_generation", "owner_generation", "result_generation")), "inertia_tensor_kg_m2": tensor, "replayed_inertia_tensor_kg_m2": tensor, "reference_frame": frame, "replayed_reference_frame": frame, "centroid_m": [0.1, 0.0, 0.0], "replayed_centroid_m": [0.1, 0.0, 0.0], "density_kg_m3": 7850.0, "replayed_density_kg_m3": 7850.0, "solid_owner": "solid:v53", "replayed_solid_owner": "solid:v53", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64}
    return {"reference": [{MATE_V53: mate, STEP_V53: step}], "measured": {}}, {"replay_identity": {HISTORY: history, MASS: mass}}


def _generations_v54(generation: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{name: generation for name in names}}


def _payloads_v54():
    locations = {"solid:a": {"translation_m": [0.0, 0.0, 0.0], "quaternion_wxyz": [1.0, 0.0, 0.0, 0.0]}}
    densities = {"solid:a": 7850.0}
    tensor = [[0.02, 0.001, 0.0], [0.001, 0.03, 0.0], [0.0, 0.0, 0.04]]
    assembly = {**_generations_v54("assembly-v54", ("density_generation", "location_generation", "center_generation", "inertia_generation", "owner_generation", "result_generation")), "solid_densities_kg_m3": densities, "result_solid_densities_kg_m3": densities, "located_solids": locations, "result_located_solids": locations, "center_of_mass_m": [0.0, 0.0, 0.0], "result_center_of_mass_m": [0.0, 0.0, 0.0], "inertia_tensor_kg_m2": tensor, "result_inertia_tensor_kg_m2": tensor, "assembly_owner": "assembly:v54", "result_assembly_owner": "assembly:v54", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64}
    sections = [{"section": "wire:1", "orientation": 1, "parameter_start": 0.0, "seam": "vertex:1"}, {"section": "wire:2", "orientation": 1, "parameter_start": 0.0, "seam": "vertex:2"}]
    topology = {"solids": 1, "shells": 1, "faces": 6, "edges": 12, "vertices": 8}
    loft = {**_generations_v54("loft-v54", ("section_generation", "parameter_generation", "seam_generation", "topology_generation", "owner_generation", "result_generation")), "section_correspondence": sections, "result_section_correspondence": sections, "resulting_topology": topology, "result_resulting_topology": topology, "shape_owner": "shape:loft-v54", "result_shape_owner": "shape:loft-v54", "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64}
    structure = {"product:root": ["part:a"], "part:a": ["solid:a"]}; colors = {"part:a": [0.2, 0.3, 0.4]}; layers = {"part:a": "layer:main"}
    step = {**_generations_v54("step-v54", ("unit_generation", "structure_generation", "color_generation", "layer_generation", "owner_generation", "result_generation")), "schema": "AP242", "replayed_schema": "AP242", "length_unit": "mm", "replayed_length_unit": "mm", "product_structure": structure, "replayed_product_structure": structure, "component_colors": colors, "replayed_component_colors": colors, "component_layers": layers, "replayed_component_layers": layers, "document_owner": "document:v54", "replayed_document_owner": "document:v54", "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64}
    vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]; triangles = [[0, 1, 2]]; orientations = [1]
    tessellation = {**_generations_v54("tess-v54", ("deflection_generation", "angle_generation", "orientation_generation", "index_generation", "owner_generation", "result_generation")), "linear_deflection_m": 1.0e-4, "replayed_linear_deflection_m": 1.0e-4, "angular_deflection_rad": 0.2, "replayed_angular_deflection_rad": 0.2, "vertices_m": vertices, "replayed_vertices_m": vertices, "triangle_indices": triangles, "replayed_triangle_indices": triangles, "triangle_orientations": orientations, "replayed_triangle_orientations": orientations, "shape_owner": "shape:tess-v54", "replayed_shape_owner": "shape:tess-v54", "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64}
    return {"reference": [{ASSEMBLY_V54: assembly, LOFT_V54: loft}], "measured": {}}, {"replay_identity": {STEP_V54: step, TESSELLATION: tessellation}}


def _payloads_v55():
    generation = "build123d-v55-test"
    generations = lambda names: {name: generation for name in names}
    topology = {"solids": 1, "shells": 1, "faces": 12, "edges": 26, "vertices": 16}
    thread = {"generation": generation, **generations(("pitch_generation", "handedness_generation", "start_generation", "topology_generation", "volume_generation", "owner_generation", "result_generation")), "pitch_m": 2.0e-3, "result_pitch_m": 2.0e-3, "handedness": "right", "result_handedness": "right", "start_count": 1, "result_start_count": 1, "thread_topology": topology, "result_thread_topology": topology, "volume_m3": 1.25e-5, "result_volume_m3": 1.25e-5, "shape_owner": "shape:thread-v55", "result_shape_owner": "shape:thread-v55", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    thickness = 1.0e-3; radius = 5.0e-3; angle = math.pi / 2.0; k = 0.4; allowance = angle * (radius + k * thickness)
    sheet = {"generation": generation, **generations(("thickness_generation", "bend_generation", "neutralaxis_generation", "pattern_generation", "area_generation", "owner_generation", "result_generation")), "thickness_m": thickness, "result_thickness_m": thickness, "inside_bend_radius_m": radius, "result_inside_bend_radius_m": radius, "bend_angle_rad": angle, "result_bend_angle_rad": angle, "neutral_axis_factor": k, "result_neutral_axis_factor": k, "bend_allowance_m": allowance, "result_bend_allowance_m": allowance, "flat_pattern_area_m2": 0.020, "result_flat_pattern_area_m2": 0.020, "folded_surface_area_m2": 0.020, "result_folded_surface_area_m2": 0.020, "shape_owner": "shape:sheet-v55", "result_shape_owner": "shape:sheet-v55", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    tolerances = [{"feature": "feature:hole-1", "kind": "position", "value": 0.05, "datum_refs": ["A", "B", "C"]}]
    pmi = {"generation": generation, **generations(("tolerance_generation", "datum_generation", "unit_generation", "product_generation", "revision_generation", "owner_generation", "result_generation")), "geometric_tolerances": tolerances, "replayed_geometric_tolerances": tolerances, "datum_frame": {"A": "face:base", "B": "face:side", "C": "axis:hole"}, "replayed_datum_frame": {"A": "face:base", "B": "face:side", "C": "axis:hole"}, "length_unit": "mm", "replayed_length_unit": "mm", "product_association": {"feature:hole-1": "part:housing"}, "replayed_product_association": {"feature:hole-1": "part:housing"}, "document_revision": "document:v55-r3", "replayed_document_revision": "document:v55-r3", "document_owner": "document:pmi-v55", "replayed_document_owner": "document:pmi-v55", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64}
    shells = [{"shell": "shell:1", "faces": ["face:1", "face:2", "face:3", "face:4"], "closed": True}]; orientations = {"face:1": 1, "face:2": 1, "face:3": 1, "face:4": 1}
    brep = {"generation": generation, **generations(("tolerance_generation", "shell_generation", "orientation_generation", "volume_generation", "owner_generation", "result_generation")), "healing_tolerance_m": 1.0e-6, "replayed_healing_tolerance_m": 1.0e-6, "sewn_shells": shells, "replayed_sewn_shells": shells, "face_orientations": orientations, "replayed_face_orientations": orientations, "closed_volume_count": 1, "replayed_closed_volume_count": 1, "volume_m3": 8.0e-5, "replayed_volume_m3": 8.0e-5, "shape_owner": "shape:brep-v55", "replayed_shape_owner": "shape:brep-v55", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64}
    return {"reference": [{THREAD: thread, SHEET: sheet}], "measured": {}}, {"replay_identity": {PMI: pmi, BREP_V55: brep}}


def _payloads_v56() -> tuple[dict[str, object], dict[str, object]]:
    generation = "build123d-v56-test"
    generations = lambda names: {name: generation for name in names}
    module = 2.0e-3
    teeth = 24
    pressure = math.radians(20.0)
    gear = {
        "generation": generation,
        **generations(("module_generation", "tooth_generation", "pressure_generation", "pitch_generation", "volume_generation", "owner_generation", "result_generation")),
        "module_m": module, "result_module_m": module,
        "tooth_count": teeth, "result_tooth_count": teeth,
        "pressure_angle_rad": pressure, "result_pressure_angle_rad": pressure,
        "pitch_diameter_m": module * teeth, "result_pitch_diameter_m": module * teeth,
        "volume_m3": 1.62e-5, "result_volume_m3": 1.62e-5,
        "shape_owner": "shape:gear-v56", "result_shape_owner": "shape:gear-v56",
        "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
    }
    frame = {"method": "parallel_transport", "samples": 65, "closed": False}
    pipe = {
        "generation": generation,
        **generations(("path_generation", "frame_generation", "twist_generation", "intersection_generation", "volume_generation", "owner_generation", "result_generation")),
        "path_length_m": 0.125, "result_path_length_m": 0.125,
        "transport_frame": frame, "result_transport_frame": frame,
        "net_twist_rad": 0.0, "result_net_twist_rad": 0.0,
        "self_intersection": False, "result_self_intersection": False,
        "volume_m3": 2.45e-5, "result_volume_m3": 2.45e-5,
        "shape_owner": "shape:pipe-v56", "result_shape_owner": "shape:pipe-v56",
        "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
    }
    transform = [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03], [0.0, 0.0, 0.0, 1.0]]
    occurrence = {"id": "occurrence:gear-1", "product": "product:gear"}
    assembly_frame = {"parent": "assembly:root", "child": "occurrence:gear-1"}
    step = {
        "generation": generation,
        **generations(("transform_generation", "unit_generation", "product_generation", "frame_generation", "owner_generation", "result_generation")),
        "occurrence_transform_4x4": transform, "replayed_occurrence_transform_4x4": transform,
        "length_unit": "m", "replayed_length_unit": "m",
        "product_occurrence": occurrence, "replayed_product_occurrence": occurrence,
        "assembly_frame": assembly_frame, "replayed_assembly_frame": assembly_frame,
        "document_owner": "document:step-v56", "replayed_document_owner": "document:step-v56",
        "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64,
    }
    mesh = {
        "generation": generation,
        **generations(("format_generation", "watertight_generation", "manifold_generation", "unit_generation", "volume_generation", "owner_generation", "result_generation")),
        "mesh_format": "3mf", "replayed_mesh_format": "3mf",
        "watertight": True, "replayed_watertight": True,
        "manifold": True, "replayed_manifold": True,
        "length_unit": "mm", "replayed_length_unit": "mm",
        "unit_scale_to_m": 1.0e-3, "replayed_unit_scale_to_m": 1.0e-3,
        "signed_volume_m3": 4.8e-5, "replayed_signed_volume_m3": 4.8e-5,
        "mesh_owner": "mesh:3mf-v56", "replayed_mesh_owner": "mesh:3mf-v56",
        "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
    }
    return {"reference": [{GEAR: gear, PIPE: pipe}], "measured": {}}, {"replay_identity": {STEP_V56: step, MESH: mesh}}
