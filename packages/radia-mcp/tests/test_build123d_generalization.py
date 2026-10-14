from __future__ import annotations

import json
from copy import deepcopy

import pytest

from radia_mcp.build123d.assembly_exchange_identity_v53 import (
    HISTORY,
    MASS,
    MATE as MATE_V53,
    STEP as STEP_V53,
    validate_public_identity as validate_public_identity_v53,
    validate_source_identity as validate_source_identity_v53,
)
from radia_mcp.build123d.assembly_replay_identity_v49 import (
    ASSEMBLY as ASSEMBLY_V49,
    BOOLEAN as BOOLEAN_V49,
    SKETCH as SKETCH_V49,
    STEP as STEP_V49,
    validate_public_identity as validate_public_identity_v49,
    validate_source_identity as validate_source_identity_v49,
)
from radia_mcp.build123d.assembly_tessellation_identity_v54 import (
    ASSEMBLY as ASSEMBLY_V54,
    LOFT as LOFT_V54,
    STEP as STEP_V54,
    TESSELLATION,
    validate_public_identity as validate_public_identity_v54,
    validate_source_identity as validate_source_identity_v54,
)
from radia_mcp.build123d.build123d_v46_identity import (
    validate_public_identity as validate_public_identity_v46,
    validate_source_identity as validate_source_identity_v46,
)
from radia_mcp.build123d.cross_artifact_cad_lineage_v47 import (
    ASSEMBLY as ASSEMBLY_V47,
    COMPOUND,
    EXTERNAL,
    ROUNDTRIP,
    validate_public_identity as validate_public_identity_v47,
    validate_source_identity as validate_source_identity_v47,
)
from radia_mcp.build123d.feature_replay_identity_v50 import (
    FEATURE,
    HEALING,
    MATE as MATE_V50,
    STL,
    validate_public_identity as validate_public_identity_v50,
    validate_source_identity as validate_source_identity_v50,
)
from radia_mcp.build123d.gear_pipe_exchange_identity_v56 import (
    GEAR,
    MESH,
    PIPE,
    STEP as STEP_V56,
    validate_public_identity as validate_public_identity_v56,
    validate_source_identity as validate_source_identity_v56,
)
from radia_mcp.build123d.mass_sweep_exchange_identity_v51 import (
    validate_public_identity as validate_public_identity_v51,
    validate_source_identity as validate_source_identity_v51,
)
from radia_mcp.build123d.selector_exchange_identity_v52 import (
    GLTF,
    SELECTOR as SELECTOR_V52,
    WORKPLANE,
    validate_public_identity as validate_public_identity_v52,
    validate_source_identity as validate_source_identity_v52,
)
from radia_mcp.build123d.semantic_cad_identity_v48 import (
    IMPORT,
    LOFT as LOFT_V48,
    SELECTOR as SELECTOR_V48,
    SUPPRESSION,
    validate_public_identity as validate_public_identity_v48,
    validate_source_identity as validate_source_identity_v48,
)
from radia_mcp.build123d.server import (
    build123d_jointed_assembly_source_replay_gate,
)
from radia_mcp.build123d.thread_sheet_identity_v55 import (
    BREP as BREP_V55,
    PMI,
    SHEET,
    THREAD,
    validate_public_identity as validate_public_identity_v55,
    validate_source_identity as validate_source_identity_v55,
)
from _build123d_generalization_payloads import (
    BOOLEAN_V46,
    PLACEMENT,
    SKETCH_V46,
    STEP_V46,
    _ASSEMBLY_V39,
    _ASSEMBLY_V41,
    _BOOLEAN_V42,
    _BOOLEAN_V44,
    _DRAFT,
    _GEAR,
    _HEAL,
    _HELIX,
    _LOCATION,
    _LOFT_V39,
    _LOFT_V41,
    _LOFT_V44,
    _OBB,
    _OFFSET,
    _REBUILD,
    _SELECTOR,
    _SHEET,
    _SHELL,
    _SKETCH,
    _STEP_V41,
    _STEP_V42,
    _STEP_V44,
    _TESSELLATION,
    _THREAD,
    _payloads_v47,
    _payloads_v48,
    _payloads_v49,
    _payloads_v50,
    _payloads_v52,
    _payloads_v53,
    _payloads_v54,
    _payloads_v55,
    _payloads_v56,
    _public,
    _public_payload,
    _public_result,
    _public_v11,
    _public_v12,
    _public_v13,
    _public_v14,
    _public_v15,
    _public_v16,
    _public_v17,
    _public_v18,
    _public_v19,
    _public_v20,
    _public_v21,
    _public_v22,
    _public_v23,
    _public_v24,
    _public_v25,
    _public_v26,
    _public_v27,
    _public_v28,
    _public_v29,
    _public_v30,
    _public_v31,
    _public_v32,
    _public_v33,
    _public_v34,
    _public_v35,
    _public_v36,
    _public_v37,
    _public_v38,
    _public_v39,
    _public_v40,
    _public_v41,
    _public_v42,
    _public_v43,
    _public_v44,
    _public_v45,
    _public_v7,
    _records,
    _source,
    _source_payload,
    _source_result,
    _source_v11,
    _source_v12,
    _source_v13,
    _source_v14,
    _source_v15,
    _source_v16,
    _source_v17,
    _source_v18,
    _source_v19,
    _source_v20,
    _source_v21,
    _source_v22,
    _source_v23,
    _source_v24,
    _source_v25,
    _source_v26,
    _source_v27,
    _source_v28,
    _source_v29,
    _source_v30,
    _source_v31,
    _source_v32,
    _source_v33,
    _source_v34,
    _source_v35,
    _source_v36,
    _source_v37,
    _source_v38,
    _source_v39,
    _source_v40,
    _source_v41,
    _source_v42,
    _source_v43,
    _source_v44,
    _source_v45,
    _source_v7,
)


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_consistent_unit_scale_error",
        "v7_public_missing_body_mass_compensation",
    ],
)
def test_generalization_v7_public(case_id):
    reference, measured = _public_v7()
    rows = measured["external_cad"]
    if case_id == "v7_public_consistent_unit_scale_error":
        for row in rows:
            row["volume"] *= 1.0e9
            row["area"] *= 1.0e6
            for key in ("min", "max", "center", "size"):
                row["bounding_box"][key] = [value * 1.0e3 for value in row["bounding_box"][key]]
            row["bounding_box"]["diagonal"] *= 1.0e3
    else:
        removed = rows.pop()
        rows[0]["volume"] += removed["volume"]
        rows[0]["area"] += removed["area"]
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"


@pytest.mark.parametrize("identity", [[], {"cad_artifacts": ["not-a-row"]}])
def test_source_identity_contract_rejects_malformed_evidence(identity):
    row = _source_v7()
    row["replay_identity"] = identity
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"


@pytest.mark.parametrize(
    ("case_id", "failed_check"),
    [
        (
            "v7_source_stale_commit_fresh_cad_digest",
            "neutral_cad_artifacts_bind_current_source_commit",
        ),
        (
            "v7_source_external_kernel_version_drift",
            "external_kernel_versions_are_replay_invariant",
        ),
    ],
)
def test_generalization_v7_source(case_id, failed_check):
    row = _source_v7()
    if case_id == "v7_source_stale_commit_fresh_cad_digest":
        row["replay_identity"]["cad_artifacts"][0]["source_commit"] = "d" * 40
    else:
        row["replay_identity"]["external_kernel"]["replay_versions"][1] = "7.9.0"
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][failed_check] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_public_mass_properties_older_than_brep",
        "v8_public_assembly_child_revision_mix",
    ],
)
def test_generalization_v8_public(case_id: str) -> None:
    reference, measured = _public()
    rows = measured["external_cad"]
    if case_id == "v8_public_mass_properties_older_than_brep":
        rows[0]["mass_property_identity"].update(
            {"brep_revision": "brep-frame-41", "brep_sha256": "f" * 64}
        )
        expected = "mass_properties_bind_current_brep_revision"
    else:
        rows[1]["assembly_identity"]["child_revisions"][
            "insert"
        ] = "child-insert-2"
        expected = "assembly_children_match_reference_revision_map"
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_source_export_precedes_source_replay",
        "v8_source_kernel_session_restarts_between_replays",
    ],
)
def test_generalization_v8_source(case_id: str) -> None:
    row = _source()
    identity = row["replay_identity"]
    if case_id == "v8_source_export_precedes_source_replay":
        identity["cad_artifacts"][0][
            "export_completed_utc"
        ] = "2026-07-16T01:59:59Z"
        expected = "neutral_cad_export_follows_source_replay"
    else:
        identity["external_kernel"]["replay_sessions"][1].update(
            {
                "session_generation": "occt-session-43",
                "process_start_utc": "2026-07-16T02:00:30Z",
            }
        )
        expected = "external_kernel_session_generation_is_continuous"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_public_center_of_mass_reference_frame_mismatch",
        "v9_public_face_adjacency_before_fillet",
    ],
)
def test_generalization_v9_public(case_id: str) -> None:
    reference, measured = _public()
    rows = measured["external_cad"]
    if case_id == "v9_public_center_of_mass_reference_frame_mismatch":
        rows[0]["mass_property_frame_identity"].update(
            {
                "frame_id": "component-local-frame",
                "transform_generation": "component-transform-9",
            }
        )
        expected = "mass_property_centers_share_reference_frames"
    else:
        rows[0]["topology_identity"].update(
            {
                "brep_revision": "brep-frame-41",
                "face_adjacency_sha256": "4" * 64,
            }
        )
        expected = "face_adjacency_matches_current_brep_revision"
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_source_import_heal_topology_digest_omitted",
        "v9_source_external_volume_before_unit_conversion",
    ],
)
def test_generalization_v9_source(case_id: str) -> None:
    row = _source()
    identity = row["replay_identity"]
    if case_id == "v9_source_import_heal_topology_digest_omitted":
        identity["topology_replay_identity"]["imports"][1].pop("topology_sha256")
        expected = "heal_and_noheal_imports_record_topology_identity"
    else:
        identity["unit_conversion_identity"].update(
            {
                "external_measurement_stage": "before_unit_conversion",
                "external_volume_unit": "mm^3",
                "declared_volume_scale_to_target": 1.0e-9,
            }
        )
        expected = "external_volume_is_measured_after_unit_conversion"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v10_public_compound_overlap_double_count_volume",
        "v10_public_center_of_mass_frame_volume_frame_mismatch",
    ],
)
def test_generalization_v10_public(case_id: str) -> None:
    reference, measured = _public()
    row = measured["external_cad"][0]
    if case_id == "v10_public_compound_overlap_double_count_volume":
        row["compound_volume_identity"].update(
            {"reported_volume_basis": "child_volume_sum", "overlap_volume": 0.125}
        )
        expected = "compound_volume_uses_physical_union_not_child_sum"
    else:
        row["placement_transform_identity"][
            "center_of_mass_transform_generation"
        ] = "assembly-transform-41"
        expected = "center_of_mass_uses_final_placement_transform"
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v10_source_boolean_result_digest_before_shape_clean",
        "v10_source_tessellation_tolerance_previous_shape",
    ],
)
def test_generalization_v10_source(case_id: str) -> None:
    row = _source()
    identity = row["replay_identity"]
    if case_id == "v10_source_boolean_result_digest_before_shape_clean":
        identity["boolean_clean_identity"].update(
            {
                "shape_clean_input_sha256": "3" * 64,
                "export_shape_generation": "shape-generation-43",
            }
        )
        expected = "boolean_export_follows_shape_clean_identity"
    else:
        identity["tessellation_identity"].update(
            {
                "tolerance_shape_generation": "shape-generation-41",
                "tessellation_generation": "tessellation-generation-41",
            }
        )
        expected = "tessellation_tolerances_belong_to_current_shape"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    ("case_id", "expected"),
    [
        (
            "v11_public_shape_healed_after_mass_properties_stale",
            "mass_properties_follow_final_healed_brep",
        ),
        (
            "v11_public_mirrored_part_inertia_tensor_frame_mismatch",
            "mirrored_inertia_tensor_uses_final_global_frame",
        ),
    ],
)
def test_generalization_v11_public(case_id: str, expected: str) -> None:
    reference, measured = _public_v11()
    row = measured["external_cad"][0]
    if case_id == "v11_public_shape_healed_after_mass_properties_stale":
        row["shape_healing_identity"].update(
            {
                "mass_property_brep_sha256": "4" * 64,
                "mass_property_shape_generation": "shape-generation-42",
            }
        )
    else:
        row["inertia_tensor_identity"].update(
            {
                "tensor_frame_id": "source-part-frame",
                "tensor_transform_generation": "source-transform-42",
            }
        )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    ("case_id", "expected"),
    [
        (
            "v11_source_step_export_tolerance_previous_kernel_session",
            "step_export_tolerances_belong_to_current_kernel_and_shape",
        ),
        (
            "v11_source_assembly_instance_uuid_reused_after_replace",
            "replacement_rotates_instance_uuid_and_rebinds_placement",
        ),
    ],
)
def test_generalization_v11_source(case_id: str, expected: str) -> None:
    row = _source_v11()
    identity = row["replay_identity"]
    if case_id == "v11_source_step_export_tolerance_previous_kernel_session":
        identity["step_export_tolerance_identity"].update(
            {
                "tolerance_kernel_session_generation": "occt-session-41",
                "tolerance_shape_generation": "shape-generation-41",
            }
        )
    else:
        component = identity["assembly_replacement_identity"]["components"][0]
        component.update(
            {
                "current_instance_uuid": component["removed_instance_uuid"],
                "placement_shape_sha256": component["removed_shape_sha256"],
            }
        )
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


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


def test_v13_public_tessellation_tolerance_length_unit_mismatch():
    reference, measured = _public_v13()
    measured["external_cad"][0]["tessellation_tolerance_unit_identity"].update({"area_evaluation_deflection_unit": "m", "area_evaluation_deflection_scale_to_m": 1.0})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["tessellated_area_uses_one_length_unit_tolerance"] is False


def test_v13_public_compound_label_topology_index_previous_boolean():
    reference, measured = _public_v13()
    measured["external_cad"][0]["compound_label_topology_identity"].update({"label_table_boolean_generation": "boolean-generation-44", "selected_subshape_parent_sha256": "4" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["compound_labels_resolve_on_final_boolean_topology"] is False


def test_v13_source_step_geometry_unit_scale_metadata_mismatch():
    row = _source_v13()
    row["replay_identity"]["step_geometry_unit_scale_identity"].update({"metadata_generation": "step-import-44", "metadata_length_unit": "mm", "metadata_scale_to_m": 1.0e-3})
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["step_geometry_and_metadata_share_one_length_unit"] is False


def test_v13_source_selector_cache_previous_shape_generation():
    row = _source_v13()
    row["replay_identity"]["selector_cache_shape_identity"].update({"selector_cache_shape_generation": "shape-generation-44", "cached_selector_query_sha256": "5" * 64})
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["selector_cache_belongs_to_active_shape_generation"] is False


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


def test_v15_public_boolean_tolerance_model_length_unit_mismatch():
    reference, measured = _public_v15()
    measured["external_cad"][0]["boolean_tolerance_length_unit_identity"].update(
        {
            "tolerance_unit": "m",
            "tolerance_scale_to_m": 1.0,
            "kernel_tolerance_m": 1.0e-6,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["boolean_tolerance_uses_one_physical_model_length_basis"] is False


def test_v15_public_nested_assembly_placement_multiplication_order_mismatch():
    reference, measured = _public_v15()
    measured["external_cad"][0]["nested_assembly_placement_order_identity"].update(
        {
            "applied_multiplication_order": "child_then_parent",
            "world_placement_generation": "assembly-46",
            "world_transform_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["nested_assembly_placements_use_parent_then_child_order"] is False


def test_v15_source_step_color_label_map_previous_topology_generation():
    row = _source_v15()
    row["replay_identity"]["step_color_label_topology_identity"].update(
        {
            "attribute_map_topology_generation": "step-topology-46",
            "attribute_face_ids": [30, 32, 33],
            "resolved_attribute_map_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["step_color_labels_follow_current_topology_generation"] is False


def test_v15_source_brep_surface_parameter_orientation_range_mismatch():
    row = _source_v15()
    row["replay_identity"]["brep_surface_parameter_orientation_identity"].update(
        {
            "surface_parameter_generation": "brep-serialization-46",
            "exported_parameter_orientation": "v_cross_u_outward",
            "exported_parameter_range_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "brep_surface_parameter_ranges_preserve_outward_orientation"
        ]
        is False
    )


def test_v16_public_mass_inertia_tensor_reference_frame_placement_generation_mismatch():
    reference, measured = _public_v16()
    measured["external_cad"][0][
        "mass_inertia_reference_frame_placement_identity"
    ].update(
        {
            "mass_property_placement_generation": "placement-49",
            "inertia_reference_frame": "local",
            "mass_property_shape_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["mass_inertia_uses_final_world_placement_frame"] is False


def test_v16_public_loft_wire_correspondence_seam_normalization_generation_mismatch():
    reference, measured = _public_v16()
    measured["external_cad"][0]["loft_wire_correspondence_seam_identity"].update(
        {
            "wire_correspondence_seam_generation": "seam-49",
            "loft_section_wire_ids": [11, 13, 12],
            "loft_wire_correspondence_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "loft_sections_use_current_seam_normalized_correspondence"
        ]
        is False
    )


def test_v16_source_step_assembly_child_parent_length_unit_metadata_mismatch():
    row = _source_v16()
    row["replay_identity"]["step_assembly_child_parent_unit_identity"].update(
        {
            "child_placement_import_generation": "step-import-49",
            "child_placement_length_unit": "m",
            "child_placement_scale_to_m": 1.0,
            "resolved_assembly_placement_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "step_assembly_child_parent_placements_share_length_unit"
        ]
        is False
    )


def test_v16_source_selector_normal_predicate_preplacement_frame_mismatch():
    row = _source_v16()
    row["replay_identity"]["selector_normal_world_frame_identity"].update(
        {
            "selector_placement_generation": "selector-placement-49",
            "evaluated_normal_frame": "local",
            "resolved_face_ids": [20, 21],
            "evaluated_normal_table_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["selector_normals_use_final_world_placement_frame"] is False


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


def test_v18_public_nested_assembly_location_transform_composition_order_mismatch():
    reference, measured = _public_v18()
    measured["external_cad"][0][
        "nested_assembly_location_transform_composition_identity"
    ].update(
        {
            "location_transform_order_generation": "transform-order-51",
            "resolved_composed_transform_sha256": ["4" * 64, "3" * 64],
            "resolved_composition_order": "child_then_parent",
            "resolved_location_tree_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "nested_assembly_locations_use_current_transform_composition_order"
    ] is False


def test_v18_public_boolean_retained_face_name_history_refine_generation_mismatch():
    reference, measured = _public_v18()
    measured["external_cad"][0][
        "boolean_retained_face_name_history_refine_identity"
    ].update(
        {
            "retained_name_refine_generation": "refine-51",
            "topology_history_refine_generation": "refine-51",
            "resolved_face_names": ["outlet", "inlet"],
            "resolved_face_ids": [42, 41],
            "resolved_topology_history_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "boolean_retained_face_names_follow_post_refine_topology_history"
    ] is False


def test_v18_source_step_occurrence_color_material_inheritance_generation_mismatch():
    row = _source_v18()
    row["replay_identity"][
        "step_occurrence_color_material_inheritance_identity"
    ].update(
        {
            "occurrence_metadata_assembly_generation": "assembly-51",
            "color_inheritance_assembly_generation": "assembly-51",
            "material_inheritance_assembly_generation": "assembly-51",
            "metadata_parent_occurrence_ids": ["root/old", "root/old/frame"],
            "imported_colors_rgb": [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
            "imported_material_names": ["steel", "aluminum"],
            "imported_occurrence_metadata_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "step_occurrence_inheritance_uses_current_assembly_generation"
    ] is False


def test_v18_source_stl_chordal_angular_tolerance_model_unit_generation_mismatch():
    row = _source_v18()
    row["replay_identity"][
        "stl_tolerance_model_length_unit_generation_identity"
    ].update(
        {
            "tolerance_conversion_model_unit_generation": "model-unit-51",
            "chordal_tolerance_length_unit": "m",
            "tessellator_chordal_tolerance_si_m": 0.01,
            "angular_tolerance_unit": "rad",
            "tessellator_angular_tolerance_rad": 5.0,
            "tessellator_tolerance_contract_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "stl_tolerances_use_current_model_length_unit_generation"
    ] is False


def test_v19_public_boolean_history_subshape_label_fillet_reorder_mismatch():
    reference, measured = _public_v19()
    measured["external_cad"][0][
        "boolean_history_subshape_label_fillet_order_identity"
    ].update(
        {
            "label_fillet_generation": "fillet-52",
            "edge_order_fillet_generation": "fillet-52",
            "resolved_subshape_labels": ["outlet", "inlet"],
            "resolved_face_ids": [52, 51],
            "resolved_fillet_edge_ids": [73, 71, 72],
            "resolved_subshape_history_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "boolean_subshape_labels_follow_current_fillet_edge_order"
    ] is False


def test_v19_public_assembly_mate_frame_unit_location_generation_mismatch():
    reference, measured = _public_v19()
    measured["external_cad"][0][
        "assembly_mate_frame_unit_location_generation_identity"
    ].update(
        {
            "mate_frame_assembly_generation": "assembly-52",
            "unit_assembly_generation": "assembly-52",
            "mate_frame_length_unit": "m",
            "resolved_local_frame_sha256": ["2" * 64, "1" * 64],
            "resolved_parent_location_sha256": ["4" * 64, "3" * 64],
            "resolved_mate_resolution_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "assembly_mates_share_current_frame_unit_and_parent_location"
    ] is False


def test_v19_source_step_import_tolerance_unit_healing_generation_mismatch():
    row = _source_v19()
    row["replay_identity"][
        "step_import_tolerance_unit_healing_generation_identity"
    ].update(
        {
            "tolerance_import_generation": "step-import-52",
            "healed_edge_map_healing_generation": "healing-52",
            "tolerance_length_unit": "m",
            "tolerance_si_m": 1.0e-5,
            "imported_healed_edge_ids": [103, 102, 101],
            "imported_healed_edge_map_sha256": "e" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "step_import_tolerances_and_healed_edges_share_current_generation"
    ] is False


def test_v19_source_tessellation_vertex_index_normal_transform_generation_mismatch():
    row = _source_v19()
    row["replay_identity"][
        "tessellation_vertex_index_normal_transform_generation_identity"
    ].update(
        {
            "vertex_location_transform_generation": "location-52",
            "normal_location_transform_generation": "location-52",
            "transformed_triangle_indices": [[0, 2, 1], [0, 3, 2]],
            "rendered_tessellation_sha256": "e" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "tessellation_vertices_indices_and_normals_use_final_transform"
    ] is False


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


def test_v21_public_boolean_result_solid_orientation_location_label_generation_mismatch():
    reference, measured = _public_v21()
    measured["external_cad"][0][
        "boolean_result_solid_orientation_location_label_generation_identity"
    ].update(
        {
            "orientation_boolean_generation": "boolean-70",
            "location_boolean_generation": "boolean-69",
            "label_boolean_generation": "boolean-68",
            "result_operand_generations": ["operand-a-70", "operand-b-70"],
            "resolved_solid_orientation": "reversed",
            "resolved_semantic_labels": ["old_body", "old_interface"],
            "resolved_result_location_sha256": "b" * 64,
            "resolved_boolean_result_sha256": "c" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_result_uses_current_orientation_location_and_labels"
    ]


def test_v21_public_tessellation_chord_angle_unit_location_generation_mismatch():
    reference, measured = _public_v21()
    measured["external_cad"][0][
        "tessellation_chord_angle_unit_location_generation_identity"
    ].update(
        {
            "tessellation_shape_generation": "shape-70",
            "metric_tessellation_generation": "tessellation-70",
            "location_tessellation_generation": "tessellation-69",
            "evaluated_chord_tolerance": 0.1,
            "evaluated_angular_tolerance_deg": 25.0,
            "evaluated_length_unit": "mm",
            "evaluated_object_location_sha256": "d" * 64,
            "evaluated_tessellation_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "tessellation_uses_current_tolerances_units_and_object_location"
    ]


def test_v21_source_step_assembly_product_id_color_location_generation_mismatch():
    row = _source_v21()
    row["replay_identity"][
        "step_assembly_product_id_color_location_generation_identity"
    ].update(
        {
            "decoder_step_export_generation": "step-export-70",
            "color_assembly_generation": "assembly-70",
            "hierarchy_assembly_generation": "assembly-69",
            "location_assembly_generation": "assembly-68",
            "decoded_product_ids": ["root", "housing", "rotor"],
            "decoded_parent_product_ids": ["", "root", "housing"],
            "decoded_component_location_sha256": ["5" * 64, "7" * 64, "6" * 64],
            "decoded_assembly_metadata_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_assembly_uses_current_product_colors_hierarchy_and_locations"
    ]


def test_v21_source_sketch_constraint_entity_id_solver_order_generation_mismatch():
    row = _source_v21()
    row["replay_identity"][
        "sketch_constraint_entity_id_solver_order_generation_identity"
    ].update(
        {
            "entity_table_sketch_generation": "sketch-70",
            "solver_order_sketch_generation": "sketch-69",
            "replay_entity_ids": [101, 104, 103],
            "replay_solver_constraint_order": [202, 201],
            "replay_constraint_entity_ids": [[101, 104], [104, 103]],
            "replay_entity_table_sha256": "b" * 64,
            "replay_constraint_table_sha256": "c" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sketch_constraints_use_current_entity_ids_and_solver_order"
    ]


def test_v22_public_joint_connector_frame_labeled_face_subshape_generation_mismatch():
    reference, measured = _public_v22()
    measured["external_cad"][0][
        "joint_connector_frame_labeled_face_subshape_generation_identity"
    ].update(
        {
            "label_table_shape_generation": "shape-80",
            "connector_shape_generation": "shape-79",
            "location_shape_generation": "shape-78",
            "resolved_labeled_face_subshape_id": "face:stale_mount",
            "resolved_labeled_face_geometry_sha256": "8" * 64,
            "evaluated_connector_axis": [0.0, 1.0, 0.0],
            "evaluated_parent_location_sha256": "9" * 64,
            "evaluated_connector_frame_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "joint_connectors_use_current_labeled_face_frame_and_parent_location"
    ]


def test_v22_public_inertia_tensor_principal_axes_density_unit_location_generation_mismatch():
    reference, measured = _public_v22()
    measured["external_cad"][0][
        "inertia_tensor_principal_axes_density_unit_location_generation_identity"
    ].update(
        {
            "density_shape_generation": "shape-80",
            "mass_property_shape_generation": "shape-79",
            "location_shape_generation": "shape-78",
            "principal_axis_shape_generation": "shape-77",
            "evaluated_density_value": 7.8,
            "evaluated_density_unit": "g/cm^3",
            "evaluated_shape_location_sha256": "b" * 64,
            "evaluated_center_of_mass": [10.0, 20.0, 30.0],
            "evaluated_inertia_tensor_sha256": "c" * 64,
            "evaluated_principal_axes_sha256": "d" * 64,
            "evaluated_mass_property_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "inertia_uses_current_density_unit_location_and_principal_axes"
    ]


def test_v22_source_step_ap242_component_transform_name_product_generation_mismatch():
    row = _source_v22()
    row["replay_identity"][
        "step_ap242_component_transform_name_product_generation_identity"
    ].update(
        {
            "decoder_step_export_generation": "step-export-80",
            "product_assembly_generation": "assembly-80",
            "name_assembly_generation": "assembly-79",
            "transform_assembly_generation": "assembly-78",
            "decoded_component_product_ids": ["root", "housing", "rotor"],
            "decoded_component_names": ["machine", "old_housing", "old_rotor"],
            "decoded_nested_transform_sha256": ["1" * 64, "3" * 64, "2" * 64],
            "decoded_ap242_product_map_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_ap242_components_use_current_products_names_and_nested_transforms"
    ]


def test_v22_source_curved_mesh_export_edge_chord_surface_label_generation_mismatch():
    row = _source_v22()
    row["replay_identity"][
        "curved_mesh_export_edge_chord_surface_label_generation_identity"
    ].update(
        {
            "mesh_shape_generation": "shape-80",
            "edge_curve_shape_generation": "shape-79",
            "surface_label_shape_generation": "shape-78",
            "metric_mesh_export_generation": "mesh-export-80",
            "exported_edge_curve_sha256": ["8" * 64, "9" * 64],
            "evaluated_chordal_tolerance": 0.1,
            "evaluated_length_unit": "mm",
            "exported_boundary_surface_labels": ["stale_outer", "stale_interface"],
            "exported_curved_mesh_sha256": "a" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "curved_mesh_export_uses_current_edges_chord_and_surface_labels"
    ]


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


def test_v26_public_fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_mismatch():
    reference, measured = _public_v26()
    measured["external_cad"][0]["fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_identity"].update(
        {"selector_feature_generation": "fillet-130", "result_feature_type": "chamfer",
         "result_edge_selector_names": ["outer-top-2", "stale-edge"], "result_persistent_edge_ids": [102, 999],
         "result_linear_tolerance_m": 1.0e-3, "result_feature_shape_sha256": "b" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fillet_chamfer_features_use_current_selectors_topology_names_tolerance_and_shape"]


def test_v26_public_mass_density_center_inertia_reference_frame_assembly_generation_mismatch():
    reference, measured = _public_v26()
    measured["external_cad"][0]["mass_density_center_inertia_reference_frame_assembly_generation_identity"].update(
        {"density_mass_generation": "mass-130", "result_density_kg_m3": 2700.0,
         "result_mass_kg": 5.4, "result_center_of_mass_m": [0.03, 0.02, 0.01],
         "result_reference_frame": "part-local", "result_mass_property_sha256": "d" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["mass_properties_use_current_density_center_inertia_frame_and_assembly"]


def test_v26_source_stl_chord_tolerance_triangle_normal_orientation_component_digest_mismatch():
    row = _source_v26()
    row["replay_identity"]["stl_chord_tolerance_triangle_normal_orientation_component_digest_generation_identity"].update(
        {"tolerance_stl_generation": "stl-130", "decoded_linear_deflection_m": 0.001,
         "decoded_triangle_count": 120, "decoded_normal_orientation": "mixed",
         "decoded_component_ids": ["base", "cover", "rotor-old"], "decoded_stl_sha256": "f" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["stl_handoff_uses_current_tolerances_triangles_normals_components_and_digests"]


def test_v26_source_builder_context_workplane_local_frame_part_identity_cache_generation_mismatch():
    row = _source_v26()
    row["replay_identity"]["builder_context_workplane_local_frame_part_identity_cache_generation_identity"].update(
        {"stack_context_generation": "builder-130", "result_context_stack": ["BuildSketch", "BuildPart", "Locations"],
         "result_workplane_normal": [0.0, 1.0, 0.0], "result_part_ids": ["base", "cover", "rotor-old"],
         "result_builder_cache_sha256": "0" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["builder_replay_uses_current_context_workplane_frame_parts_and_cache"]


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


def test_v28_public_shell_offset_face_normal_thickness_join_self_intersection_mass_property_mismatch():
    reference, measured = _public_v28()
    measured["external_cad"][0]["shell_offset_face_normal_thickness_join_self_intersection_mass_generation_identity"].update(
        {"face_shell_generation": "shell-offset-150", "result_selected_face_names": ["side-face"],
         "result_normal_direction": "inward", "result_thickness_m": 0.003,
         "result_join_mode": "intersection", "result_self_intersection": True,
         "result_shape_generation": "shell-shape-150", "result_is_valid_solid": False,
         "result_volume_m3": 0.00037, "result_mass_property_sha256": "9" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["shell_offsets_use_current_faces_normals_thickness_join_intersection_shape_and_mass"]


def test_v28_public_path_sweep_frame_transition_profile_orientation_solid_validity_volume_mismatch():
    reference, measured = _public_v28()
    measured["external_cad"][0]["path_sweep_frame_transition_profile_orientation_solid_volume_generation_identity"].update(
        {"path_sweep_generation": "path-sweep-150", "result_path_edge_names": ["path-0", "path-2", "path-1"],
         "result_moving_frame": "frenet", "result_transition_mode": "right",
         "result_profile_wire_names": ["profile-inner"], "result_profile_orientation_deg": [0.0, -12.0, -24.0],
         "result_is_valid_solid": False, "result_volume_m3": 0.00025,
         "result_sweep_shape_sha256": "a" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["path_sweeps_use_current_path_frame_transition_profile_orientation_solid_and_volume"]


def test_v28_source_step_ap242_representation_context_product_uuid_unit_color_placement_mismatch():
    row = _source_v28()
    row["replay_identity"]["step_ap242_context_product_uuid_unit_color_placement_shape_file_generation_identity"].update(
        {"context_step_generation": "ap242-import-150", "decoded_representation_context": "geometric_curve_set",
         "decoded_product_uuid_map": [["base", "33333333-3333-4333-8333-333333333333"]],
         "decoded_length_unit": "m", "decoded_shape_owner_map": [["base-shape", "cover"]],
         "decoded_step_sha256": "b" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["step_ap242_import_uses_current_context_products_units_colors_placements_shapes_and_file"]


def test_v28_source_occt_kernel_shape_hash_location_tolerance_triangulation_serialization_cache_mismatch():
    row = _source_v28()
    row["replay_identity"]["occt_kernel_shape_location_tolerance_triangulation_serialization_cache_generation_identity"].update(
        {"kernel_cache_generation": "occt-cache-150", "decoded_kernel_version": "OCCT-7.7.0",
         "decoded_shape_sha256": "d" * 64, "decoded_modeling_tolerance_m": 1.0e-4,
         "decoded_triangulation_parameters": [0.01, 0.5, 0.0],
         "decoded_serialization_format": "step", "decoded_cache_sha256": "0" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["occt_shape_cache_uses_current_kernel_shape_location_tolerance_triangulation_and_serialization"]


def test_v29_public_sheet_metal_bend_allowance_kfactor_neutral_axis_flat_pattern_area_mismatch():
    reference, measured = _public_v29()
    measured["external_cad"][0][
        "sheet_metal_bend_allowance_kfactor_neutral_axis_relief_thickness_flat_pattern_area_generation_identity"
    ].update(
        {
            "bend_sheet_generation": "sheet-metal-160",
            "pattern_sheet_generation": "sheet-metal-159",
            "result_sheet_generation": "sheet-metal-158",
            "result_bend_radius_m": 0.004,
            "result_bend_angle_deg": 80.0,
            "result_k_factor": 0.5,
            "result_neutral_axis_radius_m": 0.00475,
            "result_bend_allowance_m": 0.0066,
            "result_relief_type": "none",
            "result_relief_width_m": 0.0,
            "result_thickness_m": 0.002,
            "result_strip_width_m": 0.06,
            "result_straight_lengths_m": [0.1, 0.07],
            "result_flat_pattern_area_m2": 0.0101,
            "result_flat_pattern_wire_closed": False,
            "result_flat_pattern_sha256": "d" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sheet_metal_flat_patterns_use_current_bends_neutral_axis_relief_thickness_and_area"
    ]


def test_v29_public_joint_kinematic_loop_dof_limit_frame_closure_swept_volume_mismatch():
    reference, measured = _public_v29()
    measured["external_cad"][0][
        "joint_kinematic_loop_graph_dof_limit_connector_frame_closure_configuration_swept_volume_generation_identity"
    ].update(
        {
            "graph_joint_generation": "joint-loop-160",
            "closure_joint_generation": "joint-loop-159",
            "result_joint_generation": "joint-loop-158",
            "result_joint_graph_edges": [["ground", "revolute-a", "slider-b"]],
            "result_dof_names": ["travel_b_m", "theta_a_deg"],
            "result_dof_types": ["prismatic", "fixed"],
            "result_lower_limits": [0.02, 30.0],
            "result_upper_limits": [0.0, -30.0],
            "result_configuration_values": [0.03, 45.0],
            "result_connector_frame_sha256": "e" * 64,
            "result_loop_closure_error_m": 1.0e-3,
            "result_loop_closure_tolerance_m": 1.0e-6,
            "result_configuration_id": "stale-pose",
            "result_swept_volume_m3": 0.00051,
            "result_swept_shape_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "joint_loops_use_current_graph_dofs_limits_frames_closure_configuration_and_swept_volume"
    ]


def test_v29_source_dxf_arc_spline_layer_plane_unit_closed_wire_face_digest_mismatch():
    row = _source_v29()
    row["replay_identity"][
        "dxf_arc_spline_layer_plane_unit_closed_wire_orientation_face_digest_generation_identity"
    ].update(
        {
            "arc_dxf_generation": "dxf-face-160",
            "face_dxf_generation": "dxf-face-159",
            "result_dxf_generation": "dxf-face-158",
            "decoded_arc_parameters": [["arc-1", 0.0, 0.0, 9.0, 0.0, 170.0]],
            "decoded_spline_parameters": [["spline-1", 2, 3, "0" * 64]],
            "decoded_entity_layer_map": [["arc-1", "construction"]],
            "decoded_workplane_matrix": [
                [1.0, 0.0, 0.0, 10.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            "decoded_length_unit": "in",
            "decoded_wire_closed": False,
            "decoded_wire_orientation": "clockwise",
            "decoded_face_area_mm2": 280.0,
            "decoded_dxf_sha256": "1" * 64,
            "decoded_face_sha256": "2" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "dxf_faces_use_current_arcs_splines_layers_plane_units_closure_orientation_and_digests"
    ]


def test_v29_source_3mf_component_transform_triangle_winding_material_watertight_volume_digest_mismatch():
    row = _source_v29()
    row["replay_identity"][
        "three_mf_component_transform_triangle_winding_material_watertight_volume_unit_digest_generation_identity"
    ].update(
        {
            "component_three_mf_generation": "3mf-import-160",
            "winding_three_mf_generation": "3mf-import-159",
            "result_three_mf_generation": "3mf-import-158",
            "decoded_component_names": ["insert", "housing"],
            "decoded_component_transform_sha256": [["housing", "3" * 64]],
            "decoded_triangle_winding": "inward_clockwise",
            "decoded_material_id_map": [["housing", 2], ["insert", 1]],
            "decoded_watertight_component_names": ["housing"],
            "decoded_component_volumes_mm3": [["housing", 1000.0], ["insert", -20.0]],
            "decoded_total_volume_mm3": 980.0,
            "decoded_length_unit": "m",
            "decoded_three_mf_sha256": "4" * 64,
            "decoded_triangle_mesh_sha256": "5" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "three_mf_imports_use_current_components_transforms_winding_materials_watertight_volumes_units_and_digests"
    ]


def test_v30_public_helical_sweep_pitch_handedness_profile_frame_self_intersection_volume_mismatch():
    reference, measured = _public_v30()
    measured["external_cad"][0][
        "helical_sweep_pitch_handedness_profile_frame_turn_self_intersection_volume_centroid_shape_generation_identity"
    ].update(
        {
            "pitch_helical_generation": "helical-sweep-170",
            "profile_helical_generation": "helical-sweep-169",
            "result_helical_generation": "helical-sweep-168",
            "result_pitch_m": 0.02,
            "result_handedness": "left",
            "result_profile_frame_matrix": [[1.0, 0.0, 0.0, 0.0]],
            "result_profile_frame_sha256": "f" * 64,
            "result_turn_count": 3.5,
            "result_self_intersection": True,
            "result_volume_m3": 0.8e-6,
            "result_centroid_m": [0.01, 0.0, 0.01],
            "result_helical_shape_sha256": "0" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "helical_sweeps_use_current_pitch_handedness_profile_frame_turns_intersection_mass_and_shape"
    ]


def test_v30_public_boolean_tolerance_operand_order_volume_centroid_inertia_history_mismatch():
    reference, measured = _public_v30()
    measured["external_cad"][0][
        "boolean_tolerance_operand_order_history_volume_centroid_inertia_shape_generation_identity"
    ].update(
        {
            "tolerance_boolean_generation": "boolean-history-170",
            "history_boolean_generation": "boolean-history-169",
            "result_boolean_generation": "boolean-history-168",
            "result_operation": "fuse",
            "result_model_tolerance_m": 1.0e-3,
            "result_operand_order": ["tool-b", "body-a"],
            "result_subshape_history_sha256": "1" * 64,
            "result_volume_m3": 2.4e-6,
            "result_centroid_m": [0.02, 0.0, 0.0],
            "result_inertia_tensor_kg_m2": [[3.0e-9, 1.0e-9, 0.0]],
            "result_boolean_shape_sha256": "2" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_results_use_current_operation_tolerance_operands_history_mass_inertia_and_shape"
    ]


def test_v30_source_step_assembly_product_name_color_transform_unit_instance_digest_mismatch():
    row = _source_v30()
    row["replay_identity"][
        "step_assembly_product_color_instance_transform_unit_hierarchy_file_generation_identity"
    ].update(
        {
            "product_step_generation": "step-assembly-170",
            "transform_step_generation": "step-assembly-169",
            "result_step_generation": "step-assembly-168",
            "decoded_product_names": ["assembly", "shaft", "housing-old"],
            "decoded_product_colors_rgb": [["housing", 0.8, 0.4, 0.2]],
            "decoded_instance_order": ["shaft-1", "housing-1"],
            "decoded_instance_transform_sha256": [["housing-1", "3" * 64]],
            "decoded_length_unit": "m",
            "decoded_assembly_hierarchy": [["housing-1", "assembly"]],
            "decoded_step_sha256": "4" * 64,
            "decoded_assembly_shape_sha256": "5" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_assemblies_use_current_products_colors_instances_transforms_units_hierarchy_and_digests"
    ]


def test_v30_source_stl_facet_normal_winding_watertight_tolerance_volume_digest_mismatch():
    row = _source_v30()
    row["replay_identity"][
        "stl_facet_normal_winding_watertight_tolerance_volume_unit_file_generation_identity"
    ].update(
        {
            "normal_stl_generation": "stl-solid-170",
            "watertight_stl_generation": "stl-solid-169",
            "result_stl_generation": "stl-solid-168",
            "decoded_facet_normals": [[0.0, 0.0, -1.0]],
            "decoded_triangle_winding": "inward_clockwise",
            "decoded_unmatched_edge_count": 4,
            "decoded_merge_tolerance_m": 1.0e-3,
            "decoded_signed_volume_m3": -1.2e-6,
            "decoded_length_unit": "mm",
            "decoded_stl_sha256": "6" * 64,
            "decoded_stl_solid_sha256": "7" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "stl_solids_use_current_normals_winding_watertight_edges_tolerance_volume_units_and_digests"
    ]


def test_v31_public_assembly_occurrence_location_density_unit_mass_inertia_parallel_axis_mismatch():
    reference, measured = _public_v31()
    measured["external_cad"][0][
        "assembly_occurrence_location_density_unit_suppression_mass_center_inertia_parallel_axis_shape_generation_identity"
    ].update(
        {
            "location_assembly_generation": "assembly-mass-180",
            "density_assembly_generation": "assembly-mass-179",
            "result_assembly_generation": "assembly-mass-178",
            "result_occurrence_ids": ["shaft-1", "housing-1"],
            "result_location_matrices": [[[1.0, 0.0, 0.0, 0.0]]],
            "result_densities_kg_m3": [2.7, 7.8],
            "result_density_unit": "g/cm^3",
            "result_suppressed_occurrence_ids": ["housing-1"],
            "result_part_volumes_m3": [2.0e-4],
            "result_assembly_mass_kg": 1.56,
            "result_center_of_mass_m": [0.1, 0.0, 0.0],
            "result_inertia_reference_frame": "shaft-local-origin",
            "result_parallel_axis_applied": False,
            "result_assembly_inertia_kg_m2": [[0.01]],
            "result_assembly_shape_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "assembly_occurrences_use_current_locations_densities_units_suppression_mass_center_parallel_axis_inertia_and_shape"
    ]


def test_v31_public_loft_sweep_topology_seam_face_lineage_orientation_volume_mismatch():
    reference, measured = _public_v31()
    measured["external_cad"][0][
        "loft_sweep_profile_order_seam_guide_orientation_face_lineage_shell_volume_shape_generation_identity"
    ].update(
        {
            "profile_loft_generation": "loft-lineage-180",
            "lineage_loft_generation": "loft-lineage-179",
            "result_loft_generation": "loft-lineage-178",
            "result_profile_ids": ["section-2", "section-1", "section-0"],
            "result_profile_parameters": [0.0, 0.8, 0.5],
            "result_seam_parameters": [0.5, 0.0, 0.25],
            "result_guide_orientation": "end_to_start_left_handed",
            "result_face_lineage_pairs": [["section-0:e2", "loft:f2"]],
            "result_face_lineage_sha256": "0" * 64,
            "result_shell_closed": False,
            "result_volume_m3": 0.0,
            "result_loft_shape_sha256": "1" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "lofts_use_current_profile_order_seams_guide_face_lineage_shell_volume_and_shape"
    ]


def test_v31_source_step_ap242_occurrence_transform_length_unit_product_structure_digest_mismatch():
    row = _source_v31()
    row["replay_identity"][
        "step_ap242_representation_context_external_owner_occurrence_transform_mass_product_structure_file_generation_identity"
    ].update(
        {
            "context_step_generation": "step-ap242-context-180",
            "owner_step_generation": "step-ap242-context-179",
            "result_step_generation": "step-ap242-context-178",
            "decoded_ap_schema": "AP203",
            "decoded_representation_contexts": [["ctx-inch", "inch", "degree"]],
            "decoded_part_owners": [["housing", "part-old"]],
            "decoded_occurrence_ids": ["shaft-2", "housing-1"],
            "decoded_occurrence_transform_sha256": [["housing-1", "2" * 64]],
            "decoded_occurrence_mass_properties": [["housing-1", 2700.0, 0.0, 0.0, 0.0]],
            "decoded_product_structure_sha256": "3" * 64,
            "decoded_step_file_sha256": "4" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_ap242_occurrences_use_current_schema_context_owners_transforms_mass_structure_and_file"
    ]


def test_v31_source_stl_tessellation_chord_angle_normal_watertight_volume_tolerance_mismatch():
    row = _source_v31()
    row["replay_identity"][
        "stl_tessellation_source_brep_chord_angle_facet_component_deviation_area_volume_digest_generation_identity"
    ].update(
        {
            "source_stl_generation": "stl-tessellation-error-180",
            "deviation_stl_generation": "stl-tessellation-error-179",
            "result_stl_generation": "stl-tessellation-error-178",
            "decoded_source_brep_sha256": "5" * 64,
            "decoded_chord_tolerance_m": 1.0e-2,
            "decoded_angular_tolerance_deg": 45.0,
            "decoded_facet_count": 120,
            "decoded_connected_component_count": 3,
            "decoded_maximum_surface_deviation_m": 2.0e-2,
            "decoded_surface_area_m2": 0.03,
            "decoded_volume_m3": 1.0e-3,
            "decoded_normal_table_sha256": "6" * 64,
            "decoded_stl_file_sha256": "7" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "stl_tessellations_use_current_brep_chord_angle_facets_components_deviation_area_volume_and_digests"
    ]


def test_v32_public_boolean_fuzzy_tolerance_topology_name_face_ancestry_volume_centroid_mismatch():
    reference, measured = _public_v32()
    measured["external_cad"][0][
        "boolean_fuzzy_tolerance_topology_name_face_ancestry_count_volume_centroid_shape_generation_identity"
    ].update(
        {
            "tolerance_generation": "boolean-topology-190",
            "ancestry_generation": "boolean-topology-189",
            "result_generation": "boolean-topology-188",
            "result_fuzzy_tolerance_m": 1.0e-3,
            "result_surviving_topology_names": ["body", "face_old"],
            "result_face_ancestry": [["box:f2", "result:f3"]],
            "result_solid_count": 2,
            "result_volume_m3": 8.0e-4,
            "result_centroid_m": [0.02, 0.01, 0.0],
            "result_boolean_shape_sha256": "d" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_results_use_current_fuzzy_tolerance_topology_names_face_ancestry_count_volume_centroid_and_shape"
    ]


def test_v32_public_sweep_frenet_frame_twist_transition_self_intersection_volume_mismatch():
    reference, measured = _public_v32()
    measured["external_cad"][0][
        "sweep_frenet_frame_twist_transition_profile_orientation_self_intersection_volume_owner_shape_generation_identity"
    ].update(
        {
            "frame_generation": "sweep-frame-190",
            "intersection_generation": "sweep-frame-189",
            "result_generation": "sweep-frame-188",
            "result_frame_convention": "fixed_global",
            "result_twist_parameters_rad": [0.4, 0.2, 0.0],
            "result_transition_mode": "transformed",
            "result_profile_orientation_signs": [1, -1, 1],
            "result_self_intersection": True,
            "result_sweep_volume_m3": 0.0,
            "result_shape_owner": "stale/body2",
            "result_sweep_shape_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "swept_solids_use_current_frenet_frame_twist_transition_orientation_intersection_volume_owner_and_shape"
    ]


def test_v32_source_step_assembly_instance_transform_unit_color_material_uuid_digest_mismatch():
    row = _source_v32()
    row["replay_identity"][
        "step_assembly_instance_transform_unit_color_material_uuid_component_volume_file_generation_identity"
    ].update(
        {
            "instance_generation": "step-assembly-190",
            "material_generation": "step-assembly-189",
            "result_generation": "step-assembly-188",
            "decoded_instance_ids": ["shaft-1", "housing-1"],
            "decoded_instance_transform_sha256": [["housing-1", "f" * 64]],
            "decoded_length_unit": "mm",
            "decoded_instance_colors_rgb": [["housing-1", 1.0, 0.0, 0.0]],
            "decoded_material_labels": [["housing-1", "steel"]],
            "decoded_product_uuids": [["housing-1", "uuid-old"]],
            "decoded_component_shape_sha256": [["housing-1", "0" * 64]],
            "decoded_total_volume_m3": 1200.0,
            "decoded_step_file_sha256": "1" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_assemblies_use_current_instances_transforms_units_colors_materials_uuids_components_volume_and_file"
    ]


def test_v32_source_stl_repair_tolerance_normal_duplicate_vertex_watertight_volume_digest_mismatch():
    row = _source_v32()
    row["replay_identity"][
        "stl_repair_merge_tolerance_normal_duplicate_boundary_watertight_volume_unit_file_generation_identity"
    ].update(
        {
            "tolerance_generation": "stl-repair-190",
            "watertight_generation": "stl-repair-189",
            "result_generation": "stl-repair-188",
            "decoded_merge_tolerance_m": 1.0e-3,
            "decoded_facet_normal_sha256": "2" * 64,
            "decoded_duplicate_vertex_count": 12,
            "decoded_boundary_edge_count": 8,
            "decoded_watertight_component_count": 3,
            "decoded_repaired_volume_m3": 1.0e-3,
            "decoded_length_unit": "mm",
            "decoded_source_stl_sha256": "3" * 64,
            "decoded_repaired_stl_sha256": "4" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "stl_repairs_use_current_tolerance_normals_duplicates_boundaries_watertight_volume_unit_and_files"
    ]


def test_v33_public_loft_section_guide_parameterization_seam_orientation_smooth_volume_mismatch():
    reference, measured = _public_v33()
    measured["external_cad"][0][
        "loft_section_guide_parameterization_seam_orientation_mode_intersection_volume_shape_generation_identity"
    ].update(
        {
            "section_generation": "guided-loft-200",
            "guide_generation": "guided-loft-199",
            "result_generation": "guided-loft-198",
            "result_section_order": ["section-0", "section-2", "section-1"],
            "result_wire_parameterization_sha256": "e" * 64,
            "result_guide_intersections": [["guide-0", "section-0"]],
            "result_seam_orientation_signs": [1, -1, 1],
            "result_loft_mode": "ruled",
            "result_self_intersection": True,
            "result_loft_volume_m3": 0.0,
            "result_loft_shape_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "guided_lofts_use_current_sections_parameterization_guides_seams_mode_intersection_volume_and_shape"
    ]


def test_v33_public_mass_property_inertia_origin_density_unit_principal_axis_parallel_axis_mismatch():
    reference, measured = _public_v33()
    measured["external_cad"][0][
        "mass_property_density_unit_origin_center_principal_axis_degeneracy_parallel_axis_owner_shape_generation_identity"
    ].update(
        {
            "density_generation": "mass-inertia-200",
            "parallel_axis_generation": "mass-inertia-199",
            "result_generation": "mass-inertia-198",
            "result_density_kg_m3": 7.8,
            "result_density_unit": "g/cm^3",
            "result_mass_kg": 0.0078,
            "result_inertia_origin_m": [0.1, 0.0, 0.0],
            "result_center_of_mass_m": [0.0, 0.0, 0.0],
            "result_principal_moments_kg_m2": [0.04, 0.03, 0.02],
            "result_principal_axes": [
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            "result_degeneracy_convention": "unsorted_left_handed",
            "result_shape_owner": "stale/body-2",
            "result_mass_shape_sha256": "0" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mass_properties_use_current_density_units_origin_center_principal_axes_parallel_axis_owner_and_shape"
    ]


def test_v33_source_brep_face_pcurve_orientation_location_tolerance_serialization_digest_mismatch():
    row = _source_v33()
    row["replay_identity"][
        "brep_face_pcurve_wire_orientation_location_tolerance_surface_serializer_shape_generation_identity"
    ].update(
        {
            "pcurve_generation": "brep-semantic-roundtrip-200",
            "location_generation": "brep-semantic-roundtrip-199",
            "result_generation": "brep-semantic-roundtrip-198",
            "decoded_face_pcurve_sha256": [[1, "1" * 64]],
            "decoded_wire_orientation_signs": [[1, -1], [2, 1]],
            "decoded_nested_location_sha256": [["root/part-1", "2" * 64]],
            "decoded_edge_tolerances_m": [[11, 1.0e-3]],
            "decoded_surface_types": [[1, "bspline"], [2, "plane"]],
            "decoded_serializer_version": "occt-brep-v2",
            "decoded_brep_shape_sha256": "3" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "brep_roundtrips_use_current_pcurves_wire_orientations_locations_tolerances_surfaces_serializer_and_shape"
    ]


def test_v33_source_gltf_node_transform_triangle_material_unit_tessellation_volume_digest_mismatch():
    row = _source_v33()
    row["replay_identity"][
        "gltf_node_hierarchy_transform_winding_material_unit_tessellation_volume_file_generation_identity"
    ].update(
        {
            "node_generation": "gltf-roundtrip-200",
            "tessellation_generation": "gltf-roundtrip-199",
            "result_generation": "gltf-roundtrip-198",
            "decoded_node_hierarchy": [["root", "mesh-1"]],
            "decoded_instance_transform_sha256": [["part-1", "4" * 64]],
            "decoded_triangle_winding_sha256": "5" * 64,
            "decoded_material_assignments": [["mesh-1", "default"]],
            "decoded_length_unit": "mm",
            "decoded_linear_deflection_m": 1.0e-2,
            "decoded_angular_deflection_rad": 0.5,
            "decoded_triangle_count": 128,
            "decoded_enclosed_volume_m3": 1.2e6,
            "decoded_gltf_file_sha256": "6" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "gltf_roundtrips_use_current_hierarchy_transforms_winding_materials_units_tessellation_volume_and_file"
    ]


def test_v34_public_assembly_mate_transform_cycle_frame_closure_mass_inertia_mismatch():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
    ]
    identity.update(
        {
            "mate_generation": "assembly-mate-210",
            "mass_generation": "assembly-mate-209",
            "result_generation": "assembly-mate-208",
            "result_part_ids": ["base", "tool"],
            "result_mate_edges": [["base", "arm"], ["tool", "base"]],
            "result_mate_cycle_transform": [[1.0, 0.0, 0.0, 0.1], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            "result_frame_determinants": [1.0, -1.0, 1.0],
            "result_part_masses_kg": [2.0, 3.0, 2.0],
            "result_part_centers_m": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
            "result_assembly_mass_kg": 5.0,
            "result_assembly_center_of_mass_m": [0.0, 0.0, 0.0],
            "result_part_inertia_global_kg_m2": identity["part_inertia_local_kg_m2"],
            "result_assembly_owner": "stale/root",
            "accepted_assembly_shape_sha256": "9" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "assembly_mates_close_current_cycles_frames_mass_centers_rotated_inertia_owner_and_shape"
    ]


def test_v34_public_shell_fillet_topology_euler_thickness_volume_area_inertia_convergence_mismatch():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
    ]
    identity.update(
        {
            "topology_generation": "shell-fillet-210",
            "volume_generation": "shell-fillet-209",
            "result_generation": "shell-fillet-208",
            "result_vertex_count": 15,
            "result_edge_count": 22,
            "result_face_count": 8,
            "result_euler_characteristic": 1,
            "result_edge_face_incidence_counts": [2] * 20 + [1, 3],
            "result_nominal_wall_thickness_m": 0.004,
            "result_wall_thickness_samples_m": [0.001, 0.004, 0.006],
            "result_removed_volume_m3": 0.7,
            "result_shell_volume_m3": 0.4,
            "result_surface_area_m2": 4.0,
            "result_inertia_tensor_kg_m2": [[0.2, 0.1, 0.0], [0.0, -0.3, 0.0], [0.0, 0.0, 0.4]],
            "result_convergence_tolerances_m": [1.0e-4, 1.0e-3, 1.0e-2],
            "result_convergence_volumes_m3": [0.2, 0.3, 0.5],
            "accepted_shell_brep_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "shell_fillets_use_current_euler_manifold_thickness_volume_area_inertia_convergence_and_brep"
    ]


def test_v34_source_step_unit_entity_color_assembly_transform_shape_owner_roundtrip_mismatch():
    row = _source_v34()
    row["replay_identity"][
        "step_unit_product_entity_color_assembly_transform_shape_validity_owner_file_result_generation_identity"
    ].update(
        {
            "unit_generation": "step-semantic-roundtrip-210",
            "entity_generation": "step-semantic-roundtrip-209",
            "result_generation": "step-semantic-roundtrip-208",
            "decoded_length_unit": "mm",
            "decoded_product_entities": [[1, "old-base"]],
            "decoded_entity_colors_rgb": [[1, [0.0, 0.0, 0.0]]],
            "decoded_assembly_transform_sha256": [[1, "b" * 64]],
            "decoded_shape_count": 1,
            "decoded_solid_validity": [[1, False]],
            "decoded_source_owner": "stale/root",
            "decoded_step_file_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_roundtrips_use_current_units_products_colors_transforms_shapes_validity_owner_and_file"
    ]


def test_v34_source_brep_periodic_seam_face_orientation_edge_pcurve_manifold_digest_mismatch():
    row = _source_v34()
    row["replay_identity"][
        "brep_periodic_face_seam_edge_orientation_pcurve_tolerance_manifold_serializer_shape_result_generation_identity"
    ].update(
        {
            "seam_generation": "brep-periodic-seam-210",
            "serializer_generation": "brep-periodic-seam-209",
            "result_generation": "brep-periodic-seam-208",
            "decoded_periodic_face_ids": [1],
            "decoded_seam_edge_multiplicity": [[1, 11, 1], [2, 12, 3]],
            "decoded_face_orientation_signs": [[1, -1], [2, 1]],
            "decoded_edge_pcurve_max_deviation_m": [[11, 1.0e-2]],
            "decoded_vertex_tolerances_m": [[101, 1.0e-2]],
            "decoded_edge_face_incidence_counts": [[11, 1], [12, 3]],
            "decoded_serializer_version": "occt-brep-v2",
            "decoded_periodic_brep_shape_sha256": "d" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "periodic_brep_roundtrips_use_current_seams_orientations_pcurves_tolerances_manifold_serializer_and_shape"
    ]


def test_v34_rejects_self_consistent_nonclosing_mate_cycle():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
    ]
    cycle = [[1.0, 0.0, 0.0, 0.1], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    identity["mate_cycle_transform"] = cycle
    identity["result_mate_cycle_transform"] = cycle
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_shell_volume_without_subtraction_closure():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
    ]
    identity["shell_volume_m3"] = 0.3
    identity["result_shell_volume_m3"] = 0.3
    identity["convergence_volumes_m3"][-1] = 0.3
    identity["result_convergence_volumes_m3"][-1] = 0.3
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_step_invalid_solid():
    row = _source_v34()
    identity = row["replay_identity"][
        "step_unit_product_entity_color_assembly_transform_shape_validity_owner_file_result_generation_identity"
    ]
    identity["solid_validity"] = [[1, True], [2, False]]
    identity["decoded_solid_validity"] = [[1, True], [2, False]]
    assert _source_result(row)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_nonmanifold_periodic_seam():
    row = _source_v34()
    identity = row["replay_identity"][
        "brep_periodic_face_seam_edge_orientation_pcurve_tolerance_manifold_serializer_shape_result_generation_identity"
    ]
    identity["edge_face_incidence_counts"] = [[11, 1], [12, 3]]
    identity["decoded_edge_face_incidence_counts"] = [[11, 1], [12, 3]]
    assert _source_result(row)["status"] == "needs_attention"


def test_v35_public_mirrored_pattern_instance_handedness_transform_mass_inertia_mismatch():
    reference, measured = _public_v35()
    identity = measured["external_cad"][0]["mirrored_pattern_occurrence_handedness_transform_suppression_volume_mass_center_inertia_owner_shape_result_generation_identity"]
    identity.update({"handedness_generation": "mirrored-pattern-220", "mass_generation": "mirrored-pattern-219", "result_generation": "mirrored-pattern-218", "result_occurrence_ids": ["part:0"], "result_transform_determinants": [1.0, 1.0], "result_handedness": ["right", "right"], "result_suppressed": [False, True], "result_occurrence_volumes_m3": [1.0, 0.5], "result_occurrence_masses_kg": [2.0, 1.0], "result_occurrence_centers_m": [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], "result_occurrence_inertia_principal_kg_m2": [[0.1, 0.2, -0.3]], "result_assembly_volume_m3": 1.5, "result_assembly_mass_kg": 3.0, "result_assembly_center_m": [1.5, 0.0, 0.0], "result_assembly_owner": "stale/root", "accepted_mirrored_shape_sha256": "9" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["mirrored_patterns_use_current_occurrences_handedness_transforms_suppression_mass_center_inertia_owner_and_shape"]


def test_v35_public_offset_thicken_curvature_selfintersection_wall_volume_topology_mismatch():
    reference, measured = _public_v35()
    identity = measured["external_cad"][0]["offset_thicken_curvature_sign_selfintersection_repair_thickness_volume_topology_convergence_brep_result_generation_identity"]
    identity.update({"curvature_generation": "offset-thicken-220", "topology_generation": "offset-thicken-219", "result_generation": "offset-thicken-218", "result_minimum_curvature_radius_m": 0.005, "result_offset_m": 0.02, "result_self_intersection_count": 3, "result_repair_mode": "discard_faces", "result_wall_thickness_samples_m": [0.001, 0.02, 0.04], "result_removed_volume_m3": 0.6, "result_thickened_volume_m3": 0.5, "result_solid_count": 2, "result_shell_count": 3, "result_convergence_tolerances_m": [1.0e-6, 1.0e-5, 1.0e-4], "result_convergence_volumes_m3": [0.3, 0.4, 0.6], "accepted_thicken_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["offset_thickens_use_current_curvature_sign_intersections_repair_thickness_volume_topology_convergence_and_brep"]


def test_v35_source_step_ap242_pmi_name_color_unit_occurrence_transform_roundtrip_mismatch():
    row = _source_v35(); identity = row["replay_identity"]["step_ap242_pmi_unit_name_color_occurrence_transform_validity_owner_file_result_generation_identity"]
    identity.update({"pmi_generation": "step-ap242-pmi-220", "occurrence_generation": "step-ap242-pmi-219", "result_generation": "step-ap242-pmi-218", "decoded_schema": "AP203", "decoded_length_unit": "mm", "decoded_pmi_annotations": [[1, "linear_dimension", "in", 25.0]], "decoded_product_names": [[1, "old"]], "decoded_colors_rgb": [[1, [0.0, 0.0, 0.0]]], "decoded_occurrence_paths": [[1, "old/root"]], "decoded_transform_sha256": [[1, "b" * 64]], "decoded_solid_validity": [[1, False]], "decoded_source_owner": "stale/root", "decoded_ap242_file_sha256": "c" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["step_ap242_roundtrips_use_current_pmi_units_names_colors_occurrences_transforms_validity_owner_and_file"]


def test_v35_source_dxf_profile_arc_bulge_winding_plane_unit_extrusion_roundtrip_mismatch():
    row = _source_v35(); identity = row["replay_identity"]["dxf_profile_unit_plane_layer_arc_bulge_loop_winding_extrusion_topology_owner_file_result_generation_identity"]
    identity.update({"arc_generation": "dxf-profile-220", "extrusion_generation": "dxf-profile-219", "result_generation": "dxf-profile-218", "decoded_length_unit": "in", "decoded_sketch_plane": "XZ", "decoded_layer_names": ["0"], "decoded_arc_bulges": [[11, -0.41421356237]], "decoded_loop_winding_signs": [["outer", -1], ["hole", 1]], "decoded_closed_loop_count": 1, "decoded_extrusion_height_m": 0.1, "decoded_profile_area_m2": 0.01, "decoded_solid_count": 2, "decoded_extruded_volume_m3": 0.1, "decoded_profile_owner": "stale/profile", "decoded_dxf_file_sha256": "d" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["dxf_profiles_use_current_units_plane_layers_arcs_bulges_winding_extrusion_topology_owner_and_file"]


def test_v35_rejects_self_consistent_mirror_handedness_error():
    reference, measured = _public_v35(); identity = measured["external_cad"][0]["mirrored_pattern_occurrence_handedness_transform_suppression_volume_mass_center_inertia_owner_shape_result_generation_identity"]
    identity["handedness"] = ["right", "right"]; identity["result_handedness"] = ["right", "right"]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_offset_volume_error():
    reference, measured = _public_v35(); identity = measured["external_cad"][0]["offset_thicken_curvature_sign_selfintersection_repair_thickness_volume_topology_convergence_brep_result_generation_identity"]
    identity["thickened_volume_m3"] = 0.4; identity["result_thickened_volume_m3"] = 0.4; identity["convergence_volumes_m3"][-1] = 0.4; identity["result_convergence_volumes_m3"][-1] = 0.4
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_ap242_invalid_solid():
    row = _source_v35(); identity = row["replay_identity"]["step_ap242_pmi_unit_name_color_occurrence_transform_validity_owner_file_result_generation_identity"]
    identity["solid_validity"] = [[1, True], [2, False]]; identity["decoded_solid_validity"] = [[1, True], [2, False]]
    assert _source_result(row)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_dxf_volume_error():
    row = _source_v35(); identity = row["replay_identity"]["dxf_profile_unit_plane_layer_arc_bulge_loop_winding_extrusion_topology_owner_file_result_generation_identity"]
    identity["extruded_volume_m3"] = 0.01; identity["decoded_extruded_volume_m3"] = 0.01
    assert _source_result(row)["status"] == "needs_attention"


def test_v36_public_fillet_chamfer_edge_selection_radius_tolerance_topology_volume_owner_mismatch():
    reference, measured = _public_v36()
    identity = measured["external_cad"][0]["fillet_chamfer_edge_selection_radius_distance_tolerance_euler_volume_area_owner_brep_result_generation_identity"]
    identity.update({"selection_generation": "fillet-chamfer-contract-230", "result_feature_kind": "chamfer", "result_selected_edge_ids": [99], "result_radius_m": -0.01, "result_distance_m": 0.02, "result_modeling_tolerance_m": 0.1, "result_topology_after_v_e_f": [8, 24, 10], "result_volume_after_m3": 1.2, "result_surface_area_after_m2": -1.0, "result_shape_owner": "stale:part", "accepted_feature_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fillet_chamfer_features_use_current_edges_size_tolerance_topology_volume_area_owner_and_brep"]


def test_v36_public_loft_section_order_orientation_continuity_volume_centroid_owner_mismatch():
    reference, measured = _public_v36()
    identity = measured["external_cad"][0]["loft_section_order_orientation_guide_continuity_volume_centroid_owner_brep_result_generation_identity"]
    identity.update({"section_generation": "loft-contract-230", "result_section_ids": ["wire:z2", "wire:z0"], "result_section_parameters": [1.0, 0.0], "result_section_orientation_signs": [1, -1], "result_guide_correspondence": [[0, 3], [1, 2]], "result_continuity": "C0", "result_loft_volume_m3": -0.75, "result_loft_centroid_m": [0.0, 0.0, 2.0], "result_loft_owner": "stale:loft", "accepted_loft_brep_sha256": "b" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["lofts_use_current_sections_order_orientation_guides_continuity_volume_centroid_owner_and_brep"]


def test_v36_source_step_assembly_hierarchy_transform_name_unit_occurrence_digest_mismatch():
    row = _source_v36(); identity = row["replay_identity"]["step_assembly_hierarchy_occurrence_transform_repeated_part_name_color_unit_owner_file_result_generation_identity"]
    identity.update({"hierarchy_generation": "step-assembly-contract-230", "decoded_occurrence_paths": ["root/old"], "decoded_part_ids": ["old"], "decoded_transform_sha256": ["c" * 64], "decoded_product_names": ["old"], "decoded_colors_rgb": [[0.0, 0.0, 0.0]], "decoded_length_unit": "mm", "decoded_assembly_owner": "stale:assembly", "decoded_step_file_sha256": "d" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["step_assemblies_use_current_hierarchy_occurrences_repeated_parts_transforms_names_colors_units_owner_and_file"]


def test_v36_source_boolean_tolerance_healing_sliver_nonmanifold_history_owner_mismatch():
    row = _source_v36(); identity = row["replay_identity"]["boolean_tolerance_healing_sliver_nonmanifold_operation_history_input_output_owner_brep_result_generation_identity"]
    identity.update({"healing_generation": "boolean-history-contract-230", "decoded_operation": "cut", "decoded_fuzzy_tolerance_m": -1.0, "decoded_healing_actions": ["discard"], "decoded_sliver_face_count": 3, "decoded_nonmanifold_edge_count": 2, "decoded_input_shape_ids": ["solid:a"], "decoded_output_shape_id": "solid:old", "decoded_operation_history": [["solid:x", "solid:old"]], "decoded_input_owner": "stale:inputs", "decoded_output_owner": "stale:result", "decoded_boolean_brep_sha256": "e" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["boolean_roundtrips_use_current_operation_tolerance_healing_slivers_manifold_history_owners_and_brep"]


def test_v36_rejects_self_consistent_fillet_euler_error():
    reference, measured = _public_v36(); identity = measured["external_cad"][0]["fillet_chamfer_edge_selection_radius_distance_tolerance_euler_volume_area_owner_brep_result_generation_identity"]
    identity["topology_after_v_e_f"] = [16, 24, 9]; identity["result_topology_after_v_e_f"] = [16, 24, 9]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v36_rejects_self_consistent_boolean_nonmanifold_output():
    row = _source_v36(); identity = row["replay_identity"]["boolean_tolerance_healing_sliver_nonmanifold_operation_history_input_output_owner_brep_result_generation_identity"]
    identity["nonmanifold_edge_count"] = 1; identity["decoded_nonmanifold_edge_count"] = 1
    assert _source_result(row)["status"] == "needs_attention"


def test_v37_public_mass_properties_centroid_inertia_principal_axes_placement_density_owner_mismatch():
    reference, measured = _public_v37(); identity = measured["external_cad"][0]["mass_properties_centroid_inertia_principal_axes_placement_density_shape_brep_generation_identity"]
    identity.update({"density_generation": "mass-placement-contract-240", "result_density_kg_m3": -7800.0, "result_mass_kg": -7.8, "result_centroid_world_m": [3.0, 2.0, 1.0], "result_inertia_world_kg_m2": [[-1.0, 2.0]], "result_principal_moments_kg_m2": [3.0, 2.0, 1.0], "result_principal_axes_world": [[-1.0, 0.0, 0.0]], "result_placement_transform": [[1.0]], "result_shape_owner": "stale:shape", "accepted_shape_brep_sha256": "9" * 64})
    result = _public_result(reference, measured); assert result["status"] == "needs_attention"; assert not result["checks"]["placed_mass_properties_use_current_density_mass_centroid_inertia_principal_axes_placement_owner_and_brep"]


def test_v37_public_shell_offset_thickness_normal_side_open_faces_topology_volume_owner_mismatch():
    reference, measured = _public_v37(); identity = measured["external_cad"][0]["shell_offset_thickness_normal_side_removed_face_topology_volume_input_owner_brep_generation_identity"]
    identity.update({"thickness_generation": "shell-offset-contract-240", "result_thickness_m": -0.002, "result_offset_side": "outside", "result_face_normals": [[2.0, 0.0, 0.0]], "result_removed_face_ids": [6, 6], "result_wall_topology_v_e_f": [0, 0, 0], "result_analytical_shell_volume_m3": -0.012, "result_input_owner": "stale:input", "accepted_shell_brep_sha256": "a" * 64})
    result = _public_result(reference, measured); assert result["status"] == "needs_attention"; assert not result["checks"]["shell_offsets_use_current_thickness_side_normals_removed_faces_topology_volume_owner_and_brep"]


def test_v37_source_sketch_constraint_dof_solver_status_reference_geometry_owner_digest_mismatch():
    row = _source_v37(); identity = row["replay_identity"]["sketch_constraint_dof_solver_reference_unit_owner_source_result_generation_identity"]
    identity.update({"constraint_generation": "sketch-solve-contract-240", "solved_constraint_ids": ["distance:old"], "solved_remaining_dof": 3, "solved_solver_status": "under_constrained", "solved_reference_geometry_ids": ["stale"], "solved_length_unit": "mm", "solved_sketch_owner": "stale:sketch", "solved_sketch_source_sha256": "b" * 64, "accepted_sketch_result_sha256": "c" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"; assert not result["checks"]["sketch_solves_use_current_constraints_dof_status_references_units_owner_source_and_result"]


def test_v37_source_topological_naming_edge_face_history_ocp_version_shape_owner_mismatch():
    row = _source_v37(); identity = row["replay_identity"]["topological_naming_edge_face_history_ocp_selector_shape_feature_source_brep_generation_identity"]
    identity.update({"edge_generation": "toponame-contract-240", "replayed_edge_names": ["edge:unknown"], "replayed_face_names": ["face:old"], "replayed_operation_history": ["box", "cut"], "replayed_ocp_version": "7.7.0", "replayed_selector_result": ["face:side"], "replayed_shape_generation_id": 40, "replayed_feature_owner": "stale:feature", "replayed_feature_source_sha256": "d" * 64, "accepted_feature_brep_sha256": "e" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"; assert not result["checks"]["topological_names_use_current_edges_faces_history_ocp_selector_shape_owner_source_and_brep"]


def test_v37_rejects_self_consistent_mass_density_closure_error():
    reference, measured = _public_v37(); identity = measured["external_cad"][0]["mass_properties_centroid_inertia_principal_axes_placement_density_shape_brep_generation_identity"]; identity["mass_kg"] = identity["result_mass_kg"] = 8.0; assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v37_rejects_self_consistent_duplicate_topological_name():
    row = _source_v37(); identity = row["replay_identity"]["topological_naming_edge_face_history_ocp_selector_shape_feature_source_brep_generation_identity"]; identity["edge_names"] = identity["replayed_edge_names"] = ["edge:fillet:0", "edge:fillet:0"]; assert _source_result(row)["status"] == "needs_attention"


def test_v38_public_revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_mismatch():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_shape_brep_generation_identity"
    ]
    identity.update(
        {
            "axis_generation": "revolve-contract-257",
            "topology_generation": "revolve-contract-256",
            "result_generation": "revolve-contract-255",
            "result_axis_direction": [1.0, 0.0, 0.0],
            "result_sweep_angle_deg": 180.0,
            "result_profile_axis_crossing": True,
            "result_profile_loop_orientation": "clockwise_inward",
            "result_analytic_volume_m3": -1.0,
            "result_centroid_world_m": [1.0, 0.0, 0.0],
            "result_solid_count": 2,
            "result_boundary_genus": 0,
            "result_boundary_euler_characteristic": 2,
            "result_shape_owner": "stale:revolve",
            "accepted_revolve_brep_sha256": "9" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "revolved_solids_use_current_axis_angle_profile_orientation_pappus_volume_centroid_topology_owner_and_brep"
    ]


def test_v38_public_involute_gear_module_teeth_pressure_angle_backlash_pitch_volume_owner_mismatch():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "involute_gear_module_teeth_pressure_backlash_pitch_base_addendum_periodicity_volume_shape_brep_generation_identity"
    ]
    identity.update(
        {
            "module_generation": "involute-gear-contract-257",
            "diameter_generation": "involute-gear-contract-256",
            "result_generation": "involute-gear-contract-255",
            "result_module_m": -2.0e-3,
            "result_tooth_count": 19,
            "result_pressure_angle_deg": 45.0,
            "result_backlash_m": -1.0e-4,
            "result_pitch_diameter_m": 0.1,
            "result_base_diameter_m": 0.2,
            "result_addendum_diameter_m": 0.03,
            "result_tooth_period_angle_deg": 20.0,
            "result_gear_volume_m3": -8.0e-6,
            "result_shape_owner": "stale:gear",
            "accepted_gear_brep_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "involute_gears_use_current_module_teeth_pressure_backlash_diameters_periodicity_volume_owner_and_brep"
    ]


def test_v38_source_brep_roundtrip_tolerance_ocp_version_subshape_counts_bounds_volume_digest_mismatch():
    row = _source_v38()
    identity = row["replay_identity"][
        "brep_roundtrip_tolerance_ocp_version_subshape_count_bounds_volume_shape_owner_source_restored_generation_identity"
    ]
    identity.update(
        {
            "tolerance_generation": "brep-roundtrip-contract-257",
            "topology_generation": "brep-roundtrip-contract-256",
            "result_generation": "brep-roundtrip-contract-255",
            "restored_serialization_tolerance_m": 1.0e-2,
            "restored_ocp_version": "7.7.0",
            "restored_subshape_counts": {
                "solid": 0,
                "shell": 1,
                "face": 5,
                "edge": 12,
                "vertex": 8,
            },
            "restored_bounding_box_min_m": [0.5, 0.25, 0.1],
            "restored_bounding_box_max_m": [-0.5, -0.25, 0.0],
            "restored_volume_m3": -1.0,
            "restored_shape_owner": "stale:brep",
            "restored_source_brep_sha256": "b" * 64,
            "accepted_restored_brep_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "brep_roundtrips_use_current_tolerance_ocp_topology_bounds_volume_owner_source_and_restored_shape"
    ]


def test_v38_source_svg_path_fillrule_curve_transform_unit_wire_face_extrusion_owner_mismatch():
    row = _source_v38()
    identity = row["replay_identity"][
        "svg_path_fillrule_curve_transform_unit_wire_face_extrusion_source_digest_generation_identity"
    ]
    identity.update(
        {
            "path_generation": "svg-extrusion-contract-257",
            "unit_generation": "svg-extrusion-contract-256",
            "result_generation": "svg-extrusion-contract-255",
            "replayed_path_commands": ["M", "L"],
            "replayed_fill_rule": "evenodd",
            "replayed_curve_transform": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            "replayed_document_unit": "px",
            "replayed_document_to_meter_scale": 1.0,
            "replayed_wire_closed": False,
            "replayed_face_orientation": "clockwise_negative",
            "replayed_extrusion_volume_m3": -2.0e-5,
            "replayed_source_owner": "stale:svg",
            "accepted_svg_result_sha256": "d" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "svg_extrusions_use_current_paths_fill_transform_units_wire_face_volume_owner_and_digests"
    ]


def test_v38_rejects_self_consistent_non_pappus_volume():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_shape_brep_generation_identity"
    ]
    identity["analytic_volume_m3"] = 1.0
    identity["result_analytic_volume_m3"] = 1.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_involute_pitch_error():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "involute_gear_module_teeth_pressure_backlash_pitch_base_addendum_periodicity_volume_shape_brep_generation_identity"
    ]
    identity["pitch_diameter_m"] = 0.1
    identity["result_pitch_diameter_m"] = 0.1
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_brep_euler_error():
    row = _source_v38()
    identity = row["replay_identity"][
        "brep_roundtrip_tolerance_ocp_version_subshape_count_bounds_volume_shape_owner_source_restored_generation_identity"
    ]
    counts = {"solid": 1, "shell": 1, "face": 5, "edge": 12, "vertex": 8}
    identity["subshape_counts"] = counts
    identity["restored_subshape_counts"] = counts
    assert _source_result(row)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_svg_reflection():
    row = _source_v38()
    identity = row["replay_identity"][
        "svg_path_fillrule_curve_transform_unit_wire_face_extrusion_source_digest_generation_identity"
    ]
    reflection = [[-1.0, 0.0, 10.0], [0.0, 1.0, 20.0], [0.0, 0.0, 1.0]]
    identity["curve_transform"] = reflection
    identity["replayed_curve_transform"] = reflection
    assert _source_result(row)["status"] == "needs_attention"


def test_v39_public_loft_multiwire_hole_correspondence_section_orientation_volume_topology_mismatch():
    reference, measured = _public_v39()
    value = measured["external_cad"][0][_LOFT_V39]
    value.update({"correspondence_generation": "multiwire-loft-270", "topology_generation": "multiwire-loft-269", "result_generation": "multiwire-loft-268", "result_outer_wire_ids": [31, 21, 11], "result_inner_wire_ids": [12, 32], "result_section_order": [0, 2, 1], "result_outer_orientation": ["cw", "ccw", "ccw"], "result_inner_orientation": ["ccw", "cw", "cw"], "result_hole_continuity": [[12, 32]], "result_self_intersection_free": False, "result_volume_m3": -2.0e-3, "result_centroid_m": [1.0, 0.0, 0.0], "result_solid_count": 2, "result_boundary_euler_characteristic": 2, "result_shape_owner": "stale:loft", "accepted_loft_brep_sha256": "9" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["multiwire_lofts_use_current_wire_correspondence_orientation_holes_intersection_mass_topology_owner_and_brep"]


def test_v39_public_offset_shell_thickness_join_selfintersection_area_mass_owner_mismatch():
    reference, measured = _public_v39()
    value = measured["external_cad"][0][_OFFSET]
    value.update({"thickness_generation": "offset-shell-270", "mass_generation": "offset-shell-269", "result_generation": "offset-shell-268", "result_signed_thickness_m": 2.0e-3, "result_offset_direction": "outward", "result_join_mode": "intersection", "result_self_intersection_repaired": False, "result_outer_area_m2": 0.8, "result_inner_area_m2": 1.0, "result_enclosed_volume_m3": -1.6e-3, "result_density_kg_per_m3": 1.0, "result_mass_kg": -1.6, "result_centroid_m": [1.0, 0.0, 0.0], "result_principal_inertia_kg_m2": [-0.01, 0.01, 0.02], "result_shape_owner": "stale:offset", "accepted_offset_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["offset_shells_use_current_signed_thickness_join_repair_areas_mass_inertia_owner_and_brep"]


def test_v39_source_occt_heal_tolerance_sew_orientation_stable_name_roundtrip_mismatch():
    row = _source_v39()
    value = row["replay_identity"][_HEAL]
    value.update({"tolerance_generation": "occt-heal-270", "name_generation": "occt-heal-269", "result_generation": "occt-heal-268", "replayed_healing_tolerance_m": 1.0e-2, "replayed_sewn_shell_count": 2, "replayed_solid_count": 0, "replayed_face_orientations": [1, -1], "replayed_removed_degenerate_edge_ids": [92, 99], "replayed_stable_subshape_names": {"face:top": 6}, "replayed_roundtrip_subshape_counts": {"solid": 0, "shell": 2, "face": 5, "edge": 10}, "replayed_shape_owner": "stale:heal", "replayed_healed_brep_sha256": "b" * 64, "accepted_heal_result_sha256": "c" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["occt_heals_use_current_tolerance_sewing_orientation_degenerate_edges_names_roundtrip_owner_and_digests"]


def test_v39_source_assembly_hierarchy_location_joint_axis_collision_bom_owner_mismatch():
    row = _source_v39()
    value = row["replay_identity"][_ASSEMBLY_V39]
    value.update({"hierarchy_generation": "assembly-hierarchy-270", "location_generation": "assembly-hierarchy-269", "result_generation": "assembly-hierarchy-268", "replayed_hierarchy": {"root": ["rotor"]}, "replayed_local_locations_m": {"rotor": [1.0, 0.0, 0.0]}, "replayed_global_locations_m": {"rotor": [0.0, 0.0, 0.2]}, "replayed_joint_axis": [1.0, 1.0, 0.0], "replayed_collision_pairs": [["shaft", "rotor"]], "replayed_component_quantities": {"rotor": 2}, "replayed_bom_identity": {"rotor": "OLD-001"}, "replayed_assembly_owner": "stale:assembly", "accepted_assembly_result_sha256": "d" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["assemblies_use_current_hierarchy_locations_joint_axis_collisions_quantities_bom_owner_and_result"]


def test_v39_rejects_self_consistent_loft_hole_discontinuity():
    reference, measured = _public_v39()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_LOFT_V39]["hole_continuity"] = [[12, 32]]
            row[_LOFT_V39]["result_hole_continuity"] = [[12, 32]]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_offset_mass_error():
    reference, measured = _public_v39()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_OFFSET]["mass_kg"] = 2.0
            row[_OFFSET]["result_mass_kg"] = 2.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_heal_orientation_flip():
    row = _source_v39()
    row["replay_identity"][_HEAL]["face_orientations"] = [1, 1, -1, 1]
    row["replay_identity"][_HEAL]["replayed_face_orientations"] = [1, 1, -1, 1]
    assert _source_result(row)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_assembly_location_error():
    row = _source_v39()
    row["replay_identity"][_ASSEMBLY_V39]["global_locations_m"]["rotor"] = [0.0, 0.0, 0.2]
    row["replay_identity"][_ASSEMBLY_V39]["replayed_global_locations_m"]["rotor"] = [0.0, 0.0, 0.2]
    assert _source_result(row)["status"] == "needs_attention"


def test_v40_public_draft_mismatch():
    reference, measured = _public_v40()
    value = measured["external_cad"][0][_DRAFT]
    value.update(
        {
            "plane_generation": "draft-solid-310",
            "topology_generation": "draft-solid-309",
            "result_generation": "draft-solid-308",
            "result_neutral_plane_origin_m": [0.0, 0.0, 1.0],
            "result_neutral_plane_normal": [0.0, 1.0, 0.0],
            "result_pull_direction": [0.0, 0.0, -1.0],
            "result_signed_draft_angle_rad": -0.05235987755982988,
            "result_selected_face_ids": [1, 3],
            "result_tangent_continuity": False,
            "result_topology_signature": {"solid": 2, "face": 8},
            "result_volume_m3": -1.9e-3,
            "result_centroid_m": [1.0, 0.0, 0.0],
            "result_shape_owner": "stale:draft",
            "accepted_draft_brep_sha256": "9" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "drafts_use_current_neutral_plane_pull_angle_faces_tangency_topology_mass_owner_and_brep"
    ]


def test_v40_public_thread_mismatch():
    reference, measured = _public_v40()
    value = measured["external_cad"][0][_THREAD]
    value.update(
        {
            "pitch_generation": "modeled-thread-310",
            "profile_generation": "modeled-thread-309",
            "result_generation": "modeled-thread-308",
            "result_pitch_m": 1.0e-3,
            "result_handedness": "left",
            "result_flank_angle_rad": 0.5,
            "result_profile_type": "square",
            "result_major_diameter_m": 1.6e-2,
            "result_minor_diameter_m": 2.1e-2,
            "result_thread_length_m": 1.5e-2,
            "result_turn_count": 4.0,
            "result_runout_length_m": 3.0e-2,
            "result_self_intersection_free": False,
            "result_volume_m3": -5.0e-6,
            "result_shape_owner": "stale:thread",
            "accepted_thread_brep_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "threads_use_current_pitch_handedness_profile_diameters_runout_intersection_mass_owner_and_brep"
    ]


def test_v40_source_obb_inertia_mismatch():
    row = _source_v40()
    value = row["replay_identity"][_OBB]
    value.update(
        {
            "box_generation": "obb-inertia-310",
            "transform_generation": "obb-inertia-309",
            "result_generation": "obb-inertia-308",
            "replayed_obb_center_m": [3.0, 2.0, 1.0],
            "replayed_obb_half_extents_m": [0.5, -0.25, 0.125],
            "replayed_obb_axes": [[1.0, 0.0, 0.0]] * 3,
            "replayed_principal_moments_kg_m2": [0.025, -0.02, 0.01],
            "replayed_principal_axes": [[1.0, 0.0, 0.0]] * 3,
            "replayed_center_of_mass_world_m": [0.0, 0.0, 0.0],
            "replayed_local_to_world_transform": [[1.0, 0.0], [0.0, 1.0]],
            "replayed_length_unit": "mm",
            "replayed_inertia_unit": "g*mm^2",
            "replayed_shape_owner": "stale:obb",
            "accepted_obb_result_sha256": "b" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mass_property_replays_use_current_obb_principal_inertia_axes_com_transform_units_owner_and_result"
    ]


def test_v40_source_tessellation_mismatch():
    row = _source_v40()
    value = row["replay_identity"][_TESSELLATION]
    value.update(
        {
            "linear_generation": "tessellation-310",
            "orientation_generation": "tessellation-309",
            "result_generation": "tessellation-308",
            "replayed_linear_deflection_m": 1.0e-2,
            "replayed_angular_deflection_rad": 1.2,
            "replayed_triangle_count": 10,
            "replayed_vertex_coordinates_m": [[0.0, 0.0, 0.0]],
            "replayed_outward_orientation": False,
            "replayed_watertight": False,
            "replayed_nonmanifold_edge_count": 4,
            "replayed_stl_owner": "stale:stl",
            "replayed_source_brep_sha256": "c" * 64,
            "replayed_stl_sha256": "d" * 64,
            "accepted_tessellation_result_sha256": "e" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "tessellation_replays_use_current_deflections_triangles_vertices_orientation_watertight_stl_brep_owner_and_result"
    ]


def test_v40_rejects_self_consistent_draft_frame_misalignment():
    reference, measured = _public_v40()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_DRAFT]["pull_direction"] = [1.0, 0.0, 0.0]
            row[_DRAFT]["result_pull_direction"] = [1.0, 0.0, 0.0]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_thread_length_pitch_error():
    reference, measured = _public_v40()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_THREAD]["turn_count"] = 9.0
            row[_THREAD]["result_turn_count"] = 9.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_left_handed_obb_axes():
    row = _source_v40()
    axes = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    row["replay_identity"][_OBB]["obb_axes"] = axes
    row["replay_identity"][_OBB]["replayed_obb_axes"] = axes
    assert _source_result(row)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_duplicate_tessellation_vertex():
    row = _source_v40()
    vertices = row["replay_identity"][_TESSELLATION]["vertex_coordinates_m"]
    vertices[-1] = vertices[0]
    row["replay_identity"][_TESSELLATION]["replayed_vertex_coordinates_m"] = vertices
    assert _source_result(row)["status"] == "needs_attention"


def test_v41_public_loft_mismatch():
    reference, measured = _public_v41()
    measured["external_cad"][0][_LOFT_V41].update(
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
    row["replay_identity"][_ASSEMBLY_V41].update(
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
    row["replay_identity"][_STEP_V41].update(
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
            row[_LOFT_V41]["section_parameters"] = [0.0, 0.8, 0.7]
            row[_LOFT_V41]["result_section_parameters"] = [0.0, 0.8, 0.7]
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
    value = row["replay_identity"][_ASSEMBLY_V41]
    value["part_transforms"]["slider"] = transform
    value["replayed_part_transforms"]["slider"] = transform
    value["overlap_volumes_m3"] = [1.0e-6]
    value["replayed_overlap_volumes_m3"] = [1.0e-6]
    assert _source_result(row)["status"] == "needs_attention"


def test_v41_rejects_self_consistent_step_color_and_label_reference_error():
    row = _source_v41()
    value = row["replay_identity"][_STEP_V41]
    value["colors_rgb"]["part:1"] = [2.0, 0.0, 0.0]
    value["replayed_colors_rgb"]["part:1"] = [2.0, 0.0, 0.0]
    value["subshape_labels"]["part:9/face:1"] = "ghost"
    value["replayed_subshape_labels"]["part:9/face:1"] = "ghost"
    assert _source_result(row)["status"] == "needs_attention"


def test_v42_public_helical_sweep_mismatch():
    reference, measured = _public_v42()
    measured["external_cad"][0][_HELIX].update({"frame_generation": "helical-sweep-724", "result_pitch_m": -1.0, "result_turns": 0.0, "result_axial_rise_m": -1.0, "result_profile_frame": "fixed", "result_frenet_transport": False, "result_torsion_per_m": -1.0, "result_solid_valid": False, "result_volume_m3": -1.0, "result_centroid_m": [9.0], "result_shape_owner": "stale:helix", "accepted_helical_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["helical_sweeps_use_current_pitch_turns_frame_transport_torsion_validity_mass_owner_and_brep"]


def test_v42_public_fuzzy_boolean_mismatch():
    reference, measured = _public_v42()
    measured["external_cad"][0][_BOOLEAN_V42].update({"sliver_generation": "fuzzy-boolean-724", "result_fuzzy_tolerance_m": 1.0e-2, "result_sliver_face_count_after": 9, "result_topology_signature": {"solid": 0}, "result_solid_valid": False, "result_volume_m3": -1.0, "result_surface_area_m2": -1.0, "result_shape_owner": "stale:boolean", "accepted_boolean_brep_sha256": "b" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fuzzy_booleans_use_current_tolerance_slivers_topology_validity_mass_owner_and_brep"]


def test_v42_source_selector_cache_mismatch():
    row = _source_v42()
    row["replay_identity"][_SELECTOR].update({"renumber_generation": "selector-cache-724", "replayed_topology_renumber_map": {}, "replayed_geometry_generation_id": 724, "replayed_selector_predicate": "Area<0", "replayed_selected_feature_ids": ["face:11"], "replayed_parent_shape_owner": "stale:parent", "accepted_selector_result_sha256": "c" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["selector_caches_use_current_renumbering_geometry_predicate_features_parent_owner_and_result"]


def test_v42_source_step_assembly_mismatch():
    row = _source_v42()
    row["replay_identity"][_STEP_V42].update({"unit_generation": "step-assembly-724", "replayed_length_unit": "m", "replayed_unit_scale_to_m": 1.0, "replayed_part_names": {"part:1": "old"}, "replayed_part_colors_rgb": {"part:1": [2.0, -1.0, 0.0]}, "replayed_part_transforms_in_source_units": {}, "replayed_part_shape_ids": {}, "replayed_assembly_hierarchy": {}, "replayed_export_owner": "stale:export", "replayed_step_file_sha256": "d" * 64, "accepted_assembly_result_sha256": "e" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["step_assembly_exports_use_current_units_colors_names_transforms_shapes_hierarchy_owner_and_file"]


def test_v42_rejects_self_consistent_wrong_helical_path_length():
    reference, measured = _public_v42()
    for rows in [reference, *measured.values()]:
        for row in rows:
            value = row[_HELIX]
            wrong = 2.0 * value["path_length_m"]
            value["path_length_m"] = value["result_path_length_m"] = wrong
            value["volume_m3"] = value["result_volume_m3"] = value["profile_area_m2"] * wrong
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_oversized_fuzzy_tolerance():
    reference, measured = _public_v42()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_BOOLEAN_V42]["fuzzy_tolerance_m"] = row[_BOOLEAN_V42]["result_fuzzy_tolerance_m"] = 5.0e-5
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_selector_using_pre_renumber_features():
    row = _source_v42()
    value = row["replay_identity"][_SELECTOR]
    value["selected_feature_ids"] = value["replayed_selected_feature_ids"] = ["face:11", "face:12"]
    assert _source_result(row)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_nonrigid_step_transform():
    row = _source_v42()
    value = row["replay_identity"][_STEP_V42]
    transform = [[2, 0, 0, 100], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    value["part_transforms_in_source_units"]["part:2"] = transform
    value["replayed_part_transforms_in_source_units"]["part:2"] = transform
    assert _source_result(row)["status"] == "needs_attention"


def test_v43_public_gear_mismatch():
    reference, measured = _public_v43()
    measured["external_cad"][0][_GEAR]["result_pressure_angle_deg"] = 25.0
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["involute_gears_use_current_module_teeth_pressure_angle_backlash_volume_inertia_owner_and_brep"]


def test_v43_public_sheet_mismatch():
    reference, measured = _public_v43()
    measured["external_cad"][0][_SHEET]["result_k_factor"] = 1.2
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["sheet_metal_bends_use_current_radius_kfactor_thickness_neutral_axis_volume_area_owner_and_brep"]


def test_v43_source_location_mismatch():
    row = _source_v43()
    row["replay_identity"][_LOCATION]["replayed_translation_frame"] = "world"
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["locations_use_current_local_global_composition_rotation_order_subshape_owner_and_result"]


def test_v43_source_rebuild_mismatch():
    row = _source_v43()
    row["replay_identity"][_REBUILD]["replayed_invalidated_properties"] = []
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["parametric_rebuilds_use_current_dependencies_cache_invalidation_topology_owner_and_result"]


def test_v43_rejects_self_consistent_wrong_gear_inertia():
    reference, measured = _public_v43()
    for rows in [reference, *measured.values()]:
        for row in rows:
            wrong = [[1.0e-8, 0.0, 0.0], [0.0, 2.0e-8, 0.0], [0.0, 0.0, -3.0e-8]]
            row[_GEAR]["inertia_tensor_kg_m2"] = wrong
            row[_GEAR]["result_inertia_tensor_kg_m2"] = wrong
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v43_rejects_self_consistent_stale_parametric_cache():
    row = _source_v43()
    value = row["replay_identity"][_REBUILD]
    value["cache_key"] = value["replayed_cache_key"] = "length=0.2;radius=0.02"
    value["dependency_values"] = value["replayed_dependency_values"] = {"length_m": 0.2, "radius_m": 0.02}
    value["dependency_identity_sha256"] = value["replayed_dependency_identity_sha256"] = "a" * 64
    assert _source_result(row)["status"] == "needs_attention"


def test_v44_rejects_boolean_shell_fillet_mismatch():
    reference, measured = _public_v44()
    measured["external_cad"][0][_BOOLEAN_V44]["result_volume_m3"] = -1.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v44_rejects_loft_orientation_mismatch():
    reference, measured = _public_v44()
    measured["external_cad"][0][_LOFT_V44]["result_section_orientation"] = "inconsistent"
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v44_rejects_sketch_replay_mismatch():
    row = _source_v44()
    row["replay_identity"][_SKETCH]["replayed_constraint_order"] = ["distance:2", "coincident:1", "horizontal:3"]
    assert _source_result(row)["status"] == "needs_attention"


def test_v44_rejects_step_export_mismatch():
    row = _source_v44()
    row["replay_identity"][_STEP_V44]["replayed_length_unit"] = "m"
    assert _source_result(row)["status"] == "needs_attention"


def test_v44_rejects_malformed_numeric_identity_without_raising():
    reference, measured = _public_v44()
    measured["external_cad"][0][_BOOLEAN_V44]["volume_m3"] = {"bad": "value"}
    assert _public_result(reference, measured)["status"] == "needs_attention"

    row = _source_v44()
    row["replay_identity"][_STEP_V44]["face_count"] = {"bad": "value"}
    assert _source_result(row)["status"] == "needs_attention"


def test_v44_rejects_self_consistent_invalid_mass_properties_and_owner():
    reference, measured = _public_v44()
    identity = measured["external_cad"][0][_BOOLEAN_V44]
    identity["center_of_mass_m"] = [float("nan"), 0.0, 0.0]
    identity["result_center_of_mass_m"] = identity["center_of_mass_m"]
    identity["shape_owner"] = identity["result_shape_owner"] = ""
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v44_rejects_self_consistent_invalid_loft_geometry():
    reference, measured = _public_v44()
    identity = measured["external_cad"][0][_LOFT_V44]
    identity["section_area_m2"] = [1.0e-4, -1.0, 1.0e-4]
    identity["result_section_area_m2"] = identity["section_area_m2"]
    identity["inertia_tensor_kg_m2"] = [[-1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    identity["result_inertia_tensor_kg_m2"] = identity["inertia_tensor_kg_m2"]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v45_positive_identity_contracts():
    reference, measured = _public_v45()
    assert _public_result(reference, measured)["status"] == "ok"
    source = _source_result(_source_v45())
    assert source["status"] == "ok"
    assert source["warnings"] == []


def test_v45_rejects_boolean_brep_mismatch():
    reference, measured = _public_v45()
    measured["external_cad"][0]["boolean_fillet_shell_massproperties_volume_area_brep_owner_identity"]["result_volume_m3"] = -1.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v45_rejects_loft_frame_mismatch():
    reference, measured = _public_v45()
    measured["external_cad"][0]["loft_sweep_section_frame_tangent_continuity_inertia_export_digest_identity"]["result_section_orientation"] = "inconsistent"
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v45_rejects_sketch_cache_mismatch():
    row = _source_v45()
    row["replay_identity"]["sketch_constraint_order_plane_frame_solver_cache_shape_generation_owner_identity"]["replayed_parameter_cache_key"] = "old"
    assert _source_result(row)["status"] == "needs_attention"


def test_v45_rejects_step_topology_mismatch():
    row = _source_v45()
    row["replay_identity"]["step_units_tessellation_tolerance_face_topology_brep_export_owner_identity"]["replayed_length_unit"] = "m"
    assert _source_result(row)["status"] == "needs_attention"


def test_v46_public_and_source_positive_identity_controls():
    assert validate_public_identity_v46(_public_payload())["status"] == "ok"
    assert validate_source_identity_v46(_source_payload())["status"] == "ok"


def test_v46_public_and_source_mutations_are_rejected():
    public = _public_payload()
    public["reference"][0][PLACEMENT]["result_unit_scale_to_si"] = float("nan")
    assert validate_public_identity_v46(public)["status"] == "needs_attention"

    public = _public_payload()
    public["reference"][0][BOOLEAN_V46]["result_partial_shape_status"] = "partial"
    assert validate_public_identity_v46(public)["status"] == "needs_attention"

    source = _source_payload()
    source["replay_identity"][SKETCH_V46]["result_plane"] = "YZ"
    assert validate_source_identity_v46(source)["status"] == "needs_attention"

    source = _source_payload()
    source["replay_identity"][STEP_V46]["result_checksum_sha256"] = "6" * 64
    assert validate_source_identity_v46(source)["status"] == "needs_attention"


def test_v47_positive_replays_are_accepted() -> None:
    public, source = _payloads_v47()
    assert validate_public_identity_v47(public)["status"] == "ok"
    assert validate_source_identity_v47(source)["status"] == "ok"


def test_v47_compound_and_assembly_mutations_are_rejected() -> None:
    public, _ = _payloads_v47()
    public["reference"][0][COMPOUND]["result_child_ids"] = ["body:c", "body:a", "body:b"]
    public["reference"][0][ASSEMBLY_V47]["result_local_to_world_transform_sha256"] = "a" * 64
    assert validate_public_identity_v47(public)["status"] == "needs_attention"


def test_v47_roundtrip_and_external_reference_mutations_are_rejected() -> None:
    _, source = _payloads_v47()
    source["replay_identity"][ROUNDTRIP]["result_duplicate_uuid_count"] = 1
    source["replay_identity"][EXTERNAL]["result_dependency_cycle_count"] = 1
    assert validate_source_identity_v47(source)["status"] == "needs_attention"


def test_v48_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _payloads_v48()
    assert validate_public_identity_v48(public)["status"] == "ok"
    assert validate_source_identity_v48(source)["status"] == "ok"


def test_v48_configuration_cache_and_loft_mutations_are_rejected() -> None:
    public, _ = _payloads_v48()
    public["reference"][0][SUPPRESSION]["cache_owner"] = "build:feature-config-v48-old"
    public["reference"][0][LOFT_V48]["result_self_intersection_count"] = 1
    result = validate_public_identity_v48(public)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) == {
        "v48_feature_configuration_mass_cache_owner",
        "v48_loft_profile_seam_correspondence_owner",
    }


def test_v48_import_metadata_and_selector_mutations_are_rejected() -> None:
    _, source = _payloads_v48()
    source["replay_identity"][IMPORT]["result_source_revision"] = "source-rev-v48-old"
    source["replay_identity"][SELECTOR_V48]["result_cardinality"] = 1
    result = validate_source_identity_v48(source)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) == {
        "v48_import_unit_metadata_subshape_owner",
        "v48_selector_cardinality_witness_history",
    }


def test_v49_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _payloads_v49()
    assert validate_public_identity_v49(public)["status"] == "ok"
    assert validate_source_identity_v49(source)["status"] == "ok"


def test_v49_public_assembly_and_sketch_mutations_are_rejected() -> None:
    public, _ = _payloads_v49()
    public["reference"][0][ASSEMBLY_V49]["result_assembly_mass_kg"] = 13.0
    public["reference"][0][SKETCH_V49]["result_remaining_dof"] = 2
    result = validate_public_identity_v49(public)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v49_source_step_and_boolean_mutations_are_rejected() -> None:
    _, source = _payloads_v49()
    source["replay_identity"][STEP_V49]["result_ap_schema"] = "AP203"
    source["replay_identity"][BOOLEAN_V49]["result_cache_shape_sha256"] = "a" * 64
    result = validate_source_identity_v49(source)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v50_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _payloads_v50()
    assert validate_public_identity_v50(public)["status"] == "ok"
    assert validate_source_identity_v50(source)["status"] == "ok"


def test_v50_public_mutations_are_rejected() -> None:
    public, _ = _payloads_v50()
    public["reference"][0][MATE_V50]["result_remaining_dof"] = 5
    public["reference"][0][FEATURE]["result_fillet_radius_m"] = 0.02
    result = validate_public_identity_v50(public)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v50_source_mutations_are_rejected() -> None:
    _, source = _payloads_v50()
    source["replay_identity"][HEALING]["result_solid_count"] = 0
    source["replay_identity"][STL]["result_length_unit"] = "mm"
    result = validate_source_identity_v50(source)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v50_self_consistent_invalid_frames_and_feature_history_are_rejected() -> None:
    public, _ = _payloads_v50()
    for row in [public["reference"][0], public["measured"]["cad"][0]]:
        row[MATE_V50]["mate_frames"]["shaft"] = row[MATE_V50]["result_mate_frames"]["shaft"] = [0, 0, 0.02, 2, 0, 0, 0]
        row[FEATURE]["topology_history"] = row[FEATURE]["result_topology_history"] = {"edge:11": ["edge:31"]}
    assert validate_public_identity_v50(public)["status"] == "needs_attention"


def test_v50_self_consistent_open_shell_and_bad_stl_normals_are_rejected() -> None:
    _, source = _payloads_v50()
    healing = source["replay_identity"][HEALING]
    healing["solid_count"] = healing["result_solid_count"] = 0
    stl = source["replay_identity"][STL]
    stl["normal_counts"] = stl["result_normal_counts"] = {"outward": 12000, "inward": 480, "degenerate": 0}
    assert validate_source_identity_v50(source)["status"] == "needs_attention"


def test_v51_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _records()
    assert validate_public_identity_v51(public)["status"] == "ok"
    assert validate_source_identity_v51(source)["status"] == "ok"


def test_v51_public_mutations_are_rejected() -> None:
    public, _ = _records()
    public["reference"][0]["mass_properties_frame_inertia_parallel_axis_density_shape_owner_identity"]["result_density_kg_m3"] = 2700.0
    public["reference"][0]["sweep_profile_path_trihedron_transition_selfintersection_history_owner_identity"]["result_self_intersections"] = ["edge:99"]
    assert validate_public_identity_v51(public)["status"] == "needs_attention"


def test_v51_source_mutations_are_rejected() -> None:
    _, source = _records()
    source["replay_identity"]["step_external_reference_occurrence_name_schema_unit_color_owner_identity"]["result_step_schema"] = "AP203"
    source["replay_identity"]["brep_occt_version_location_precision_triangulation_cache_owner_identity"]["result_occt_version"] = "7.8.0"
    assert validate_source_identity_v51(source)["status"] == "needs_attention"


def test_v51_invalid_canonical_records_are_rejected() -> None:
    public, source = _records()
    public["reference"][0]["mass_properties_frame_inertia_parallel_axis_density_shape_owner_identity"]["shifted_inertia_kg_m2"] = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.25]]
    source["replay_identity"]["step_external_reference_occurrence_name_schema_unit_color_owner_identity"]["occurrence_colors"] = {"occurrence:1": [2.0, 0.0, 0.0]}
    assert validate_public_identity_v51(public)["status"] == "needs_attention"
    assert validate_source_identity_v51(source)["status"] == "needs_attention"


def test_v52_positive_public_and_source_replays_are_accepted():
    public, source = _payloads_v52()
    assert validate_public_identity_v52(public)["status"] == "ok"
    assert validate_source_identity_v52(source)["status"] == "ok"


def test_v52_public_mutations_are_rejected():
    public, _ = _payloads_v52(); value = deepcopy(public); value["reference"][0][SELECTOR_V52]["result_topology_revision"] = "topology:stale"
    assert validate_public_identity_v52(value)["status"] == "needs_attention"


def test_v52_source_mutations_are_rejected():
    _, source = _payloads_v52(); value = deepcopy(source); value["replay_identity"][GLTF]["replayed_axis_convention"] = "Z_up_left_handed"
    assert validate_source_identity_v52(value)["status"] == "needs_attention"


def test_v52_invalid_canonical_records_are_rejected():
    public, _ = _payloads_v52(); value = deepcopy(public); value["reference"][0][WORKPLANE]["wire_closed"] = False; value["reference"][0][WORKPLANE]["result_wire_closed"] = False
    assert validate_public_identity_v52(value)["status"] == "needs_attention"


def test_v53_positive_public_and_source_replays_are_accepted():
    public, source = _payloads_v53()
    assert validate_public_identity_v53(public)["status"] == "ok"
    assert validate_source_identity_v53(source)["status"] == "ok"


def test_v53_frozen_mutations_are_rejected():
    public, source = _payloads_v53(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][MATE_V53]["result_handedness"] = "left"
    public["reference"][0][STEP_V53]["result_length_unit"] = "inch"
    source["replay_identity"][HISTORY]["replayed_face_ancestry"] = {"face:r1": ["face:stale"]}
    source["replay_identity"][MASS]["replayed_density_kg_m3"] = 2700.0
    assert validate_public_identity_v53(public)["status"] == "needs_attention"
    assert validate_source_identity_v53(source)["status"] == "needs_attention"


def test_v53_self_consistent_nonphysical_records_are_rejected():
    public, source = _payloads_v53(); public = deepcopy(public); source = deepcopy(source)
    bad_colors = {"component:a": [1.2, 0.0, 0.0], "component:b": [0.0, 0.0, 0.0]}
    public["reference"][0][STEP_V53]["component_colors"] = public["reference"][0][STEP_V53]["result_component_colors"] = bad_colors
    bad_tensor = [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.05]]
    source["replay_identity"][MASS]["inertia_tensor_kg_m2"] = source["replay_identity"][MASS]["replayed_inertia_tensor_kg_m2"] = bad_tensor
    assert validate_public_identity_v53(public)["status"] == "needs_attention"
    assert validate_source_identity_v53(source)["status"] == "needs_attention"


def test_v54_positive_public_and_source_identities_are_accepted():
    public, source = _payloads_v54()
    assert validate_public_identity_v54(public)["status"] == "ok"
    assert validate_source_identity_v54(source)["status"] == "ok"


def test_v54_frozen_mutations_are_rejected():
    public, source = _payloads_v54(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][ASSEMBLY_V54]["result_center_of_mass_m"] = [1.0, 0.0, 0.0]
    public["reference"][0][LOFT_V54]["result_section_correspondence"] = []
    source["replay_identity"][STEP_V54]["replayed_schema"] = "AP203"
    source["replay_identity"][TESSELLATION]["replayed_triangle_indices"] = [[0, 1, 9]]
    assert validate_public_identity_v54(public)["status"] == "needs_attention"
    assert validate_source_identity_v54(source)["status"] == "needs_attention"


def test_v54_self_consistent_nonphysical_records_are_rejected():
    public, source = _payloads_v54(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][ASSEMBLY_V54]["located_solids"]["solid:a"]["quaternion_wxyz"] = [2.0, 0.0, 0.0, 0.0]
    public["reference"][0][ASSEMBLY_V54]["result_located_solids"] = public["reference"][0][ASSEMBLY_V54]["located_solids"]
    public["reference"][0][LOFT_V54]["resulting_topology"]["edges"] = 11
    public["reference"][0][LOFT_V54]["result_resulting_topology"] = public["reference"][0][LOFT_V54]["resulting_topology"]
    source["replay_identity"][STEP_V54]["component_colors"]["part:a"] = [1.2, 0.0, 0.0]
    source["replay_identity"][STEP_V54]["replayed_component_colors"] = source["replay_identity"][STEP_V54]["component_colors"]
    source["replay_identity"][TESSELLATION]["triangle_indices"] = source["replay_identity"][TESSELLATION]["replayed_triangle_indices"] = [[0, 1, 9]]
    assert validate_public_identity_v54(public)["status"] == "needs_attention"
    assert validate_source_identity_v54(source)["status"] == "needs_attention"


def test_v54_malformed_nested_values_reject_without_raising():
    public, source = _payloads_v54(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][ASSEMBLY_V54]["solid_densities_kg_m3"] = {"solid:a": [7850.0]}
    source["replay_identity"][TESSELLATION]["triangle_indices"] = [[[0], 1, 2]]
    source["replay_identity"][TESSELLATION]["triangle_orientations"] = [[1]]
    assert validate_public_identity_v54(public)["status"] == "needs_attention"
    assert validate_source_identity_v54(source)["status"] == "needs_attention"


def test_v55_positive_public_and_source_identities_are_accepted():
    public, source = _payloads_v55()
    assert validate_public_identity_v55(public)["status"] == "ok"
    assert validate_source_identity_v55(source)["status"] == "ok"


def test_v55_frozen_mutations_are_rejected():
    public, source = _payloads_v55(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][THREAD]["result_pitch_m"] = 3.0e-3
    public["reference"][0][SHEET]["result_neutral_axis_factor"] = 1.5
    source["replay_identity"][PMI]["replayed_length_unit"] = "inch"
    source["replay_identity"][BREP_V55]["replayed_closed_volume_count"] = 0
    assert validate_public_identity_v55(public)["status"] == "needs_attention"
    assert validate_source_identity_v55(source)["status"] == "needs_attention"


def test_v55_self_consistent_bad_topology_or_bend_allowance_is_rejected():
    public, _ = _payloads_v55(); public = deepcopy(public)
    public["reference"][0][THREAD]["thread_topology"]["edges"] = 23
    public["reference"][0][THREAD]["result_thread_topology"] = public["reference"][0][THREAD]["thread_topology"]
    public["reference"][0][SHEET]["bend_allowance_m"] = public["reference"][0][SHEET]["result_bend_allowance_m"] = 0.1
    assert validate_public_identity_v55(public)["status"] == "needs_attention"


def test_v55_self_consistent_dangling_pmi_or_open_brep_is_rejected():
    _, source = _payloads_v55(); source = deepcopy(source)
    source["replay_identity"][PMI]["product_association"] = source["replay_identity"][PMI]["replayed_product_association"] = {"feature:other": "part:housing"}
    source["replay_identity"][BREP_V55]["sewn_shells"][0]["closed"] = False
    source["replay_identity"][BREP_V55]["replayed_sewn_shells"] = source["replay_identity"][BREP_V55]["sewn_shells"]
    assert validate_source_identity_v55(source)["status"] == "needs_attention"


def test_v55_numeric_sha256_values_are_rejected():
    public, source = _payloads_v55()
    numeric_digest = int("9" * 64)
    rows = [*public["reference"], source["replay_identity"]]
    for container in rows:
        for row in container.values():
            for name in tuple(row):
                if name.endswith("_sha256"):
                    row[name] = numeric_digest
    assert validate_public_identity_v55(public)["status"] == "needs_attention"
    assert validate_source_identity_v55(source)["status"] == "needs_attention"


def test_v56_positive_public_and_source_identities_are_accepted() -> None:
    public, source = _payloads_v56()
    assert validate_public_identity_v56(public)["status"] == "ok"
    assert validate_source_identity_v56(source)["status"] == "ok"


def test_v56_frozen_mutations_are_rejected() -> None:
    public, source = _payloads_v56()
    public = deepcopy(public)
    source = deepcopy(source)
    public["reference"][0][GEAR]["result_pitch_diameter_m"] = 0.060
    public["reference"][0][PIPE]["result_self_intersection"] = True
    source["replay_identity"][STEP_V56]["replayed_length_unit"] = "inch"
    source["replay_identity"][MESH]["replayed_signed_volume_m3"] = -4.8e-5
    assert validate_public_identity_v56(public)["status"] == "needs_attention"
    assert validate_source_identity_v56(source)["status"] == "needs_attention"


def test_v56_self_consistent_geometry_contradictions_are_rejected() -> None:
    public, _ = _payloads_v56()
    public = deepcopy(public)
    gear = public["reference"][0][GEAR]
    gear["pitch_diameter_m"] = gear["result_pitch_diameter_m"] = 0.060
    pipe = public["reference"][0][PIPE]
    pipe["self_intersection"] = pipe["result_self_intersection"] = True
    assert validate_public_identity_v56(public)["status"] == "needs_attention"


def test_v56_self_consistent_exchange_contradictions_are_rejected() -> None:
    _, source = _payloads_v56()
    source = deepcopy(source)
    step = source["replay_identity"][STEP_V56]
    step["occurrence_transform_4x4"] = step["replayed_occurrence_transform_4x4"] = [[1.0]]
    mesh = source["replay_identity"][MESH]
    mesh["unit_scale_to_m"] = mesh["replayed_unit_scale_to_m"] = 0.0254
    assert validate_source_identity_v56(source)["status"] == "needs_attention"


def test_v56_numeric_digests_are_rejected() -> None:
    public, source = _payloads_v56()
    numeric_digest = int("1" * 64)
    for identity in (GEAR, PIPE):
        contract = public["reference"][0][identity]
        contract["result_sha256"] = numeric_digest
        contract["accepted_result_sha256"] = numeric_digest
    for identity in (STEP_V56, MESH):
        contract = source["replay_identity"][identity]
        contract["result_sha256"] = numeric_digest
        contract["accepted_result_sha256"] = numeric_digest
    assert validate_public_identity_v56(public)["status"] == "needs_attention"
    assert validate_source_identity_v56(source)["status"] == "needs_attention"


def test_v56_malformed_public_rows_are_not_silently_ignored() -> None:
    public, _ = _payloads_v56()
    public["reference"].append(42)
    assert validate_public_identity_v56(public)["status"] == "needs_attention"


def test_v56_unhashable_exchange_enums_reject_without_raising() -> None:
    _, source = _payloads_v56()
    step = source["replay_identity"][STEP_V56]
    step["length_unit"] = step["replayed_length_unit"] = []
    mesh = source["replay_identity"][MESH]
    mesh["mesh_format"] = mesh["replayed_mesh_format"] = []
    assert validate_source_identity_v56(source)["status"] == "needs_attention"
