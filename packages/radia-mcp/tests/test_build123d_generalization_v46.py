from __future__ import annotations

from copy import deepcopy

from radia_mcp.build123d.build123d_v46_identity import (
    validate_public_identity,
    validate_source_identity,
)


PLACEMENT = "placement_unit_scale_inertia_coordinate_frame_nonfinite_identity"
BOOLEAN = "boolean_partial_shape_fillet_recovery_topology_identity"
SKETCH = "sketch_solver_partial_constraint_warning_plane_unit_identity"
STEP = "step_import_partial_face_tolerance_coordinate_frame_checksum_identity"


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
    row = {PLACEMENT: placement, BOOLEAN: boolean}
    return {"reference": [row], "measured": {"headless": [deepcopy(row)]}}


def _source_payload() -> dict[str, object]:
    return {
        "replay_identity": {
            SKETCH: {
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
            STEP: {
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


def test_v46_public_and_source_positive_identity_controls():
    assert validate_public_identity(_public_payload())["status"] == "ok"
    assert validate_source_identity(_source_payload())["status"] == "ok"


def test_v46_public_and_source_mutations_are_rejected():
    public = _public_payload()
    public["reference"][0][PLACEMENT]["result_unit_scale_to_si"] = float("nan")
    assert validate_public_identity(public)["status"] == "needs_attention"

    public = _public_payload()
    public["reference"][0][BOOLEAN]["result_partial_shape_status"] = "partial"
    assert validate_public_identity(public)["status"] == "needs_attention"

    source = _source_payload()
    source["replay_identity"][SKETCH]["result_plane"] = "YZ"
    assert validate_source_identity(source)["status"] == "needs_attention"

    source = _source_payload()
    source["replay_identity"][STEP]["result_checksum_sha256"] = "6" * 64
    assert validate_source_identity(source)["status"] == "needs_attention"


def test_v46_case_ids_are_frozen():
    case_ids = (
        "v46_public_placement_unit_scale_inertia_coordinate_frame_nonfinite_mismatch",
        "v46_public_boolean_partial_shape_fillet_failure_recovery_topology_mismatch",
        "v46_source_tool_sketch_solver_partial_constraint_warning_plane_unit_mismatch",
        "v46_source_tool_step_import_partial_face_tolerance_coordinate_frame_checksum_mismatch",
    )
    assert len(case_ids) == 4
    assert all(case_id.startswith("v46_") for case_id in case_ids)
