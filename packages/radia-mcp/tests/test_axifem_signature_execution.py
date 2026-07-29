from __future__ import annotations

import copy
import hashlib
import json

from radia_mcp.fem.axifem_signature_execution import (
    validate_axifem_signature_execution,
)


FEATURES = {
    "arc_geometry": True,
    "circuits": False,
    "external_region": False,
    "moving_band": False,
    "nonlinear_bh": False,
    "periodic_boundary": False,
    "point_properties": False,
}


def _sha(value: dict) -> str:
    payload = dict(value)
    payload.pop("evidence_payload_sha256", None)
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _packet() -> dict:
    rows = []
    for index in range(27):
        rows.append(
            {
                "case_id": f"axisymmetric_dc_case_{index + 1:02d}",
                "private_record_sha256": f"{index + 1:064x}",
                "solver_lane": "axifem_henrotte",
                "formulation": "axisymmetric_magnetostatic_Aphi",
                "requested_features": dict(FEATURES),
                "mapped_features": dict(FEATURES),
                "mesh": {
                    "format": "netgen_vol",
                    "sha256": f"{index + 101:064x}",
                    "element_count": 128,
                    "vertex_count": 81,
                    "polygon_count": 2,
                },
                "display_mesh": {
                    "format": "gmsh_msh_4_1",
                    "sha256": f"{index + 201:064x}",
                },
                "result": {
                    "execution_status": "passed",
                    "solution_norm": 1.0,
                    "field_l2_sq_t2_m3": 0.1,
                    "relative_algebraic_residual": 1.0e-12,
                    "residual_limit": 1.0e-7,
                    "finite_nonzero_solution": True,
                    "nonlinear_executed": False,
                },
                "source_model_fidelity": {
                    "geometry": "source_profile_polygonized",
                    "material": "source_region_materials",
                    "stimulus": "source-faithful",
                    "boundary_operator": "source_homogeneous_dirichlet",
                    "required_features_mapped": True,
                },
                "point_source_evidence": {
                    "method": "vertex_dirac_ring_current",
                    "source_count": 0,
                    "embedded_vertex_count": 0,
                    "axis_annihilated_count": 0,
                    "max_vertex_distance_m": 0.0,
                    "max_weak_load_identity_abs_error_a_m": 0.0,
                    "all_point_properties_exact": False,
                },
                "point_potential_evidence": {
                    "method": "vertex_essential_aphi",
                    "constraint_count": 0,
                    "embedded_vertex_count": 0,
                    "max_vertex_distance_m": 0.0,
                    "max_constraint_abs_error_wb_per_m": 0.0,
                },
                "bh_curve_evidence": [],
                "boundary_operator_evidence": {
                    "source_boundary_object_count": 1,
                    "mixed_boundaries": [],
                    "periodic_pairs": [],
                    "constraint_residual": {
                        "constraint_count": 0,
                        "constraint_component_count": 0,
                        "point_constraint_count": 0,
                        "reduced_dof_count": 81,
                        "max_identification_abs_error": 0.0,
                        "max_point_constraint_abs_error": 0.0,
                        "relative_reduced_residual": 1.0e-12,
                    },
                    "dual_boundary": None,
                    "mapped": True,
                },
                "external_region_evidence": {
                    "mapped": True,
                    "region_count": 0,
                    "parameters_m": {},
                    "sample_factors": [],
                },
                "readiness_class": "source_faithful",
                "migration_disposition": "migrate_now",
                "migration_blockers": [],
                "timing_seconds": {
                    "parse_and_semantics": 0.1,
                    "mesh_and_vol_roundtrip": 0.2,
                    "solve_and_verify": 0.3,
                    "total": 0.7,
                },
            }
        )
    packet = {
        "schema": "radia.axifem-signature-execution.v1",
        "executed_at_utc": "2026-07-29T12:00:00+00:00",
        "execution_version": {
            "producer": "fixture",
            "producer_version": "24",
            "radia_mcp": "1.4.21",
            "ngsolve": "6.2.2604",
            "gmsh": "4.14.0",
        },
        "case_count": 27,
        "execution_passed_count": 27,
        "source_faithful_solver_ready_count": 27,
        "validation_surrogate_count": 0,
        "records": rows,
        "source_solver_launched": False,
        "retirement_ready_from_this_manifest_alone": False,
        "axisymmetric_dc_lane_retirement_ready": True,
        "retirement_ready": True,
        "classification_complete": True,
    }
    packet["evidence_payload_sha256"] = _sha(packet)
    return packet


def test_accepts_only_the_complete_27_case_retirement_packet() -> None:
    faithful = validate_axifem_signature_execution(_packet())

    assert faithful["status"] == "accepted"
    assert faithful["retirement_ready"] is True

    incomplete = _packet()
    incomplete["records"].pop()
    incomplete["case_count"] = 26
    incomplete["execution_passed_count"] = 26
    incomplete["source_faithful_solver_ready_count"] = 26
    incomplete["retirement_ready"] = False
    incomplete["evidence_payload_sha256"] = _sha(incomplete)
    result = validate_axifem_signature_execution(incomplete)
    assert result["status"] == "rejected"
    assert result["checks"]["all_27_signature_representatives_present"] is False


def test_rejects_false_faithful_claim_path_leak_and_residual_drift() -> None:
    packet = _packet()
    packet["records"][0]["source_model_fidelity"]["stimulus"] = (
        "manufactured-current-density"
    )
    packet["records"][0]["result"]["relative_algebraic_residual"] = 1.0
    packet["debug_path"] = r"\\private-host\private-share\case.vol"
    packet["evidence_payload_sha256"] = _sha(packet)

    result = validate_axifem_signature_execution(packet)

    assert result["status"] == "rejected"
    assert result["checks"]["public_boundary_has_no_paths"] is False
    assert result["checks"]["record_contracts"] is False


def test_rejects_periodic_constraint_count_laundering() -> None:
    packet = _packet()
    row = packet["records"][0]
    row["requested_features"]["periodic_boundary"] = True
    row["mapped_features"]["periodic_boundary"] = True
    row["source_model_fidelity"]["boundary_operator"] = (
        "source_homogeneous_dirichlet+signed_periodic_trace"
    )
    row["boundary_operator_evidence"]["source_boundary_object_count"] = 2
    row["boundary_operator_evidence"]["periodic_pairs"] = [
        {
            "boundary_property_index": 1,
            "trace_kind": "arc",
            "phase": -1.0,
            "source_object_count": 2,
            "vertex_pair_count": 2,
            "edge_pair_count": 1,
            "constraint_count": 3,
            "coordinate_tolerance_m": 1.0e-8,
        }
    ]
    row["boundary_operator_evidence"]["constraint_residual"][
        "constraint_count"
    ] = 2
    row["boundary_operator_evidence"]["constraint_residual"][
        "constraint_component_count"
    ] = 2
    packet["evidence_payload_sha256"] = _sha(packet)

    result = validate_axifem_signature_execution(packet)
    assert result["status"] == "rejected"
    assert result["checks"]["record_contracts"] is False


def test_rejects_external_region_identity_laundering() -> None:
    packet = _packet()
    row = packet["records"][0]
    row["requested_features"]["external_region"] = True
    row["mapped_features"]["external_region"] = True
    row["external_region_evidence"] = {
        "mapped": True,
        "region_count": 1,
        "parameters_m": {
            "z0": 0.0,
            "outer_radius": 0.1,
            "inner_radius": -0.2,
        },
        "sample_factors": [1.25],
        "coefficient_identity": "nu_external=nu_material*(r2+(z-z0)2)*Ri/Ro3",
    }
    packet["evidence_payload_sha256"] = _sha(packet)

    result = validate_axifem_signature_execution(packet)
    assert result["status"] == "rejected"
    assert result["checks"]["record_contracts"] is False


def test_rejects_dual_average_identity_drift() -> None:
    packet = _packet()
    row = packet["records"][0]
    base = row["boundary_operator_evidence"]["constraint_residual"]
    row["source_model_fidelity"]["boundary_operator"] = (
        "source_dual_boundary_average"
    )
    row["boundary_operator_evidence"]["dual_boundary"] = "source_dual"
    row["boundary_operator_evidence"]["constraint_residual"] = {
        "method": "dual_component_residuals",
        "natural": copy.deepcopy(base),
        "essential": copy.deepcopy(base),
        "average_identity_max_abs_error": 1.0e-6,
    }
    packet["evidence_payload_sha256"] = _sha(packet)

    result = validate_axifem_signature_execution(packet)
    assert result["status"] == "rejected"
    assert result["checks"]["record_contracts"] is False


def test_fem_server_registers_axifem_signature_execution_gate() -> None:
    from radia_mcp.fem.server import mcp

    names = {tool.name for tool in mcp._tool_manager.list_tools()}
    assert "fem_axifem_signature_execution_gate" in names
