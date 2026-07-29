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


def _packet(*, faithful: bool = True) -> dict:
    mapped = dict(FEATURES)
    fidelity = {
        "geometry": "source_profile_polygonized",
        "material": "source_region_materials",
        "stimulus": "source-faithful" if faithful else "manufactured-current-density",
        "boundary_operator": "source_homogeneous_dirichlet"
        if faithful
        else "homogeneous_outer_truncation",
        "required_features_mapped": faithful,
    }
    row = {
        "case_id": "axisymmetric_dc_case_01",
        "private_record_sha256": "a" * 64,
        "solver_lane": "axifem_henrotte",
        "formulation": "axisymmetric_magnetostatic_Aphi",
        "requested_features": dict(FEATURES),
        "mapped_features": mapped,
        "mesh": {
            "format": "netgen_vol",
            "sha256": "b" * 64,
            "element_count": 128,
            "vertex_count": 81,
            "polygon_count": 2,
        },
        "display_mesh": {"format": "gmsh_msh_4_1", "sha256": "c" * 64},
        "result": {
            "execution_status": "passed",
            "solution_norm": 1.0,
            "field_l2_sq_t2_m3": 0.1,
            "relative_algebraic_residual": 1.0e-12,
            "residual_limit": 1.0e-7,
            "finite_nonzero_solution": True,
            "nonlinear_executed": False,
        },
        "source_model_fidelity": fidelity,
        "readiness_class": "source_faithful" if faithful else "validation_surrogate",
        "timing_seconds": {
            "parse_and_semantics": 0.1,
            "mesh_and_vol_roundtrip": 0.2,
            "solve_and_verify": 0.3,
            "total": 0.7,
        },
    }
    packet = {
        "schema": "radia.axifem-signature-execution.v1",
        "executed_at_utc": "2026-07-29T12:00:00+00:00",
        "execution_version": {
            "producer": "fixture",
            "producer_version": "20",
            "radia_mcp": "1.4.19",
            "ngsolve": "6.2.2604",
            "gmsh": "4.14.0",
        },
        "case_count": 1,
        "execution_passed_count": 1,
        "source_faithful_solver_ready_count": int(faithful),
        "validation_surrogate_count": int(not faithful),
        "records": [row],
        "source_solver_launched": False,
        "retirement_ready": faithful,
    }
    packet["evidence_payload_sha256"] = _sha(packet)
    return packet


def test_accepts_execution_and_keeps_surrogate_distinct_from_retirement() -> None:
    faithful = validate_axifem_signature_execution(_packet(faithful=True))
    surrogate = validate_axifem_signature_execution(_packet(faithful=False))

    assert faithful["status"] == "accepted"
    assert faithful["retirement_ready"] is True
    assert surrogate["status"] == "accepted"
    assert surrogate["retirement_ready"] is False
    assert surrogate["validation_surrogate_count"] == 1


def test_rejects_false_faithful_claim_path_leak_and_residual_drift() -> None:
    packet = _packet(faithful=False)
    packet["records"][0]["readiness_class"] = "source_faithful"
    packet["records"][0]["result"]["relative_algebraic_residual"] = 1.0
    packet["debug_path"] = r"\\private-host\private-share\case.vol"
    packet["evidence_payload_sha256"] = _sha(packet)

    result = validate_axifem_signature_execution(packet)

    assert result["status"] == "rejected"
    assert result["checks"]["public_boundary_has_no_paths"] is False
    assert result["checks"]["record_contracts"] is False


def test_fem_server_registers_axifem_signature_execution_gate() -> None:
    from radia_mcp.fem.server import mcp

    names = {tool.name for tool in mcp._tool_manager.list_tools()}
    assert "fem_axifem_signature_execution_gate" in names
