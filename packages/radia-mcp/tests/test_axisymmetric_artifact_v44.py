from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.femm_v44_identity import validate_public_identity


_PROMOTED_CASE_IDS = (
    "v44_public_axisymmetric_actuator_flux_force_coenergy_gap_derivative_axisfactor_mesh_mismatch",
    "v44_public_electrostatic_fringe_capacitance_charge_energy_interface_axisfactor_mesh_mismatch",
)


def _identity():
    generation = "g-844"
    fields = {
        "generation": generation,
        "flux_generation": generation,
        "force_generation": generation,
        "coenergy_generation": generation,
        "gap_generation": generation,
        "derivative_generation": generation,
        "axis_factor_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "flux_linkage_wb_turn": 0.12,
        "result_flux_linkage_wb_turn": 0.12,
        "force_n": 12.0,
        "result_force_n": 12.0,
        "coenergy_j": 2.1,
        "result_coenergy_j": 2.1,
        "gap_m": 0.002,
        "result_gap_m": 0.002,
        "coenergy_force_derivative_n": 12.0,
        "result_coenergy_force_derivative_n": 12.0,
        "axisymmetric_factor": 2.0 * math.pi,
        "result_axisymmetric_factor": 2.0 * math.pi,
        "mesh_owner": "mesh:g-844",
        "result_mesh_owner": "mesh:g-844",
        "result_sha256": "a" * 64,
        "accepted_result_sha256": "a" * 64,
    }
    return {"v44_axisymmetric_actuator_flux_force_coenergy_gap_derivative_axisfactor_mesh_mismatch": fields}


def test_axisymmetric_actuator_identity_accepts_closed_artifact():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_axisymmetric_actuator_identity_rejects_stale_result():
    identity = _identity()
    identity[next(iter(identity))]["result_force_n"] = -12.0
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())
