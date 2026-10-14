from __future__ import annotations

from radia_mcp.radia_ngsolve.motor_v44_identity import validate_public_identity


_PROMOTED_CASE_IDS = (
    "v44_public_ipmsm_torque_ripple_radialforce_modal_power_efficiency_energy_mesh_mismatch",
    "v44_public_inductionmotor_slip_rotorloss_torque_current_power_heat_energy_mismatch",
)


def _identity():
    return {
        "v44_ipmsm_torque_ripple_radialforce_modal_power_efficiency_energy_mesh_mismatch": {
            "generation": "g", "torque_generation": "g", "radial_force_generation": "g", "modal_generation": "g", "power_generation": "g", "efficiency_generation": "g", "energy_generation": "g", "mesh_generation": "g", "result_generation": "g",
            "torque_ripple_rms_nm": 0.8, "result_torque_ripple_rms_nm": 0.8, "radial_force_space_orders": [6, 12], "result_radial_force_space_orders": [6, 12], "modal_excitation_n": 0.3, "result_modal_excitation_n": 0.3, "electromagnetic_power_w": 1000.0, "result_electromagnetic_power_w": 1000.0, "mechanical_power_w": 950.0, "result_mechanical_power_w": 950.0, "efficiency": 0.95, "result_efficiency": 0.95, "energy_closure_residual": 1.0e-10, "result_energy_closure_residual": 1.0e-10, "mesh_owner": "mesh:g", "result_mesh_owner": "mesh:g", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v44_inductionmotor_slip_rotorloss_torque_current_power_heat_energy_mismatch": {
            "generation": "h", "slip_generation": "h", "rotor_loss_generation": "h", "torque_generation": "h", "current_generation": "h", "power_generation": "h", "heat_generation": "h", "energy_generation": "h", "mesh_generation": "h", "result_generation": "h",
            "slip": 0.05, "result_slip": 0.05, "rotor_copper_loss_w": 50.0, "result_rotor_copper_loss_w": 50.0, "torque_nm": 20.0, "result_torque_nm": 20.0, "stator_current_a": 10.0, "result_stator_current_a": 10.0, "electrical_power_w": 1000.0, "result_electrical_power_w": 1000.0, "mechanical_power_w": 900.0, "result_mechanical_power_w": 900.0, "heat_loss_w": 50.0, "result_heat_loss_w": 50.0, "energy_closure_residual": 1.0e-10, "result_energy_closure_residual": 1.0e-10, "mesh_owner": "mesh:h", "result_mesh_owner": "mesh:h", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


def test_motor_v44_public_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_motor_v44_public_identity_rejects_efficiency_mismatch():
    identity = _identity()
    identity[next(iter(identity))]["result_efficiency"] = 1.1
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())
