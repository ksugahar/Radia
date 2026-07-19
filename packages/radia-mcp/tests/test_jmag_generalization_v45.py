from __future__ import annotations

from radia_mcp.radia_ngsolve.motor_v44_identity import validate_public_identity


_PROMOTED_CASE_IDS = (
    "v45_public_ipmsm_torque_ripple_radial_force_modal_power_efficiency_energy_mesh_owner_mismatch",
    "v45_public_induction_machine_slip_rotor_loss_torque_heat_power_energy_result_mismatch",
)


def _identity():
    generation = "test-845"
    return {
        "v45_public_ipmsm_torque_ripple_radial_force_modal_power_efficiency_energy_mesh_owner_mismatch": {
            "generation": generation, **{key: generation for key in ("torque_generation", "radial_force_generation", "modal_generation", "power_generation", "efficiency_generation", "energy_generation", "mesh_generation", "result_generation")},
            "torque_ripple_rms_nm": 0.8, "result_torque_ripple_rms_nm": 0.8, "radial_force_space_orders": [6, 12], "result_radial_force_space_orders": [6, 12], "modal_excitation_n": 0.3, "result_modal_excitation_n": 0.3, "electromagnetic_power_w": 1000.0, "result_electromagnetic_power_w": 1000.0, "mechanical_power_w": 950.0, "result_mechanical_power_w": 950.0, "efficiency": 0.95, "result_efficiency": 0.95, "energy_closure_residual": 1e-10, "result_energy_closure_residual": 1e-10, "mesh_owner": "mesh:test", "result_mesh_owner": "mesh:test", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v45_public_induction_machine_slip_rotor_loss_torque_heat_power_energy_result_mismatch": {
            "generation": generation, **{key: generation for key in ("slip_generation", "rotor_loss_generation", "torque_generation", "heat_generation", "power_generation", "energy_generation", "result_generation")},
            "slip": 0.05, "result_slip": 0.05, "rotor_copper_loss_w": 50.0, "result_rotor_copper_loss_w": 50.0, "torque_nm": 20.0, "result_torque_nm": 20.0, "heat_loss_w": 50.0, "result_heat_loss_w": 50.0, "electrical_power_w": 1000.0, "result_electrical_power_w": 1000.0, "mechanical_power_w": 900.0, "result_mechanical_power_w": 900.0, "energy_closure_residual": 1e-10, "result_energy_closure_residual": 1e-10, "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


def test_v45_jmag_public_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_v45_jmag_public_identity_rejects_torque_mutation():
    identity = _identity()
    identity["v45_public_ipmsm_torque_ripple_radial_force_modal_power_efficiency_energy_mesh_owner_mismatch"]["result_efficiency"] = 1.1
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())
