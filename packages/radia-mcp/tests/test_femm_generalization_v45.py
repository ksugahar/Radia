from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.femm_v44_identity import validate_public_identity


_PROMOTED_CASE_IDS = (
    "v45_public_axisymmetric_force_coenergy_derivative_gap_coordinate_axisfactor_mesh_owner_mismatch",
    "v45_public_electrostatic_fringe_charge_energy_capacitance_interface_flux_axisfactor_mismatch",
)


def _identity():
    generation = "test-845"
    return {
        "v45_public_axisymmetric_force_torque_coenergy_stress_contour_owner_mismatch": {
            "generation": generation,
            **{key: generation for key in ("force_generation", "torque_generation", "coenergy_generation", "stress_contour_generation", "axis_factor_generation", "mesh_generation", "result_generation")},
            "force_method": "weighted_stress_tensor", "result_force_method": "weighted_stress_tensor", "torque_method": "airgap_contour", "result_torque_method": "airgap_contour",
            "force_n": [12.0, -3.0], "result_force_n": [12.0, -3.0], "torque_nm": 0.8, "result_torque_nm": 0.8, "coenergy_j": 2.1, "result_coenergy_j": 2.1,
            "axisymmetric_factor": 2.0 * math.pi, "result_axisymmetric_factor": 2.0 * math.pi, "contour_owner": "contour:test-845", "result_contour_owner": "contour:test-845", "mesh_owner": "mesh:test-845", "result_mesh_owner": "mesh:test-845",
            "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v45_public_electrostatic_fringe_charge_energy_capacitance_interface_flux_axisfactor_mismatch": {
            "generation": generation,
            **{key: generation for key in ("charge_generation", "energy_generation", "capacitance_generation", "interface_generation", "axis_factor_generation", "mesh_generation", "result_generation")},
            "voltage_v": 100.0, "result_voltage_v": 100.0, "capacitance_f": 1.0e-9, "result_capacitance_f": 1.0e-9, "charge_c": 1.0e-7, "result_charge_c": 1.0e-7,
            "stored_energy_j": 5.0e-6, "result_stored_energy_j": 5.0e-6, "interface_flux_c": 1.0e-7, "result_interface_flux_c": 1.0e-7, "axisymmetric_factor": 2.0 * math.pi, "result_axisymmetric_factor": 2.0 * math.pi,
            "mesh_owner": "mesh:test-845", "result_mesh_owner": "mesh:test-845", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


def test_v45_public_femm_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_v45_public_femm_identity_rejects_force_owner_mutation():
    identity = _identity()
    identity["v45_public_axisymmetric_force_torque_coenergy_stress_contour_owner_mismatch"]["result_mesh_owner"] = "stale:mesh"
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())


def test_v45_public_femm_identity_rejects_fringe_energy_mutation():
    identity = _identity()
    identity["v45_public_electrostatic_fringe_charge_energy_capacitance_interface_flux_axisfactor_mismatch"]["result_stored_energy_j"] = -1.0
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())
