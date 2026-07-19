from __future__ import annotations

import math

from test_magnetic_force_generalization_v42 import _identity_v42
from test_magnetic_force_generalization_v31 import _gate


_SOLENOID = "axisymmetric_solenoid_flux_inductance_force_coenergy_axisfactor_mesh_generation_identity"
_DIELECTRIC = "dielectric_interface_capacitance_charge_flux_energy_reciprocity_mesh_generation_identity"
_PROMOTED_CASE_IDS = (
    "v43_public_axisymmetric_solenoid_flux_inductance_force_coenergy_axisfactor_mesh_mismatch",
    "v43_public_dielectric_interface_capacitance_charge_flux_energy_reciprocity_mesh_mismatch",
)


def _identity_v43():
    identity = _identity_v42()
    generation = "axisym-solenoid-843"
    identity[_SOLENOID] = {
        "solenoid_generation": generation,
        **{key: generation for key in ("flux_generation", "inductance_generation", "force_generation", "coenergy_generation", "axis_factor_generation", "mesh_generation", "result_generation")},
        "current_a": 3.0, "result_current_a": 3.0, "mean_radius_m": 0.02, "result_mean_radius_m": 0.02,
        "flux_linkage_wb_turn": 0.12, "result_flux_linkage_wb_turn": 0.12, "inductance_h": 0.04, "result_inductance_h": 0.04,
        "axial_force_n": 12.0, "result_axial_force_n": 12.0, "coenergy_force_derivative_n": 12.0, "result_coenergy_force_derivative_n": 12.0,
        "two_pi_r_factor_m": 2.0 * math.pi * 0.02, "result_two_pi_r_factor_m": 2.0 * math.pi * 0.02,
        "mesh_owner": "mesh:axisym-solenoid-843", "result_mesh_owner": "mesh:axisym-solenoid-843",
        "solenoid_result_sha256": "5" * 64, "accepted_solenoid_result_sha256": "5" * 64,
    }
    generation = "dielectric-interface-843"
    identity[_DIELECTRIC] = {
        "dielectric_generation": generation,
        **{key: generation for key in ("capacitance_generation", "charge_generation", "flux_generation", "energy_generation", "interface_generation", "mesh_generation", "result_generation")},
        "voltage_v": 100.0, "result_voltage_v": 100.0, "capacitance_f": 1.0e-9, "result_capacitance_f": 1.0e-9,
        "conductor_charge_c": 1.0e-7, "result_conductor_charge_c": 1.0e-7, "normal_displacement_flux_c": 1.0e-7, "result_normal_displacement_flux_c": 1.0e-7,
        "stored_energy_j": 5.0e-6, "result_stored_energy_j": 5.0e-6, "interface_continuity_residual_c": 1.0e-12, "result_interface_continuity_residual_c": 1.0e-12,
        "reciprocity_residual_f": 1.0e-12, "result_reciprocity_residual_f": 1.0e-12,
        "mesh_owner": "mesh:dielectric-interface-843", "result_mesh_owner": "mesh:dielectric-interface-843",
        "dielectric_result_sha256": "6" * 64, "accepted_dielectric_result_sha256": "6" * 64,
    }
    return identity


def test_v43_public_positive_contracts():
    assert _gate(_identity_v43())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v43_public_solenoid_mismatch():
    identity = _identity_v43()
    identity[_SOLENOID]["result_inductance_h"] = 0.08
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["axisymmetric_solenoids_close_flux_inductance_coenergy_force_two_pi_r_mesh_and_result"]


def test_v43_public_dielectric_mismatch():
    identity = _identity_v43()
    identity[_DIELECTRIC]["result_stored_energy_j"] = -1.0
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["dielectric_interfaces_close_capacitance_charge_flux_energy_reciprocity_mesh_and_result"]


def test_v43_public_rejects_self_consistent_wrong_two_pi_r_factor():
    identity = _identity_v43()
    record = identity[_SOLENOID]
    record["two_pi_r_factor_m"] = record["result_two_pi_r_factor_m"] = 1.0
    assert _gate(identity)["status"] == "needs_attention"
