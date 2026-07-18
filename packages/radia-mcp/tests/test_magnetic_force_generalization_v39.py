from __future__ import annotations

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v38 import _identity_v38


_NONLINEAR = "nonlinear_magnetic_circuit_bh_flux_linkage_coenergy_incremental_inductance_force_power_mesh_result_generation_identity"
_ELECTROSTATIC = "electrostatic_capacitance_matrix_charge_energy_reciprocity_gauge_conductor_mesh_result_generation_identity"
_PROMOTED_CASE_IDS = (
    "v39_public_nonlinear_magnetic_circuit_bh_coenergy_incremental_inductance_force_mismatch",
    "v39_public_electrostatic_capacitance_matrix_charge_energy_reciprocity_gauge_mismatch",
)


def _identity_v39():
    identity = _identity_v38()
    generation = "nonlinear-magnetic-circuit-271"
    identity[_NONLINEAR] = {
        "magnetic_generation": generation,
        **{key: generation for key in ("bh_generation", "flux_generation", "coenergy_generation", "inductance_generation", "force_generation", "power_generation", "mesh_generation", "result_generation")},
        "bh_curve_h_a_per_m": [0.0, 100.0, 500.0, 1000.0],
        "result_bh_curve_h_a_per_m": [0.0, 100.0, 500.0, 1000.0],
        "bh_curve_b_t": [0.0, 0.5, 1.2, 1.45],
        "result_bh_curve_b_t": [0.0, 0.5, 1.2, 1.45],
        "operating_h_a_per_m": 500.0,
        "result_operating_h_a_per_m": 500.0,
        "operating_b_t": 1.2,
        "result_operating_b_t": 1.2,
        "terminal_current_a": 10.0,
        "result_terminal_current_a": 10.0,
        "flux_linkage_wb_turn": 2.0e-2,
        "result_flux_linkage_wb_turn": 2.0e-2,
        "magnetic_energy_j": 8.0e-2,
        "result_magnetic_energy_j": 8.0e-2,
        "coenergy_j": 1.2e-1,
        "result_coenergy_j": 1.2e-1,
        "current_increment_a": 5.0e-1,
        "result_current_increment_a": 5.0e-1,
        "flux_linkage_increment_wb_turn": 1.0e-3,
        "result_flux_linkage_increment_wb_turn": 1.0e-3,
        "incremental_inductance_h": 2.0e-3,
        "result_incremental_inductance_h": 2.0e-3,
        "virtual_displacement_m": 1.0e-3,
        "result_virtual_displacement_m": 1.0e-3,
        "coenergy_increment_j": 2.0e-2,
        "result_coenergy_increment_j": 2.0e-2,
        "virtual_work_force_n": 20.0,
        "result_virtual_work_force_n": 20.0,
        "terminal_voltage_v": 5.0,
        "result_terminal_voltage_v": 5.0,
        "terminal_power_w": 50.0,
        "result_terminal_power_w": 50.0,
        "mesh_owner": "magnetic:mesh-271",
        "accepted_mesh_owner": "magnetic:mesh-271",
        "nonlinear_result_sha256": "1" * 64,
        "accepted_nonlinear_result_sha256": "1" * 64,
    }
    generation = "electrostatic-capacitance-271"
    identity[_ELECTROSTATIC] = {
        "electrostatic_generation": generation,
        **{key: generation for key in ("matrix_generation", "charge_generation", "energy_generation", "reciprocity_generation", "gauge_generation", "conductor_generation", "mesh_generation", "result_generation")},
        "capacitance_matrix_f": [[2.0e-11, -2.0e-11], [-2.0e-11, 2.0e-11]],
        "result_capacitance_matrix_f": [[2.0e-11, -2.0e-11], [-2.0e-11, 2.0e-11]],
        "terminal_voltage_v": [100.0, 0.0],
        "result_terminal_voltage_v": [100.0, 0.0],
        "terminal_charge_c": [2.0e-9, -2.0e-9],
        "result_terminal_charge_c": [2.0e-9, -2.0e-9],
        "field_energy_j": 1.0e-7,
        "result_field_energy_j": 1.0e-7,
        "reciprocity_residual_f": 0.0,
        "result_reciprocity_residual_f": 0.0,
        "reference_gauge": "conductor:2=0V",
        "result_reference_gauge": "conductor:2=0V",
        "conductor_owner": "electrostatic:conductors-271",
        "accepted_conductor_owner": "electrostatic:conductors-271",
        "mesh_owner": "electrostatic:mesh-271",
        "accepted_mesh_owner": "electrostatic:mesh-271",
        "electrostatic_result_sha256": "2" * 64,
        "accepted_electrostatic_result_sha256": "2" * 64,
    }
    return identity


def test_v39_public_positive_nonlinear_and_electrostatic_closure():
    assert _gate(_identity_v39())["status"] == "ok"


def test_v39_public_nonlinear_magnetic_circuit_bh_coenergy_incremental_inductance_force_mismatch():
    identity = _identity_v39()
    identity[_NONLINEAR].update({"bh_generation": "nonlinear-magnetic-circuit-270", "force_generation": "nonlinear-magnetic-circuit-269", "result_generation": "nonlinear-magnetic-circuit-268", "result_operating_h_a_per_m": 1000.0, "result_operating_b_t": -1.2, "result_flux_linkage_wb_turn": -2.0e-2, "result_magnetic_energy_j": -8.0e-2, "result_coenergy_j": -1.2e-1, "result_current_increment_a": -5.0e-1, "result_flux_linkage_increment_wb_turn": -1.0e-3, "result_incremental_inductance_h": -2.0e-3, "result_coenergy_increment_j": -2.0e-2, "result_virtual_work_force_n": -20.0, "result_terminal_power_w": -50.0, "accepted_mesh_owner": "stale:mesh", "accepted_nonlinear_result_sha256": "a" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["nonlinear_magnetic_circuits_close_bh_flux_coenergy_incremental_inductance_force_power_mesh_and_result"]


def test_v39_public_electrostatic_capacitance_matrix_charge_energy_reciprocity_gauge_mismatch():
    identity = _identity_v39()
    identity[_ELECTROSTATIC].update({"matrix_generation": "electrostatic-capacitance-270", "gauge_generation": "electrostatic-capacitance-269", "result_generation": "electrostatic-capacitance-268", "result_capacitance_matrix_f": [[2.0e-11, 1.0e-11], [-2.0e-11, 2.0e-11]], "result_terminal_charge_c": [1.0e-9, 1.0e-9], "result_field_energy_j": -1.0e-7, "result_reciprocity_residual_f": 1.0e-11, "result_reference_gauge": "floating", "accepted_conductor_owner": "stale:conductors", "accepted_mesh_owner": "stale:mesh", "accepted_electrostatic_result_sha256": "b" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["electrostatic_capacitance_closes_matrix_charge_energy_reciprocity_gauge_owners_and_result"]


def test_v39_public_rejects_self_consistent_wrong_coenergy():
    identity = _identity_v39()
    identity[_NONLINEAR]["coenergy_j"] = 0.2
    identity[_NONLINEAR]["result_coenergy_j"] = 0.2
    assert _gate(identity)["status"] == "needs_attention"


def test_v39_public_accepts_nonlinear_energy_coenergy_partition():
    identity = _identity_v39()
    row = identity[_NONLINEAR]
    assert row["magnetic_energy_j"] != row["coenergy_j"]
    assert row["coenergy_j"] != 0.5 * row["terminal_current_a"] * row["flux_linkage_wb_turn"]
    assert _gate(identity)["status"] == "ok"


def test_v39_public_rejects_self_consistent_nonsymmetric_capacitance():
    identity = _identity_v39()
    matrix = [[2.0e-11, -1.0e-11], [-2.0e-11, 2.0e-11]]
    identity[_ELECTROSTATIC]["capacitance_matrix_f"] = matrix
    identity[_ELECTROSTATIC]["result_capacitance_matrix_f"] = matrix
    assert _gate(identity)["status"] == "needs_attention"
