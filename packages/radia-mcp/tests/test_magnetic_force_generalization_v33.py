from __future__ import annotations

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v32 import _identity_v32


_PROMOTED_CASE_IDS = (
    "v33_public_electrostatic_capacitance_matrix_reciprocity_charge_neutrality_energy_mismatch",
    "v33_public_axisymmetric_heat_flux_conduction_convection_source_2pir_balance_mismatch",
)


def _identity_v33():
    identity = _identity_v32()
    generation = "electrostatic-capacitance-201"
    matrix = [
        [2e-12, -1e-12, -1e-12],
        [-1e-12, 2e-12, -1e-12],
        [-1e-12, -1e-12, 2e-12],
    ]
    identity[
        "electrostatic_capacitance_conductor_reciprocity_neutrality_voltage_charge_energy_unit_mesh_owner_result_identity"
    ] = {
        "electrostatic_generation": generation,
        **{
            key: generation
            for key in (
                "conductor_generation",
                "reciprocity_generation",
                "neutrality_generation",
                "voltage_generation",
                "charge_generation",
                "energy_generation",
                "unit_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "conductor_order": ["left", "shield", "right"],
        "result_conductor_order": ["left", "shield", "right"],
        "capacitance_matrix_f": matrix,
        "result_capacitance_matrix_f": [row[:] for row in matrix],
        "voltages_v": [1.0, 0.0, -1.0],
        "result_voltages_v": [1.0, 0.0, -1.0],
        "charges_c": [3e-12, 0.0, -3e-12],
        "result_charges_c": [3e-12, 0.0, -3e-12],
        "electrostatic_energy_j": 3e-12,
        "result_electrostatic_energy_j": 3e-12,
        "capacitance_unit": "F",
        "result_capacitance_unit": "F",
        "charge_unit": "C",
        "result_charge_unit": "C",
        "voltage_unit": "V",
        "result_voltage_unit": "V",
        "energy_unit": "J",
        "result_energy_unit": "J",
        "electrostatic_mesh_sha256": "1" * 64,
        "result_electrostatic_mesh_sha256": "1" * 64,
        "result_owner": "electrostatic/capacitance-matrix-201",
        "accepted_result_owner": "electrostatic/capacitance-matrix-201",
        "electrostatic_result_sha256": "2" * 64,
        "accepted_electrostatic_result_sha256": "2" * 64,
    }
    generation = "axisymmetric-heat-balance-201"
    identity[
        "axisymmetric_heat_conduction_convection_source_boundary_flux_weight_temperature_mesh_owner_result_identity"
    ] = {
        "heat_generation": generation,
        **{
            key: generation
            for key in (
                "conduction_generation",
                "convection_generation",
                "source_generation",
                "boundary_generation",
                "weight_generation",
                "temperature_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "conduction_w": 40.0,
        "result_conduction_w": 40.0,
        "convection_w": 30.0,
        "result_convection_w": 30.0,
        "volume_source_w": 100.0,
        "result_volume_source_w": 100.0,
        "boundary_flux_w": 30.0,
        "result_boundary_flux_w": 30.0,
        "axisymmetric_weight": "2*pi*r",
        "result_axisymmetric_weight": "2*pi*r",
        "temperature_reference_k": 293.15,
        "result_temperature_reference_k": 293.15,
        "heat_balance_tolerance_w": 1e-9,
        "heat_mesh_sha256": "3" * 64,
        "result_heat_mesh_sha256": "3" * 64,
        "heat_result_owner": "axisymmetric/thermal-body-201",
        "result_heat_result_owner": "axisymmetric/thermal-body-201",
        "heat_result_sha256": "4" * 64,
        "accepted_heat_result_sha256": "4" * 64,
    }
    return identity


def test_v33_public_positive_electrostatic_and_heat_closure():
    assert _gate(_identity_v33())["status"] == "ok"


def test_v33_public_electrostatic_capacitance_matrix_reciprocity_charge_neutrality_energy_mismatch():
    identity = _identity_v33()
    identity[
        "electrostatic_capacitance_conductor_reciprocity_neutrality_voltage_charge_energy_unit_mesh_owner_result_identity"
    ].update(
        {
            "reciprocity_generation": "electrostatic-capacitance-200",
            "energy_generation": "electrostatic-capacitance-199",
            "result_generation": "electrostatic-capacitance-198",
            "result_conductor_order": ["right", "shield", "left"],
            "result_capacitance_matrix_f": [[2e-12, 1e-12, -1e-12], [-2e-12, 2e-12, 0.0], [-1e-12, -1e-12, 3e-12]],
            "result_voltages_v": [-1.0, 0.0, 1.0],
            "result_charges_c": [3e-12, 1e-12, -2e-12],
            "result_electrostatic_energy_j": -3e-12,
            "result_capacitance_unit": "pF",
            "result_charge_unit": "nC",
            "result_electrostatic_mesh_sha256": "9" * 64,
            "accepted_result_owner": "stale/result",
            "accepted_electrostatic_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrostatic_capacitance_closes_conductors_reciprocity_neutrality_charge_energy_units_mesh_owner_and_result"
    ]


def test_v33_public_axisymmetric_heat_flux_conduction_convection_source_2pir_balance_mismatch():
    identity = _identity_v33()
    identity[
        "axisymmetric_heat_conduction_convection_source_boundary_flux_weight_temperature_mesh_owner_result_identity"
    ].update(
        {
            "source_generation": "axisymmetric-heat-balance-200",
            "weight_generation": "axisymmetric-heat-balance-199",
            "result_generation": "axisymmetric-heat-balance-198",
            "result_conduction_w": 20.0,
            "result_convection_w": 10.0,
            "result_boundary_flux_w": 5.0,
            "result_axisymmetric_weight": "1",
            "result_temperature_reference_k": 20.0,
            "result_heat_mesh_sha256": "b" * 64,
            "result_heat_result_owner": "planar/old",
            "accepted_heat_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_heat_closes_conduction_convection_source_flux_weight_temperature_mesh_owner_and_result"
    ]
