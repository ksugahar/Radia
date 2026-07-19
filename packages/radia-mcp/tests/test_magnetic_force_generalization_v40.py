from __future__ import annotations

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v39 import _identity_v39


_HEAT = "heat_conduction_contact_convection_flux_temperature_energy_mesh_result_generation_identity"
_CURRENT = "current_flow_electrode_voltage_current_resistance_joule_power_reciprocity_conductor_mesh_result_generation_identity"
_PROMOTED_CASE_IDS = (
    "v40_public_heat_conduction_contact_resistance_convection_flux_temperature_energy_mismatch",
    "v40_public_current_flow_electrode_voltage_current_resistance_joule_power_reciprocity_mismatch",
)


def _identity_v40():
    identity = _identity_v39()
    generation = "heat-contact-311"
    area = 1.0e-2
    length = 2.0e-2
    conductivity = 15.0
    contact_resistance_area = 1.0e-4
    convection = 25.0
    hot = 373.15
    ambient = 293.15
    contact_resistance = contact_resistance_area / area
    convection_resistance = 1.0 / (convection * area)
    heat_rate = (hot - ambient) / (
        length / (conductivity * area) + contact_resistance + convection_resistance
    )
    identity[_HEAT] = {
        "heat_generation": generation,
        **{
            key: generation
            for key in (
                "conductivity_generation", "contact_generation", "convection_generation",
                "flux_generation", "temperature_generation", "energy_generation",
                "mesh_generation", "result_generation",
            )
        },
        "conductivity_w_per_m_k": conductivity,
        "result_conductivity_w_per_m_k": conductivity,
        "conduction_length_m": length,
        "result_conduction_length_m": length,
        "boundary_area_m2": area,
        "result_boundary_area_m2": area,
        "contact_resistance_m2_k_per_w": contact_resistance_area,
        "result_contact_resistance_m2_k_per_w": contact_resistance_area,
        "convection_coefficient_w_per_m2_k": convection,
        "result_convection_coefficient_w_per_m2_k": convection,
        "hot_temperature_k": hot,
        "result_hot_temperature_k": hot,
        "ambient_temperature_k": ambient,
        "result_ambient_temperature_k": ambient,
        "boundary_heat_flux_w_per_m2": heat_rate / area,
        "result_boundary_heat_flux_w_per_m2": heat_rate / area,
        "interface_temperature_jump_k": heat_rate * contact_resistance,
        "result_interface_temperature_jump_k": heat_rate * contact_resistance,
        "convection_surface_temperature_k": ambient + heat_rate * convection_resistance,
        "result_convection_surface_temperature_k": ambient + heat_rate * convection_resistance,
        "total_heat_rate_w": heat_rate,
        "result_total_heat_rate_w": heat_rate,
        "energy_balance_residual_w": 0.0,
        "result_energy_balance_residual_w": 0.0,
        "mesh_owner": "heat:mesh-311",
        "accepted_mesh_owner": "heat:mesh-311",
        "heat_result_sha256": "5" * 64,
        "accepted_heat_result_sha256": "5" * 64,
    }
    generation = "current-flow-311"
    identity[_CURRENT] = {
        "current_generation": generation,
        **{
            key: generation
            for key in (
                "electrode_generation", "voltage_generation", "current_generation_id",
                "resistance_generation", "joule_generation", "power_generation",
                "reciprocity_generation", "conductor_generation", "mesh_generation",
                "result_generation",
            )
        },
        "electrode_voltage_v": [10.0, 0.0],
        "result_electrode_voltage_v": [10.0, 0.0],
        "terminal_current_a": [2.0, -2.0],
        "result_terminal_current_a": [2.0, -2.0],
        "effective_resistance_ohm": 5.0,
        "result_effective_resistance_ohm": 5.0,
        "joule_loss_w": 20.0,
        "result_joule_loss_w": 20.0,
        "terminal_power_w": 20.0,
        "result_terminal_power_w": 20.0,
        "reciprocity_residual_ohm": 0.0,
        "result_reciprocity_residual_ohm": 0.0,
        "conductor_owner": "current:conductors-311",
        "accepted_conductor_owner": "current:conductors-311",
        "mesh_owner": "current:mesh-311",
        "accepted_mesh_owner": "current:mesh-311",
        "current_result_sha256": "6" * 64,
        "accepted_current_result_sha256": "6" * 64,
    }
    return identity


def test_v40_public_positive_heat_and_current_flow_closure():
    assert _gate(_identity_v40())["status"] == "ok"


def test_v40_public_heat_contact_mismatch():
    identity = _identity_v40()
    identity[_HEAT].update(
        {
            "contact_generation": "heat-contact-310",
            "energy_generation": "heat-contact-309",
            "result_generation": "heat-contact-308",
            "result_conductivity_w_per_m_k": -15.0,
            "result_contact_resistance_m2_k_per_w": -1.0e-4,
            "result_convection_coefficient_w_per_m2_k": 0.0,
            "result_boundary_heat_flux_w_per_m2": -1.0,
            "result_interface_temperature_jump_k": -1.0,
            "result_convection_surface_temperature_k": 400.0,
            "result_total_heat_rate_w": -20.0,
            "result_energy_balance_residual_w": 5.0,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_heat_result_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "heat_contact_models_close_conduction_contact_convection_flux_temperature_energy_mesh_and_result"
    ]


def test_v40_public_current_flow_mismatch():
    identity = _identity_v40()
    identity[_CURRENT].update(
        {
            "electrode_generation": "current-flow-310",
            "power_generation": "current-flow-309",
            "result_generation": "current-flow-308",
            "result_electrode_voltage_v": [0.0, 10.0],
            "result_terminal_current_a": [2.0, 2.0],
            "result_effective_resistance_ohm": -5.0,
            "result_joule_loss_w": -20.0,
            "result_terminal_power_w": 0.0,
            "result_reciprocity_residual_ohm": 1.0,
            "accepted_conductor_owner": "stale:conductor",
            "accepted_mesh_owner": "stale:mesh",
            "accepted_current_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "current_flow_models_close_electrodes_currents_resistance_joule_terminal_power_reciprocity_owners_and_result"
    ]


def test_v40_public_rejects_self_consistent_wrong_contact_jump():
    identity = _identity_v40()
    identity[_HEAT]["interface_temperature_jump_k"] = 1.0
    identity[_HEAT]["result_interface_temperature_jump_k"] = 1.0
    assert _gate(identity)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_current_nonconservation():
    identity = _identity_v40()
    identity[_CURRENT]["terminal_current_a"] = [2.0, -1.5]
    identity[_CURRENT]["result_terminal_current_a"] = [2.0, -1.5]
    assert _gate(identity)["status"] == "needs_attention"
