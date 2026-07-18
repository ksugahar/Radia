from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v30 import _v30

_PROMOTED_CASE_IDS = (
    "v31_public_behavioral_source_event_discontinuity_timestep_derivative_charge_energy_mismatch",
    "v31_public_touchstone_reference_impedance_frequency_unit_port_order_passivity_mismatch",
)


def _v31():
    payload = _v30()
    positive = payload["metrics"]["positive"]
    generation = "behavioral-event-181"
    positive[
        "behavioral_source_event_timestep_derivative_charge_energy_initial_owner_result_identity"
    ] = {
        "behavioral_generation_id": generation,
        "event_behavioral_generation_id": generation,
        "timestep_behavioral_generation_id": generation,
        "derivative_behavioral_generation_id": generation,
        "charge_behavioral_generation_id": generation,
        "energy_behavioral_generation_id": generation,
        "initial_behavioral_generation_id": generation,
        "owner_behavioral_generation_id": generation,
        "result_behavioral_generation_id": generation,
        "event_time_s": 1.0e-6,
        "result_event_time_s": 1.0e-6,
        "time_grid_s": [0.0, 0.5e-6, 1.0e-6, 1.5e-6, 2.0e-6],
        "result_time_grid_s": [0.0, 0.5e-6, 1.0e-6, 1.5e-6, 2.0e-6],
        "derivative_convention": "right_limit_after_event",
        "result_derivative_convention": "right_limit_after_event",
        "initial_state": {"voltage_v": 0.0, "current_a": 0.0},
        "result_initial_state": {"voltage_v": 0.0, "current_a": 0.0},
        "integrated_charge_c": 1.0e-6,
        "result_integrated_charge_c": 1.0e-6,
        "source_energy_j": 5.0e-6,
        "result_source_energy_j": 5.0e-6,
        "waveform_owner_sha256": "1" * 64,
        "result_waveform_owner_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "touchstone-network-181"
    positive[
        "touchstone_impedance_frequency_parameter_port_complex_passivity_file_result_identity"
    ] = {
        "touchstone_generation_id": generation,
        "impedance_touchstone_generation_id": generation,
        "frequency_touchstone_generation_id": generation,
        "parameter_touchstone_generation_id": generation,
        "port_touchstone_generation_id": generation,
        "complex_touchstone_generation_id": generation,
        "passivity_touchstone_generation_id": generation,
        "file_touchstone_generation_id": generation,
        "result_touchstone_generation_id": generation,
        "reference_impedance_ohm": 50.0,
        "result_reference_impedance_ohm": 50.0,
        "frequency_unit": "GHz",
        "result_frequency_unit": "GHz",
        "parameter_type": "S",
        "result_parameter_type": "S",
        "port_order": [1, 2],
        "result_port_order": [1, 2],
        "complex_format": "RI",
        "result_complex_format": "RI",
        "maximum_singular_value": 0.98,
        "result_maximum_singular_value": 0.98,
        "passivity_tolerance": 1.0e-9,
        "touchstone_file_sha256": "3" * 64,
        "parsed_touchstone_file_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v31_public_positive_behavioral_source_and_touchstone_identities():
    assert ideal_transformer_identity_gate(_v31())["status"] == "ok"


def test_v31_public_behavioral_source_event_discontinuity_timestep_derivative_charge_energy_mismatch():
    payload = _v31()
    identity = payload["metrics"]["positive"][
        "behavioral_source_event_timestep_derivative_charge_energy_initial_owner_result_identity"
    ]
    identity.update(
        {
            "event_behavioral_generation_id": "behavioral-event-180",
            "energy_behavioral_generation_id": "behavioral-event-179",
            "result_event_time_s": 1.2e-6,
            "result_time_grid_s": [0.0, 0.7e-6, 1.4e-6, 2.0e-6],
            "result_derivative_convention": "centered_across_jump",
            "result_initial_state": {"voltage_v": 5.0, "current_a": 1.0},
            "result_integrated_charge_c": -1.0e-6,
            "result_source_energy_j": 8.0e-6,
            "result_waveform_owner_sha256": "9" * 64,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "behavioral_sources_use_current_event_grid_derivative_charge_energy_initial_owner_and_result"
    ]


def test_v31_public_touchstone_reference_impedance_frequency_unit_port_order_passivity_mismatch():
    payload = _v31()
    identity = payload["metrics"]["positive"][
        "touchstone_impedance_frequency_parameter_port_complex_passivity_file_result_identity"
    ]
    identity.update(
        {
            "impedance_touchstone_generation_id": "touchstone-network-180",
            "passivity_touchstone_generation_id": "touchstone-network-179",
            "result_reference_impedance_ohm": 75.0,
            "result_frequency_unit": "MHz",
            "result_parameter_type": "Y",
            "result_port_order": [2, 1],
            "result_complex_format": "DB",
            "result_maximum_singular_value": 1.2,
            "parsed_touchstone_file_sha256": "b" * 64,
            "accepted_result_sha256": "c" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "touchstone_networks_use_current_impedance_units_parameters_ports_complex_passivity_file_and_result"
    ]
