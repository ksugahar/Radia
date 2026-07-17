from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v23 import _summary_v23


_PROMOTED_CASE_IDS = (
    "v24_public_deembedding_reference_plane_phase_causality_passivity_grid_mismatch",
    "v24_public_field_circuit_cosim_port_sign_impedance_power_balance_generation_mismatch",
)


def _summary_v24():
    summary = _summary_v23()
    for index, row in enumerate(summary["runs"]):
        deembed_generation = f"deembed-{101 + index}"
        row[
            "deembedding_reference_plane_phase_causality_passivity_grid_generation_identity"
        ] = {
            "deembedding_generation": deembed_generation,
            "reference_plane_deembedding_generation": deembed_generation,
            "phase_deembedding_generation": deembed_generation,
            "causality_deembedding_generation": deembed_generation,
            "passivity_deembedding_generation": deembed_generation,
            "frequency_grid_deembedding_generation": deembed_generation,
            "result_deembedding_generation": deembed_generation,
            "port_mode_ids": ["P1:M1", "P2:M1"],
            "result_port_mode_ids": ["P1:M1", "P2:M1"],
            "reference_plane_offsets_m": [0.001, 0.001],
            "result_reference_plane_offsets_m": [0.001, 0.001],
            "frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "unwrapped_phase_rad": [-0.2, -0.35, -0.5],
            "result_unwrapped_phase_rad": [-0.2, -0.35, -0.5],
            "causality_check_passed": True,
            "result_causality_check_passed": True,
            "passivity_max_singular_values": [0.82, 0.84, 0.86],
            "result_passivity_max_singular_values": [0.82, 0.84, 0.86],
            "deembedded_network_sha256": "1" * 64,
            "result_deembedded_network_sha256": "1" * 64,
        }
        cosim_generation = f"field-circuit-{101 + index}"
        row[
            "field_circuit_cosim_port_sign_impedance_power_balance_generation_identity"
        ] = {
            "cosim_generation": cosim_generation,
            "field_port_cosim_generation": cosim_generation,
            "circuit_port_cosim_generation": cosim_generation,
            "sign_cosim_generation": cosim_generation,
            "impedance_cosim_generation": cosim_generation,
            "power_balance_cosim_generation": cosim_generation,
            "result_cosim_generation": cosim_generation,
            "port_id": "P1:M1",
            "result_port_id": "P1:M1",
            "current_sign_convention": "positive_into_field_port",
            "result_current_sign_convention": "positive_into_field_port",
            "voltage_reference": "positive_to_negative_terminal",
            "result_voltage_reference": "positive_to_negative_terminal",
            "port_voltage_ri_v": [10.0, 2.0],
            "result_port_voltage_ri_v": [10.0, 2.0],
            "port_current_ri_a": [0.2, -0.04],
            "result_port_current_ri_a": [0.2, -0.04],
            "port_impedance_ri_ohm": [46.15384615384615, 19.23076923076923],
            "result_port_impedance_ri_ohm": [46.15384615384615, 19.23076923076923],
            "phasor_amplitude_convention": "rms",
            "result_phasor_amplitude_convention": "rms",
            "field_absorbed_power_w": 1.92,
            "circuit_delivered_power_w": 1.92,
            "power_balance_residual_w": 0.0,
            "result_power_balance_residual_w": 0.0,
            "cosim_result_sha256": "2" * 64,
            "reported_cosim_result_sha256": "2" * 64,
        }
    return summary


def test_v24_public_positive_deembedding_and_field_circuit_cosim() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v24())["status"] == "ok"


def test_v24_public_deembedding_phase_causality_passivity_grid_mismatch() -> None:
    summary = _summary_v24()
    summary["runs"][0][
        "deembedding_reference_plane_phase_causality_passivity_grid_generation_identity"
    ].update(
        {
            "reference_plane_deembedding_generation": "deembed-100",
            "phase_deembedding_generation": "deembed-99",
            "causality_deembedding_generation": "deembed-98",
            "passivity_deembedding_generation": "deembed-97",
            "frequency_grid_deembedding_generation": "deembed-96",
            "result_reference_plane_offsets_m": [0.002, -0.001],
            "result_frequency_grid_hz": [1.0e9, 1.6e9, 2.0e9],
            "result_unwrapped_phase_rad": [-0.2, 5.9, -0.5],
            "result_causality_check_passed": False,
            "result_passivity_max_singular_values": [0.82, 1.12, 0.86],
            "result_deembedded_network_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "deembedded_network_uses_current_planes_phase_causality_passivity_and_grid"
    ]


def test_v24_public_field_circuit_cosim_sign_impedance_power_mismatch() -> None:
    summary = _summary_v24()
    summary["runs"][0][
        "field_circuit_cosim_port_sign_impedance_power_balance_generation_identity"
    ].update(
        {
            "field_port_cosim_generation": "field-circuit-100",
            "circuit_port_cosim_generation": "field-circuit-99",
            "sign_cosim_generation": "field-circuit-98",
            "impedance_cosim_generation": "field-circuit-97",
            "power_balance_cosim_generation": "field-circuit-96",
            "result_current_sign_convention": "positive_out_of_field_port",
            "result_voltage_reference": "negative_to_positive_terminal",
            "result_port_impedance_ri_ohm": [-46.15384615384615, -19.23076923076923],
            "circuit_delivered_power_w": -1.2,
            "result_power_balance_residual_w": 3.16,
            "reported_cosim_result_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "field_circuit_cosim_uses_current_sign_impedance_and_power_balance"
    ]
