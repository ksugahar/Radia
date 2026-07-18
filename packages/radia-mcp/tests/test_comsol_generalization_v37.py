from __future__ import annotations

import math

from test_comsol_generalization_v36 import (
    _summary,
    _with_v36_force_and_modal_identity,
    gate,
)

_PROMOTED_CASE_IDS = (
    "v37_public_capacitance_matrix_charge_energy_gauge_reciprocity_terminal_owner_mismatch",
    "v37_public_thermoelastic_harmonic_phase_loss_work_temperature_displacement_mesh_mismatch",
)


def _with_v37_capacitance_and_thermoelastic_identity(summary: dict) -> dict:
    summary = _with_v36_force_and_modal_identity(summary)
    generation = "capacitance-closure-241"
    matrix = [[2.0e-12, -2.0e-12], [-2.0e-12, 2.0e-12]]
    summary[
        "capacitance_matrix_charge_energy_gauge_reciprocity_terminal_mesh_result_generation_identity"
    ] = {
        "capacitance_generation": generation,
        **{
            key: generation
            for key in (
                "matrix_generation",
                "charge_generation",
                "energy_generation",
                "gauge_generation",
                "reciprocity_generation",
                "terminal_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "terminal_names": ["terminal_1", "reference_2"],
        "result_terminal_names": ["terminal_1", "reference_2"],
        "reference_terminal": "reference_2",
        "result_reference_terminal": "reference_2",
        "capacitance_matrix_f": matrix,
        "result_capacitance_matrix_f": matrix,
        "terminal_potential_v": [1.0, 0.0],
        "result_terminal_potential_v": [1.0, 0.0],
        "terminal_charge_c": [2.0e-12, -2.0e-12],
        "result_terminal_charge_c": [2.0e-12, -2.0e-12],
        "stored_energy_j": 1.0e-12,
        "result_stored_energy_j": 1.0e-12,
        "reciprocity_tolerance": 1.0e-12,
        "result_reciprocity_tolerance": 1.0e-12,
        "terminal_owner": "comp1/es/terminals-241",
        "accepted_terminal_owner": "comp1/es/terminals-241",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "thermoelastic-harmonic-241"
    heat = [10.0, 2.0]
    summary[
        "thermoelastic_harmonic_heat_phase_temperature_displacement_work_loss_frequency_mesh_result_generation_identity"
    ] = {
        "thermoelastic_generation": generation,
        **{
            key: generation
            for key in (
                "heat_generation",
                "phase_generation",
                "temperature_generation",
                "displacement_generation",
                "work_generation",
                "loss_generation",
                "frequency_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "frequency_hz": 1000.0,
        "result_frequency_hz": 1000.0,
        "heat_source_complex_w": heat,
        "result_heat_source_complex_w": heat,
        "heat_source_phase_rad": math.atan2(heat[1], heat[0]),
        "result_heat_source_phase_rad": math.atan2(heat[1], heat[0]),
        "temperature_complex_k": [5.0, 1.0],
        "result_temperature_complex_k": [5.0, 1.0],
        "displacement_complex_m": [1.0e-6, -2.0e-7],
        "result_displacement_complex_m": [1.0e-6, -2.0e-7],
        "thermal_expansion_work_j": 2.0e-3,
        "result_thermal_expansion_work_j": 2.0e-3,
        "mechanical_loss_j": 1.0e-4,
        "result_mechanical_loss_j": 1.0e-4,
        "loss_convention": "positive_dissipated_per_cycle",
        "result_loss_convention": "positive_dissipated_per_cycle",
        "mesh_owner": "comp1/mesh1:thermoelastic-241",
        "accepted_mesh_owner": "comp1/mesh1:thermoelastic-241",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v37_public_positive_capacitance_and_thermoelastic_contracts() -> None:
    assert gate(_with_v37_capacitance_and_thermoelastic_identity(_summary()))["status"] == "ok"


def test_v37_public_capacitance_matrix_charge_energy_gauge_reciprocity_terminal_owner_mismatch() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "capacitance_matrix_charge_energy_gauge_reciprocity_terminal_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "matrix_generation": "capacitance-closure-240",
            "result_reference_terminal": "terminal_1",
            "result_capacitance_matrix_f": [[2.0e-12, 1.0e-12], [-2.0e-12, 1.0e-12]],
            "result_terminal_charge_c": [1.0e-12, 1.0e-12],
            "result_stored_energy_j": -1.0e-12,
            "accepted_terminal_owner": "comp1/es/old-terminals",
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "capacitance_results_use_current_matrix_charge_energy_gauge_reciprocity_terminals_mesh_and_result"
    ]


def test_v37_public_thermoelastic_harmonic_phase_loss_work_temperature_displacement_mesh_mismatch() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "thermoelastic_harmonic_heat_phase_temperature_displacement_work_loss_frequency_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "phase_generation": "thermoelastic-harmonic-240",
            "result_frequency_hz": -1000.0,
            "result_heat_source_phase_rad": -2.0,
            "result_temperature_complex_k": [-5.0, -1.0],
            "result_displacement_complex_m": [-1.0e-6, 2.0e-7],
            "result_thermal_expansion_work_j": -2.0e-3,
            "result_mechanical_loss_j": -1.0e-4,
            "accepted_mesh_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thermoelastic_harmonics_use_current_heat_phase_temperature_displacement_work_loss_frequency_mesh_and_result"
    ]


def test_v37_public_rejects_self_consistent_nonzero_capacitance_row_sum() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "capacitance_matrix_charge_energy_gauge_reciprocity_terminal_mesh_result_generation_identity"
    ]
    wrong = [[2.0e-12, -1.0e-12], [-1.0e-12, 2.0e-12]]
    identity["capacitance_matrix_f"] = wrong
    identity["result_capacitance_matrix_f"] = wrong
    assert gate(summary)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_negative_thermoelastic_loss() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "thermoelastic_harmonic_heat_phase_temperature_displacement_work_loss_frequency_mesh_result_generation_identity"
    ]
    identity["mechanical_loss_j"] = -1.0e-4
    identity["result_mechanical_loss_j"] = -1.0e-4
    assert gate(summary)["status"] == "needs_attention"
