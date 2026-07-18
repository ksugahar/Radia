from __future__ import annotations

import math

from test_comsol_generalization_v35 import (
    _summary,
    _with_v35_multirate_and_adjoint_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v36_public_magnetostatic_virtual_work_coenergy_force_displacement_mesh_owner_mismatch",
    "v36_public_acoustic_modal_participation_effective_mass_damping_reconstruction_mismatch",
)


def _modal_response(
    response_frequencies: list[float],
    mode_frequencies: list[float],
    damping: list[float],
    participation: list[float],
    probe_factors: list[float],
) -> list[list[float]]:
    values = []
    for frequency in response_frequencies:
        omega = 2.0 * math.pi * frequency
        response = 0.0j
        for mode_hz, zeta, factor, probe in zip(
            mode_frequencies, damping, participation, probe_factors, strict=True
        ):
            omega_mode = 2.0 * math.pi * mode_hz
            response += probe * factor / complex(
                omega_mode**2 - omega**2,
                2.0 * zeta * omega_mode * omega,
            )
        values.append([response.real, response.imag])
    return values


def _with_v36_force_and_modal_identity(summary: dict) -> dict:
    summary = _with_v35_multirate_and_adjoint_identity(summary)
    generation = "virtual-work-coenergy-231"
    displacement = [-1.0e-4, 0.0, 1.0e-4]
    coenergy = [0.4997, 0.5, 0.5003]
    force = (coenergy[2] - coenergy[0]) / (displacement[2] - displacement[0])
    summary[
        "magnetostatic_virtual_work_coenergy_force_displacement_current_mesh_frame_solution_result_generation_identity"
    ] = {
        "force_generation": generation,
        **{
            key: generation
            for key in (
                "displacement_generation",
                "coenergy_generation",
                "current_generation",
                "mesh_generation",
                "frame_generation",
                "solution_generation",
                "result_generation",
            )
        },
        "displacement_m": displacement,
        "result_displacement_m": displacement,
        "coenergy_j": coenergy,
        "result_coenergy_j": coenergy,
        "held_source_convention": "constant_current",
        "result_held_source_convention": "constant_current",
        "force_sign_convention": "positive_dcoenergy_dx",
        "result_force_sign_convention": "positive_dcoenergy_dx",
        "central_coenergy_force_n": force,
        "result_force_n": force,
        "force_tolerance_n": 1.0e-9,
        "result_force_tolerance_n": 1.0e-9,
        "coordinate_frame": "stationary_cartesian_x",
        "result_coordinate_frame": "stationary_cartesian_x",
        "displaced_mesh_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "result_displaced_mesh_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "force_solution_owner": "std1/sol1:parametric_displacement",
        "result_force_solution_owner": "std1/sol1:parametric_displacement",
        "force_result_sha256": "4" * 64,
        "accepted_force_result_sha256": "4" * 64,
    }

    generation = "acoustic-modal-participation-231"
    frequencies = [100.0, 160.0]
    masses = [2.0, 1.5]
    participation = [0.5, 0.4]
    damping = [0.01, 0.02]
    probes = [1.0, 0.8]
    response_frequencies = [90.0, 120.0, 180.0]
    response = _modal_response(
        response_frequencies, frequencies, damping, participation, probes
    )
    effective_masses = [
        factor**2 * mass
        for factor, mass in zip(participation, masses, strict=True)
    ]
    summary[
        "acoustic_modal_normalization_effective_mass_participation_damping_frequency_reconstruction_mesh_result_generation_identity"
    ] = {
        "modal_generation": generation,
        **{
            key: generation
            for key in (
                "normalization_generation",
                "mass_generation",
                "participation_generation",
                "damping_generation",
                "frequency_generation",
                "reconstruction_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "normalization": "peak_displacement",
        "result_normalization": "peak_displacement",
        "mode_frequency_hz": frequencies,
        "result_mode_frequency_hz": frequencies,
        "modal_mass_kg": masses,
        "result_modal_mass_kg": masses,
        "participation_factor": participation,
        "result_participation_factor": participation,
        "effective_modal_mass_kg": effective_masses,
        "result_effective_modal_mass_kg": effective_masses,
        "damping_ratio": damping,
        "result_damping_ratio": damping,
        "probe_mode_factor": probes,
        "result_probe_mode_factor": probes,
        "response_frequency_hz": response_frequencies,
        "result_response_frequency_hz": response_frequencies,
        "probe_response_complex": response,
        "result_probe_response_complex": response,
        "response_tolerance": 1.0e-12,
        "result_response_tolerance": 1.0e-12,
        "modal_mesh_sha256": "5" * 64,
        "result_modal_mesh_sha256": "5" * 64,
        "modal_result_sha256": "6" * 64,
        "accepted_modal_result_sha256": "6" * 64,
    }
    return summary


def test_v36_public_positive_force_and_modal_contracts() -> None:
    assert gate(_with_v36_force_and_modal_identity(_summary()))["status"] == "ok"


def test_v36_public_magnetostatic_virtual_work_coenergy_force_displacement_mesh_owner_mismatch() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    summary[
        "magnetostatic_virtual_work_coenergy_force_displacement_current_mesh_frame_solution_result_generation_identity"
    ].update(
        {
            "coenergy_generation": "virtual-work-coenergy-230",
            "mesh_generation": "virtual-work-coenergy-229",
            "result_generation": "virtual-work-coenergy-228",
            "result_displacement_m": [-2.0e-4, 0.0, 1.0e-4],
            "result_coenergy_j": [0.4999, 0.5, 0.5001],
            "result_held_source_convention": "constant_voltage",
            "result_force_sign_convention": "negative_dcoenergy_dx",
            "result_force_n": -3.0,
            "result_coordinate_frame": "rotor_cylindrical",
            "result_displaced_mesh_sha256": ["a" * 64, "b" * 64, "c" * 64],
            "result_force_solution_owner": "std_old/sol0",
            "accepted_force_result_sha256": "d" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetostatic_force_uses_current_displacement_coenergy_source_sign_mesh_frame_owner_and_result"
    ]


def test_v36_public_acoustic_modal_participation_effective_mass_damping_reconstruction_mismatch() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    summary[
        "acoustic_modal_normalization_effective_mass_participation_damping_frequency_reconstruction_mesh_result_generation_identity"
    ].update(
        {
            "normalization_generation": "acoustic-modal-participation-230",
            "result_generation": "acoustic-modal-participation-228",
            "result_normalization": "mass_normalized",
            "result_mode_frequency_hz": [160.0, 100.0],
            "result_modal_mass_kg": [1.0, 1.0],
            "result_participation_factor": [0.4, -0.5],
            "result_effective_modal_mass_kg": [0.16, 0.25],
            "result_damping_ratio": [-0.02, 0.01],
            "result_response_frequency_hz": [90.0, 125.0, 180.0],
            "result_probe_response_complex": [[0.0, 0.0]] * 3,
            "result_modal_mesh_sha256": "e" * 64,
            "accepted_modal_result_sha256": "f" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "acoustic_modes_use_current_normalization_mass_participation_damping_reconstruction_mesh_and_result"
    ]


def test_v36_public_rejects_self_consistent_wrong_coenergy_derivative() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    identity = summary[
        "magnetostatic_virtual_work_coenergy_force_displacement_current_mesh_frame_solution_result_generation_identity"
    ]
    identity["coenergy_j"] = [0.4999, 0.5, 0.5001]
    identity["result_coenergy_j"] = [0.4999, 0.5, 0.5001]
    assert gate(summary)["status"] == "needs_attention"


def test_v36_public_rejects_self_consistent_wrong_effective_mass() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    identity = summary[
        "acoustic_modal_normalization_effective_mass_participation_damping_frequency_reconstruction_mesh_result_generation_identity"
    ]
    identity["effective_modal_mass_kg"] = [0.6, 0.24]
    identity["result_effective_modal_mass_kg"] = [0.6, 0.24]
    assert gate(summary)["status"] == "needs_attention"
