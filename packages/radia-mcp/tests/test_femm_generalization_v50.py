from copy import deepcopy

from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v50 import CIRCUIT, IRON, validate_public_identity


PROMOTED_CASE_IDS = {
    "v50_public_complex_current_phasor_circuit_series_parallel_turns_depth_owner_mismatch",
    "v50_public_iron_loss_steinmetz_frequency_flux_waveform_component_owner_mismatch",
}


def _identity() -> dict[str, object]:
    generation = "femm-v50"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    waveform = [0.0, 0.7, 1.2, 0.7, 0.0, -0.7, -1.2, -0.7, 0.0]
    components = {"hysteresis_w": 1.8, "eddy_w": 0.9, "excess_w": 0.3, "total_w": 3.0}
    return {
        CIRCUIT: {
            "generation": generation, "current_generation": generation, "connection_generation": generation,
            "turn_generation": generation, "depth_generation": generation, "owner_generation": generation,
            "result_generation": generation, "current_phasor_a": [3.0, -1.0], "result_current_phasor_a": [3.0, -1.0],
            "winding_connection": "series", "result_winding_connection": "series", "parallel_paths": 1,
            "result_parallel_paths": 1, "turns": 120, "result_turns": 120, "depth_m": 0.04, "result_depth_m": 0.04,
            "circuit_owner": "circuit:coil-v50", "result_circuit_owner": "circuit:coil-v50", **result,
        },
        IRON: {
            "generation": generation, "coefficient_generation": generation, "frequency_generation": generation,
            "waveform_generation": generation, "component_generation": generation, "owner_generation": generation,
            "result_generation": generation, "steinmetz_coefficients": {"kh": 0.021, "ke": 0.00018, "alpha": 1.62},
            "result_steinmetz_coefficients": {"kh": 0.021, "ke": 0.00018, "alpha": 1.62},
            "frequency_hz": 400.0, "result_frequency_hz": 400.0, "flux_density_waveform_t": waveform,
            "result_flux_density_waveform_t": waveform, "loss_components": components, "result_loss_components": components,
            "material_owner": "material:lamination-v50", "result_material_owner": "material:lamination-v50", **result,
        },
    }


def test_v50_positive_harmonic_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v50_circuit_and_iron_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[CIRCUIT]["result_parallel_paths"] = 2
    identity[IRON]["result_frequency_hz"] = 50.0
    assert not all(validate_public_identity(identity).values())


def test_v50_self_consistent_invalid_waveform_and_loss_closure_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[IRON]["flux_density_waveform_t"] = identity[IRON]["result_flux_density_waveform_t"] = [0.0, 1.0, 0.0]
    identity[IRON]["loss_components"] = identity[IRON]["result_loss_components"] = {
        "hysteresis_w": 1.8, "eddy_w": 0.9, "excess_w": 0.3, "total_w": 8.0
    }
    assert validate_public_identity(identity)["v50_steinmetz_frequency_waveform_components_material_owner"] is False
