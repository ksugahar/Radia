from copy import deepcopy

from radia_mcp.radia_ngsolve.network_artifact_identity_v51 import GROUP_DELAY, S_PARAMETER, validate_public_v51_identity


PROMOTED_CASE_IDS = {
    "v51_public_sparameter_renormalization_complex_zref_modal_impedance_wavebasis_port_owner_mismatch",
    "v51_public_group_delay_unwrap_frequency_derivative_smoothing_window_trace_owner_mismatch",
}


def _payload() -> dict[str, object]:
    generation = "network-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    frequencies = [1.0e9, 1.5e9, 2.0e9, 2.5e9, 3.0e9]
    phase = [0.0, -18.0, -36.0, -54.0, -72.0]
    delay = [1.0e-10] * 5
    run = {
        S_PARAMETER: {
            "generation": generation, "renormalization_generation": generation, "zref_generation": generation,
            "modal_generation": generation, "wave_generation": generation, "owner_generation": generation,
            "result_generation": generation, "renormalized": True, "result_renormalized": True,
            "complex_reference_impedance_ohm": [50.0, 5.0], "result_complex_reference_impedance_ohm": [50.0, 5.0],
            "modal_impedances_ohm": {"port1_mode1": [48.0, 2.0], "port2_mode1": [52.0, -1.0]},
            "result_modal_impedances_ohm": {"port1_mode1": [48.0, 2.0], "port2_mode1": [52.0, -1.0]},
            "wave_basis": "power_waves", "result_wave_basis": "power_waves", "port_owner": "port:network-v51",
            "result_port_owner": "port:network-v51", **result,
        },
        GROUP_DELAY: {
            "generation": generation, "unwrap_generation": generation, "frequency_generation": generation,
            "derivative_generation": generation, "smoothing_generation": generation, "owner_generation": generation,
            "result_generation": generation, "frequency_hz": frequencies, "result_frequency_hz": frequencies,
            "unwrapped_phase_deg": phase, "result_unwrapped_phase_deg": phase, "group_delay_s": delay,
            "result_group_delay_s": delay, "derivative_definition": "minus_dphi_rad_domega",
            "result_derivative_definition": "minus_dphi_rad_domega", "smoothing_window_points": 5,
            "result_smoothing_window_points": 5, "trace_owner": "trace:s21-v51", "result_trace_owner": "trace:s21-v51", **result,
        },
    }
    return {"runs": [deepcopy(run), deepcopy(run)]}


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_v51_identity(_payload()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    payload = _payload()
    payload["runs"][0][S_PARAMETER].update({"result_wave_basis": "pseudo_waves", "result_port_owner": "port:stale"})
    payload["runs"][0][GROUP_DELAY].update({"result_derivative_definition": "dphi_deg_df", "result_trace_owner": "trace:stale"})
    assert not all(validate_public_v51_identity(payload).values())


def test_v51_self_consistent_wrong_network_semantics_are_rejected() -> None:
    payload = _payload()
    for run in payload["runs"]:
        run[S_PARAMETER]["wave_basis"] = run[S_PARAMETER]["result_wave_basis"] = "pseudo_waves"
        run[GROUP_DELAY]["smoothing_window_points"] = run[GROUP_DELAY]["result_smoothing_window_points"] = 4
    assert not all(validate_public_v51_identity(payload).values())
