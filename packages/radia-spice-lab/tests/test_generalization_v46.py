from ltspice_converter.ltspice_v46_gates import validate_ltspice_v46_identity


PROMOTED_CASE_IDS = (
    "v46_public_buck_failed_timestep_raw_truncation_measure_nan_inf_mismatch",
    "v46_public_noise_monte_carlo_seed_sample_filter_psd_bin_partial_mismatch",
)


def _positive():
    return {
        "buck_v46_failed_timestep_raw_truncation_measure_nan_inf_identity": {
            "generation_id": "buck-v46",
            "result_generation_id": "buck-v46",
            "failed_timestep_index": 7,
            "result_failed_timestep_index": 7,
            "raw_truncated": False,
            "result_raw_truncated": False,
            "raw_point_count": 5000,
            "result_raw_point_count": 5000,
            "expected_raw_point_count": 5000,
            "result_expected_raw_point_count": 5000,
            "measure_status": "complete",
            "result_measure_status": "complete",
            "nonfinite_count": 0,
            "result_nonfinite_count": 0,
            "analysis_status": "completed",
            "result_analysis_status": "completed",
            "waveform_owner": "waveform/current",
            "result_waveform_owner": "waveform/current",
            "raw_sha256": "a" * 64,
            "result_raw_sha256": "a" * 64,
            "result_sha256": "b" * 64,
            "accepted_result_sha256": "b" * 64,
        },
        "noise_v46_monte_carlo_seed_sample_filter_psd_bin_partial_identity": {
            "generation_id": "noise-v46",
            "result_generation_id": "noise-v46",
            "random_seed": 20260719,
            "result_random_seed": 20260719,
            "sample_count": 64,
            "result_sample_count": 64,
            "sample_filter": "discard_first_8",
            "result_sample_filter": "discard_first_8",
            "psd_bin_hz": 1000.0,
            "result_psd_bin_hz": 1000.0,
            "partial_sweep": False,
            "result_partial_sweep": False,
            "sweep_status": "complete",
            "result_sweep_status": "complete",
            "nonfinite_count": 0,
            "result_nonfinite_count": 0,
            "measure_owner": "measure/current",
            "result_measure_owner": "measure/current",
            "raw_sha256": "c" * 64,
            "result_raw_sha256": "c" * 64,
            "result_sha256": "d" * 64,
            "accepted_result_sha256": "d" * 64,
        },
    }


def test_v46_public_replay_identity():
    assert validate_ltspice_v46_identity(_positive()) is True


def test_v46_public_rejects_nonfinite_count_mutation():
    value = _positive()
    value["buck_v46_failed_timestep_raw_truncation_measure_nan_inf_identity"]["result_nonfinite_count"] = 1
    assert validate_ltspice_v46_identity(value) is False
