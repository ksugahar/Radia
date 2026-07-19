from __future__ import annotations

from radia_mcp.radia_ngsolve.cst_v46_identity import validate_public_v46_identity


def _payload():
    generation = "test-846"
    return {"runs": [{
        "v46_public_time_domain_port_wave_impedance_unit_scale_partial_trace_mismatch": {
            "generation": generation, **{key: generation for key in ("time_generation", "trace_generation", "impedance_generation", "unit_scale_generation", "partial_generation", "result_generation")},
            "time_s": [0.0, 1.0e-10, 2.0e-10], "result_time_s": [0.0, 1.0e-10, 2.0e-10], "port_wave_trace_v": [0.0, 0.5, 1.0], "result_port_wave_trace_v": [0.0, 0.5, 1.0], "wave_impedance_unit": "ohm", "result_wave_impedance_unit": "ohm", "unit_scale_to_si": 1.0, "result_unit_scale_to_si": 1.0, "partial_trace": False, "result_partial_trace": False, "port_owner": "port:test", "result_port_owner": "port:test", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v46_public_field_monitor_coordinate_frame_sampling_window_nan_inf_mismatch": {
            "generation": generation, **{key: generation for key in ("frame_generation", "sampling_generation", "finite_generation", "monitor_generation", "result_generation")},
            "coordinate_frame": "global_cartesian", "result_coordinate_frame": "global_cartesian", "sampling_window_s": [0.0, 1.0e-9], "result_sampling_window_s": [0.0, 1.0e-9], "field_samples": [0.0, 0.2, 0.4], "result_field_samples": [0.0, 0.2, 0.4], "nonfinite_sample_count": 0, "result_nonfinite_sample_count": 0, "monitor_owner": "monitor:test", "result_monitor_owner": "monitor:test", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }]}


def test_v46_public_cst_identity_accepts_closed_artifacts():
    checks = validate_public_v46_identity(_payload())
    assert checks and all(checks.values())


def test_v46_public_cst_identity_rejects_unit_partial_and_nonfinite_mutations():
    payload = _payload()
    payload["runs"][0]["v46_public_time_domain_port_wave_impedance_unit_scale_partial_trace_mismatch"]["result_partial_trace"] = True
    payload["runs"][0]["v46_public_field_monitor_coordinate_frame_sampling_window_nan_inf_mismatch"]["result_nonfinite_sample_count"] = 1
    checks = validate_public_v46_identity(payload)
    assert checks and not all(checks.values())
