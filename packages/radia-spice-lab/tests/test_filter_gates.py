from ltspice_converter.filter_gates import sallen_key_filter_family_gate


def _rows():
    return [
        {"id": "low_q", "R1_ohm": 10e3, "R2_ohm": 10e3,
         "C1_F": 10e-9, "C2_F": 10e-9, "dc_gain": 1.0,
         "peak_gain": 1.0, "minus3dB_frequency_Hz": 1024.718,
         "gain_at_100k": 1e-4, "ideal_complex_relative_l2_to_20k": 1e-4},
        {"id": "high_q", "R1_ohm": 10e3, "R2_ohm": 10e3,
         "C1_F": 40e-9, "C2_F": 10e-9, "dc_gain": 1.0,
         "peak_gain": 1.15, "minus3dB_frequency_Hz": 1019.0,
         "gain_at_100k": 1e-4, "ideal_complex_relative_l2_to_20k": 1e-4},
    ]


def _set_exact_cutoffs(rows):
    computed = sallen_key_filter_family_gate(rows, max_minus3db_relative_error=1.0)["rows"]
    for source, item in zip(rows, computed):
        source["minus3dB_frequency_Hz"] = item["ideal_minus3dB_frequency_Hz"]


def test_sallen_key_family_accepts_cutoff_identity_and_q_peaking():
    rows = _rows()
    _set_exact_cutoffs(rows)
    result = sallen_key_filter_family_gate(rows)
    assert result["status"] == "ok"
    assert result["checks"]["peaking_behavior_matches_quality_factor"] is True


def test_sallen_key_family_rejects_scalar_cutoff_match_with_bad_complex_response():
    rows = _rows()
    _set_exact_cutoffs(rows)
    rows[1]["ideal_complex_relative_l2_to_20k"] = 0.2
    rows[1]["peak_gain"] = 0.9
    result = sallen_key_filter_family_gate(rows)
    assert result["status"] == "needs_attention"
    assert result["checks"]["complex_response_matches_ideal"] is False
    assert result["checks"]["peaking_behavior_matches_quality_factor"] is False
