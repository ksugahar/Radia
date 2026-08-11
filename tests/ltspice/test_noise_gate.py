from __future__ import annotations

from copy import deepcopy

from radia.ltspice.mcp_server import rc_thermal_noise_psd_gate as mcp_gate
from radia.ltspice.noise_gate import rc_thermal_noise_psd_gate


def good_summary() -> dict:
    first = {
        "capacitance_f": 1.0e-6,
        "point_count": 601,
        "low_frequency_density": 9.1038688e-10,
        "analytic_low_frequency_density": 9.1038646e-10,
        "numeric_psd_integrated_rms": 6.4305206e-8,
        "analytic_finite_band_rms": 6.4302346e-8,
        "noise_measure_integrated": 6.4304144e-8,
        "direct_density_integral": 1.8671358e-5,
    }
    second = {
        "capacitance_f": 0.6366198e-6,
        "point_count": 601,
        "low_frequency_density": 9.1038693e-10,
        "analytic_low_frequency_density": 9.1038649e-10,
        "numeric_psd_integrated_rms": 8.0550823e-8,
        "analytic_finite_band_rms": 8.0547248e-8,
        "noise_measure_integrated": 8.0549495e-8,
        "direct_density_integral": 2.7273679e-5,
    }
    return {
        "analysis": "noise",
        "rms_integration": "sqrt_integral_density_squared_df",
        "measure_integ_semantics": "noise_rms_not_ordinary_integral",
        "resistance_ohm": 50.0,
        "temperature_k": 300.15,
        "frequency_start_hz": 1.0,
        "frequency_stop_hz": 1.0e6,
        "units": {
            "spectral_density": "V/sqrt(Hz)",
            "rms_noise": "V",
            "ordinary_density_integral": "V*sqrt(Hz)",
        },
        "cases": [first, second],
    }


def test_noise_gate_accepts_density_rms_and_capacitance_pair():
    result = rc_thermal_noise_psd_gate(good_summary())
    assert result["status"] == "ok"
    assert mcp_gate(good_summary())["status"] == "ok"


def test_noise_gate_rejects_ordinary_integral_semantics():
    bad = deepcopy(good_summary())
    bad["measure_integ_semantics"] = "ordinary_integral"
    result = rc_thermal_noise_psd_gate(bad)
    assert result["status"] == "needs_attention"
    assert "noise_measure_context_semantics_explicit" in result["issues"]


def test_noise_gate_rejects_density_used_as_rms():
    bad = deepcopy(good_summary())
    bad["cases"][0]["noise_measure_integrated"] = bad["cases"][0][
        "direct_density_integral"
    ]
    result = rc_thermal_noise_psd_gate(bad)
    assert result["status"] == "needs_attention"
    assert "noise_measure_matches_psd_rms" in result["issues"]


def test_noise_gate_rejects_dimensionally_wrong_units():
    bad = deepcopy(good_summary())
    bad["units"]["spectral_density"] = "V"
    result = rc_thermal_noise_psd_gate(bad)
    assert result["status"] == "needs_attention"
    assert "spectral_and_integrated_units_distinct" in result["issues"]
