import copy

from ltspice_converter.mcp_server import monte_carlo_tolerance_family_gate as mcp_gate
from ltspice_converter.monte_carlo_gate import monte_carlo_tolerance_family_gate


def good_summary():
    return {
        "distribution": "independent_uniform_symmetric",
        "tolerance_fraction": 0.05,
        "sample_count": 1000,
        "independent_parts_per_equivalent": 5,
        "units": {"output": "V", "tolerance": "fraction"},
        "resistor_equivalent": {
            "single": {"nominal": 1.0, "mean": 1.000108, "standard_deviation": 0.029297},
            "series": {"nominal": 10.0, "mean": 9.999780, "standard_deviation": 0.131306},
            "parallel": {"nominal": 10.0, "mean": 9.993008, "standard_deviation": 0.131303},
        },
        "divider": {
            "single_per_arm": {"nominal": 5.0, "mean": 4.999696, "standard_deviation": 0.100901},
            "multi_per_arm": {"nominal": 5.0, "mean": 5.000352, "standard_deviation": 0.045686},
        },
        "gate_tolerances": {
            "minimum_samples": 500,
            "mean_relative_error": 0.01,
            "relative_sigma_theory_error": 0.12,
            "series_parallel_relative_sigma_difference": 0.05,
            "root_n_reduction_relative_error": 0.12,
        },
    }


def test_accepts_uniform_tolerance_root_n_family():
    result = monte_carlo_tolerance_family_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["relative_sigmas_match_uniform_theory"] is True
    assert mcp_gate(good_summary())["status"] == "ok"


def test_rejects_missing_multi_part_variance_reduction():
    bad = copy.deepcopy(good_summary())
    bad["divider"]["multi_per_arm"]["standard_deviation"] = 0.100901
    result = monte_carlo_tolerance_family_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["divider_family_follows_root_n_reduction"] is False
    assert result["checks"]["relative_sigmas_match_uniform_theory"] is False


def test_rejects_distribution_without_independence_contract():
    bad = good_summary()
    bad["distribution"] = "unspecified"
    result = monte_carlo_tolerance_family_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["uniform_independent_tolerance_declared"] is False
