import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def _load(name):
    return json.loads((HERE / f"{name}.json").read_text(encoding="utf-8"))


def test_excitation_invariant_field_evidence():
    data = _load("excitation_invariant_field")
    derived = {
        "optimized_direction_drift_below_flat": (
            data["optimized"]["D_dir"] < data["flat_cut"]["D_dir"]
        ),
        "invariance_improvement_gt_5": data["invariance_factor"] > 5.0,
        "linear_control_direction_drift_lt_1e_12": (
            max(data["flat_invariance_curve"]["D_dir_linear_control"])
            < 1e-12
        ),
    }

    assert data["schema"] == (
        "radia.validation.clebsch_hodograph.excitation_invariant_field.v1"
    )
    assert data["checks"] == derived
    assert all(derived.values())


def test_hodograph_bending_evidence():
    data = _load("hodograph_bending_sy")
    derived = {
        "fem_field_relative_error_lt_2e_4": data["fem_verify"]["By_rel_err"] < 2e-4,
        "pole_match_error_lt_2e_6": data["fem_verify"]["pole_match_err"] < 2e-6,
        "optimized_peak_lt_1p02": data["design_curve"]["best_peak"] < 1.02,
    }

    assert data["schema"] == (
        "radia.validation.clebsch_hodograph.hodograph_bending_sy.v1"
    )
    assert data["checks"] == derived
    assert all(derived.values())


def test_hodograph_feasibility_evidence():
    data = _load("hodograph_feasibility_2d")
    derived = {
        "fem_field_relative_error_lt_2e_4": data["fem_verify"]["By_rel_err"] < 2e-4,
        "pole_match_error_lt_2e_6": data["fem_verify"]["pole_match_err"] < 2e-6,
        "vonmises_pde_residual_lt_1e_3": data["vonmises_chart"]["pde_resid"] < 1e-3,
        "vonmises_chart_single_valued": data["vonmises_chart"]["single_valued"],
    }

    assert data["schema"] == (
        "radia.validation.clebsch_hodograph.hodograph_feasibility_2d.v1"
    )
    assert data["checks"] == derived
    assert all(derived.values())


def test_edge_focusing_tracking_evidence():
    data = _load("edge_focusing_tracking")
    summary = data["summary"]
    derived = {
        "finest_width_slope_within_0p02": abs(summary["finest_w_slope"] - 1.0) < 0.02,
        "max_relative_error_vs_enge_lt_0p01": summary["max_rel_err_vs_enge"] < 0.01,
        "beta0_baseline_abs_lt_0p02": abs(summary["beta0_baseline"]) < 0.02,
    }

    assert data["schema"] == (
        "radia.validation.clebsch_hodograph.edge_focusing_tracking.v1"
    )
    assert data["checks"] == derived
    assert all(derived.values())


def test_edge_focusing_fem_evidence():
    data = _load("edge_focusing_fem_results")
    headline = data["headline"]
    cross_check = data["hdiv_vim_cross_check"]
    scaled = [row["rho_m"] * row["dK_in_fem"] for row in data["rho_sweep"]]
    spread = (max(scaled) - min(scaled)) / (sum(scaled) / len(scaled))
    derived = {
        "rho_scaled_focusing_spread_lt_0p02": spread < 0.02,
        "matched_engine_relative_difference_lt_0p01": (
            cross_check["engine_agreement"]["rel_diff_hdiv2"] < 0.01
        ),
        "fem_to_model_ratio_between_0p9_and_0p97": (
            0.9 < headline["fem_over_model"] < 0.97
        ),
        "spurious_beta0_term_abs_lt_0p002": (
            abs(headline["spurious_T1_in_beta0"]) < 0.002
        ),
    }

    assert data["schema"] == (
        "radia.validation.clebsch_hodograph.edge_focusing_fem.v1"
    )
    assert data["checks"] == derived
    assert all(derived.values())
