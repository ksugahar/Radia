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
