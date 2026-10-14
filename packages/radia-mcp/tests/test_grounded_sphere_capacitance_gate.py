import copy
import json

import pytest

from radia_mcp.radia_ngsolve.grounded_sphere_capacitance_gate import (
    grounded_sphere_capacitance_convergence_gate,
)
from radia_mcp.radia_ngsolve.server import (
    grounded_sphere_capacitance_convergence_gate as mcp_gate,
)


ANALYTIC = 4.475021350377041e-9


def _row(name, nodes, elements, voltage, charge_capacitance, volume_capacitance, corrected_capacitance):
    voltage_squared_half = 0.5 * voltage * voltage
    return {
        "case": name,
        "conductor": [voltage, charge_capacitance * voltage],
        "stored_energy_J": volume_capacitance * voltage_squared_half,
        "mixed_boundary_energy_J": (
            corrected_capacitance - volume_capacitance
        )
        * voltage_squared_half,
        "node_count": nodes,
        "element_count": elements,
    }


def _summary():
    fine = (4.4768183860739025e-9, 4.4065354321089306e-9, 4.476841341968186e-9)
    return {
        "problem_contract": {
            "analysis": "electrostatics",
            "problem_type": "axisymmetric",
            "length_units": "meters",
            "sphere_radius_m": 25.0,
            "sphere_center_height_m": 35.0,
            "ground_plane_voltage_V": 0.0,
            "outer_radius_m": 150.0,
            "open_boundary_c0_F_per_m2": 1.1805583756826701e-13,
            "open_boundary_asymptotic_order": 2,
        },
        "analytic": {"capacitance_F": ANALYTIC},
        "cases": [
            _row("coarse", 3920, 7433, 100.0, 4.483762721130678e-9, 4.413308881518958e-9, 4.4837811889149265e-9),
            _row("medium", 5248, 10085, 100.0, 4.481293406514265e-9, 4.41090042552528e-9, 4.481311830307388e-9),
            _row("fine", 18633, 36668, 100.0, *fine),
            _row("fine_repeat", 18633, 36668, 100.0, *fine),
            _row("fine_negative", 18633, 36668, -100.0, *fine),
        ],
    }


def test_accepts_image_series_refinement_boundary_energy_replay_and_sign_covariance():
    result = grounded_sphere_capacitance_convergence_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["analytic_relative_errors"][-1] < 5.0e-4
    assert result["metrics"]["minimum_uncorrected_energy_relative_error"] > 1.0e-2
    assert result["metrics"]["maximum_corrected_energy_relative_error"] < 1.0e-5


def test_rejects_missing_mixed_boundary_energy():
    bad = _summary()
    bad["cases"][2]["mixed_boundary_energy_J"] = 0.0
    result = grounded_sphere_capacitance_convergence_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["mixed_boundary_energy_restores_charge_energy_identity"] is False
    assert result["checks"]["electrostatic_capacitances_and_energies_are_positive"] is False


def test_rejects_stale_voltage_sign_reversal():
    bad = copy.deepcopy(_summary())
    bad["cases"][4]["conductor"][1] *= -0.8
    result = grounded_sphere_capacitance_convergence_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["voltage_reversal_has_odd_charge_and_even_energy"] is False


def test_rejects_positive_charge_at_negative_voltage():
    bad = copy.deepcopy(_summary())
    bad["cases"][4]["conductor"][1] = abs(bad["cases"][4]["conductor"][1])
    result = grounded_sphere_capacitance_convergence_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["voltage_reversal_has_odd_charge_and_even_energy"] is False


def test_rejects_a_sphere_intersecting_the_ground_plane():
    bad = _summary()
    bad["problem_contract"]["sphere_center_height_m"] = 20.0
    with pytest.raises(ValueError, match="must exceed"):
        grounded_sphere_capacitance_convergence_gate(bad)


def test_mcp_wrapper_returns_json_and_handles_invalid_input():
    good = json.loads(mcp_gate(json.dumps(_summary())))
    assert good["status"] == "ok"
    bad = json.loads(mcp_gate("[]"))
    assert bad["status"] == "invalid_input"
    bad_tolerance = json.loads(
        mcp_gate(json.dumps(_summary()), max_corrected_energy_relative_error=-1.0)
    )
    assert bad_tolerance["status"] == "invalid_input"


@pytest.mark.parametrize(
    "case_id",
    ["problem_type", "open_order", "case_order", "replay_charge", "mixed_energy"],
)
def test_counterfactual_curriculum90_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "problem_type":
        bad["problem_contract"]["problem_type"] = "planar"
    elif case_id == "open_order":
        bad["problem_contract"]["open_boundary_asymptotic_order"] = 1
    elif case_id == "case_order":
        bad["cases"][3]["case"] = "repeat_unknown"
    elif case_id == "replay_charge":
        bad["cases"][3]["conductor"][1] *= 1.1
    else:
        bad["cases"][2]["mixed_boundary_energy_J"] = 0.0
    result = json.loads(mcp_gate(json.dumps(bad)))
    assert result["status"] in {"needs_attention", "invalid_input"}


def test_generalization_v3s_rejects_open_boundary_not_enclosing_conductor():
    bad = copy.deepcopy(_summary())
    bad["problem_contract"]["outer_radius_m"] = 50.0
    result = json.loads(mcp_gate(json.dumps(bad)))
    assert result["status"] in {"needs_attention", "invalid_input"}


@pytest.mark.parametrize(
    "case_id",
    ["v4_ground_voltage", "v4_length_units", "v4_analytic_capacitance", "v4_refinement_elements", "v4_negative_stored_energy"],
)
def test_counterfactual_curriculum90_v4_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v4_ground_voltage":
        bad["problem_contract"]["ground_plane_voltage_V"] = 1.0
    elif case_id == "v4_length_units":
        bad["problem_contract"]["length_units"] = "millimeters"
    elif case_id == "v4_analytic_capacitance":
        bad["analytic"]["capacitance_F"] *= 1.2
    elif case_id == "v4_refinement_elements":
        bad["cases"][2]["element_count"] = bad["cases"][1]["element_count"]
    else:
        bad["cases"][2]["stored_energy_J"] = -1.0
    result = json.loads(mcp_gate(json.dumps(bad)))
    assert result["status"] in {"needs_attention", "invalid_input"}


def test_generalization_v5_rejects_stale_radius_dependent_analytic_value():
    bad = copy.deepcopy(_summary())
    bad["problem_contract"]["sphere_radius_m"] *= 1.5
    result = json.loads(mcp_gate(json.dumps(bad)))
    assert result["status"] in {"needs_attention", "invalid_input"}


@pytest.mark.parametrize(
    "case_id",
    [
        "v6_public_voltage_reversal_energy_drift",
        "v6_public_charge_energy_disagreement",
    ],
)
def test_generalization_v6_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v6_public_voltage_reversal_energy_drift":
        bad["cases"][4]["stored_energy_J"] *= 1.01
    else:
        bad["cases"][2]["conductor"][1] *= 1.02
    assert grounded_sphere_capacitance_convergence_gate(bad)["status"] == "needs_attention"
