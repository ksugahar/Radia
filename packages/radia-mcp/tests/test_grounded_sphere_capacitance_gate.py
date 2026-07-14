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
