import copy
import json

from radia_mcp.radia_ngsolve.server import moving_conductor_eddy_brake_gate


def good() -> dict:
    return {
        "units": {
            "time": "s",
            "displacement": "m",
            "velocity": "m/s",
            "force": "N",
            "loss": "W",
        },
        "observable_contract": {
            "force_family": "lorentz_force_on_moving_conductor",
            "loss_family": "joule_loss_in_moving_conductor",
            "force_value_kind": "absolute_magnitude",
        },
        "energy_balance_contract": {
            "magnetic_energy_rate_available": False,
            "mechanical_work_vs_joule_heat": "diagnostic_only",
        },
        "rows": [
            {
                "time_s": 0.1 * index,
                "displacement_m": 0.2 * index,
                "velocity_m_s": 2.0,
                "lorentz_force_n": 3.0 + index,
                "lorentz_force_parts_n": [1.0, 2.0 + index, 0.0],
                "joule_loss_w": 2.0 + 0.5 * index,
                "joule_loss_parts_w": [1.0, 1.0 + 0.5 * index, 0.0],
            }
            for index in range(7)
        ],
    }


def test_accepts_motion_and_exact_force_loss_decompositions():
    result = json.loads(moving_conductor_eddy_brake_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_force_decomposition_relative_error"] == 0.0
    assert result["metrics"]["maximum_loss_decomposition_relative_error"] == 0.0


def test_rejects_force_namespace_mix_and_bad_velocity():
    payload = copy.deepcopy(good())
    payload["rows"][3]["velocity_m_s"] = 2000.0
    payload["rows"][4]["lorentz_force_parts_n"][0] = 99.0
    result = json.loads(moving_conductor_eddy_brake_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["displacement_derivative_matches_velocity"] is False
    assert result["checks"]["lorentz_force_decomposition_closes"] is False


def test_rejects_unqualified_mechanical_equals_joule_claim():
    payload = good()
    payload["energy_balance_contract"] = {
        "magnetic_energy_rate_available": False,
        "mechanical_work_vs_joule_heat": "must_equal",
    }
    result = json.loads(moving_conductor_eddy_brake_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["missing_magnetic_energy_term_acknowledged"] is False
