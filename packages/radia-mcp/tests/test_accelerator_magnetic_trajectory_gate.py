import copy
import json

from radia_mcp.accelerator.magnetic_trajectory_gate import evaluate_magnetic_trajectory_pair
from radia_mcp.accelerator.server import accelerator_magnetic_trajectory_pair_gate


def _trajectory(pid, x_final):
    return {
        "emission_id": 0,
        "particle_id": pid,
        "point_count": 4,
        "position_initial": [float(pid), 0.0, 0.0],
        "position_final": [x_final, 0.0, 2.0],
        "speed_final_m_per_s": 1.0e6,
        "energy_final_ev": 1.0e5,
    }


def _summary():
    return {
        "position_unit": "mm",
        "speed_unit": "m/s",
        "energy_unit": "eV",
        "magnetic_off": {
            "source_current_a": 1.0,
            "collision_current_a": 1.0,
            "collision_power_w": 1.0e5,
            "boundary_hit_count": 3,
            "trajectories": [_trajectory(i, float(i)) for i in range(3)],
        },
        "magnetic_on": {
            "source_current_a": 1.0,
            "collision_current_a": 1.0,
            "collision_power_w": 1.0e5 * (1.0 + 1.0e-8),
            "boundary_hit_count": 3,
            "trajectories": [_trajectory(i, float(i) + 1.0e-3) for i in range(3)],
        },
    }


def test_magnetic_trajectory_gate_accepts_deflection_without_work():
    result = evaluate_magnetic_trajectory_pair(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert json.loads(accelerator_magnetic_trajectory_pair_gate(json.dumps(_summary())))["status"] == "ok"


def test_magnetic_trajectory_gate_rejects_energy_and_current_drift():
    bad = copy.deepcopy(_summary())
    bad["magnetic_on"]["trajectories"][1]["energy_final_ev"] *= 1.02
    bad["magnetic_on"]["collision_current_a"] = 0.9
    result = evaluate_magnetic_trajectory_pair(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["magnetic_field_preserves_energy"] is False
    assert result["checks"]["source_collision_current_closes"] is False


def test_magnetic_trajectory_gate_rejects_no_deflection_and_id_mismatch():
    bad = copy.deepcopy(_summary())
    bad["magnetic_on"]["trajectories"] = copy.deepcopy(bad["magnetic_off"]["trajectories"])
    bad["magnetic_on"]["trajectories"][2]["particle_id"] = 99
    result = evaluate_magnetic_trajectory_pair(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["trajectory_ids_match"] is False
    assert result["checks"]["magnetic_field_changes_trajectory"] is False
