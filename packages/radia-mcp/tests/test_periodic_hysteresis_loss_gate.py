import copy
import json
import math

from radia_mcp.magnetic_materials.server import (
    periodic_hysteresis_loss_energy_gate,
)


def good() -> dict:
    rows = []
    for index in range(31):
        time_s = 0.2 * index
        startup = 1.0 + 0.2 * math.exp(-time_s) if time_s < 2.0 else 1.0
        power = startup * (15000.0 + 30000.0 * math.sin(math.pi * time_s))
        rows.append(
            {
                "time_s": time_s,
                "joule_loss_w": 0.0,
                "hysteresis_loss_w": power,
                "iron_loss_w": power,
                "hysteresis_part_w": power,
                "hysteresis_total_w": power,
            }
        )
    return {
        "contract": {
            "time_unit": "s",
            "power_unit": "W",
            "cycle_period_s": 2.0,
            "single_part_total": True,
        },
        "rows": rows,
    }


def test_accepts_periodic_positive_energy_with_negative_power_intervals():
    result = json.loads(periodic_hysteresis_loss_energy_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["checks"]["instantaneous_power_has_return_interval"] is True
    assert result["metrics"]["latest_cycle_energy_j"] > 0.0


def test_rejects_nonperiodic_final_cycle_and_negative_cycle_energy():
    payload = copy.deepcopy(good())
    for row in payload["rows"][-11:]:
        row["hysteresis_loss_w"] *= -1.0
        row["iron_loss_w"] = row["hysteresis_loss_w"]
        row["hysteresis_part_w"] = row["hysteresis_loss_w"]
        row["hysteresis_total_w"] = row["hysteresis_loss_w"]
    result = json.loads(periodic_hysteresis_loss_energy_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["steady_waveform_repeats"] is False
    assert result["checks"]["cycle_energies_are_positive"] is False


def test_rejects_decomposition_and_part_total_mismatch():
    payload = copy.deepcopy(good())
    payload["rows"][10]["iron_loss_w"] += 100.0
    payload["rows"][20]["hysteresis_total_w"] += 100.0
    result = json.loads(periodic_hysteresis_loss_energy_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["iron_loss_decomposition_closes"] is False
    assert result["checks"]["single_part_total_closes"] is False
