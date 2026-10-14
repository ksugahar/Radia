import copy
import json

from radia_mcp.radia_ngsolve.loss_temperature_coupling_gate import loss_temperature_coupling_gate
from radia_mcp.radia_ngsolve.server import loss_temperature_coupling_gate as mcp_gate


TARGET_LOSS = [48618.25146117232, 43261.47171118402, 41289.46387123719, 39712.68803144368, 38998.80481709342]
AUX_LOSS = [10581.004756770177, 10633.959619230067, 10657.89742605562, 10672.403994127058, 10680.354683960557]
ACTIVE_POWER = [64138.719730656885, 58153.63022360937, 55864.84999756735, 54053.470026155446, 53176.329304105704]


def _magnetic():
    return [
        {
            "target_loss_w": target,
            "auxiliary_loss_w": auxiliary,
            "other_loss_w": 0.0,
            "total_loss_w": target + auxiliary,
            "active_input_power_w": active,
        }
        for target, auxiliary, active in zip(TARGET_LOSS, AUX_LOSS, ACTIVE_POWER)
    ]


def _thermal():
    temperatures = [20.0, 35.658026449388615, 49.14605889297928, 61.673723906569045, 73.50035392202635, 84.91956512101599]
    return [
        {
            "time_s": 0.2 * index,
            "target_heat_source_w": 0.0 if index == 0 else TARGET_LOSS[index - 1] / 96.0,
            "average_temperature_c": temperature,
        }
        for index, temperature in enumerate(temperatures)
    ]


def test_loss_temperature_coupling_accepts_scaled_live_shape_and_dispatches():
    result = loss_temperature_coupling_gate(_magnetic(), _thermal(), loss_to_heat_scale=1.0 / 96.0)
    assert result["status"] == "ok"
    assert result["metrics"]["time_step_s"] == 0.2
    assert json.loads(mcp_gate(_magnetic(), _thermal(), 1.0 / 96.0))["status"] == "ok"


def test_loss_temperature_coupling_rejects_stale_heat_row_and_nonmonotone_temperature():
    thermal = copy.deepcopy(_thermal())
    thermal[3]["target_heat_source_w"] *= 0.8
    thermal[4]["average_temperature_c"] = thermal[3]["average_temperature_c"] - 1.0
    result = loss_temperature_coupling_gate(_magnetic(), thermal, loss_to_heat_scale=1.0 / 96.0)
    assert result["status"] == "needs_attention"
    assert result["checks"]["loss_to_heat_scale_matches"] is False
    assert result["checks"]["average_temperature_strictly_increases"] is False
