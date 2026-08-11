from radia.ltspice.mcp_server import cockcroft_walton_stage_gate as mcp_gate
from radia.ltspice.voltage_multiplier_gate import cockcroft_walton_stage_gate


def _args():
    return dict(
        vin_peak_v=9.9999961853,
        stage1_avg_v=18.7194127452,
        stage2_avg_v=37.45995233,
        stage2_previous_avg_v=37.4292950183,
        stage1_ripple_vpp=0.0225086212158,
        stage2_ripple_vpp=0.0301704406738,
        load_ohm=1e6,
        load_avg_a=3.74599523391e-5,
        source_power_delivered_w=0.00319411230065,
        load_power_w=0.00140324809082,
    )


def test_accepts_live_two_stage_multiplier_measurements():
    result = cockcroft_walton_stage_gate(**_args())
    assert result["status"] == "ok"
    assert mcp_gate(**_args())["status"] == "ok"


def test_rejects_unsettled_or_collapsed_second_stage():
    result = cockcroft_walton_stage_gate(
        **{**_args(), "stage2_avg_v": 25.0, "stage2_previous_avg_v": 20.0}
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["two_stage_voltage_law"] is False
    assert result["checks"]["late_windows_settled"] is False


def test_rejects_impossible_real_power_gain():
    result = cockcroft_walton_stage_gate(
        **{**_args(), "source_power_delivered_w": 0.001, "load_power_w": 0.002}
    )
    assert result["checks"]["source_power_covers_load"] is False
    assert result["checks"]["efficiency_physical"] is False
