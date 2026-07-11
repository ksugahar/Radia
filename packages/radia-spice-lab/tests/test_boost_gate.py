from ltspice_converter.boost_gate import boost_converter_steady_state_gate
from ltspice_converter.mcp_server import boost_converter_steady_state_gate as mcp_gate
def good(): return {"metrics":{"input_voltage_v":2,"output_average_v":5,"output_ripple_vpp":.1,"input_power_delivered_w":10,"load_power_w":8,"inductor_voltage_average_v":.001,"capacitor_current_average_a":1e-5,"charge_balance_limit_a":1e-3,"tail_previous_relative_drift":1e-4}}
def test_ok_and_dispatch(): assert boost_converter_steady_state_gate(good())["status"]=="ok" and mcp_gate(good())["status"]=="ok"
def test_rejects_active_and_drifting():
 row=good(); row["metrics"]["load_power_w"]=12; row["metrics"]["tail_previous_relative_drift"]=.2
 result=boost_converter_steady_state_gate(row); assert result["status"]=="needs_attention"
