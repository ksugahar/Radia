"""Solver-neutral switched boost-converter steady-state gate."""
def boost_converter_steady_state_gate(summary: dict) -> dict:
    m=summary.get("metrics") or {}; req=["input_voltage_v","output_average_v","output_ripple_vpp","input_power_delivered_w","load_power_w","inductor_voltage_average_v","capacitor_current_average_a","tail_previous_relative_drift"]
    try: v={k:float(m[k]) for k in req}
    except (KeyError,TypeError,ValueError) as exc: raise ValueError("missing boost metrics") from exc
    eff=v["load_power_w"]/v["input_power_delivered_w"] if v["input_power_delivered_w"]>0 else -1
    checks={"boosts_voltage":v["output_average_v"]>1.5*v["input_voltage_v"],"tail_is_periodic":v["tail_previous_relative_drift"]<.01,
            "ripple_is_bounded":0<v["output_ripple_vpp"]<.5*v["output_average_v"],"power_is_passive":0<eff<1,
            "inductor_volt_second_balances":abs(v["inductor_voltage_average_v"])<.02*v["input_voltage_v"],
            "capacitor_charge_balances":abs(v["capacitor_current_average_a"])<float(m.get("charge_balance_limit_a",1e-3))}
    return {"policy":"boost_converter_steady_state_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[k for k,x in checks.items() if not x],"metrics":{"boost_ratio":v["output_average_v"]/v["input_voltage_v"],"efficiency":eff},"lesson":"Gate a switched converter only after periodic steady state, power passivity, inductor volt-second balance, and capacitor charge balance all close."}
