"""Solver-neutral half-wave capacitor-input rectifier gate."""
def half_wave_rectifier_gate(vin_peak_v, frequency_hz, capacitance_f, load_ohm, vout_avg_v, vout_pp_v, diode_avg_a):
    vals=[float(x) for x in (vin_peak_v,frequency_hz,capacitance_f,load_ohm,vout_avg_v,vout_pp_v,diode_avg_a)]
    if any(x<=0 for x in vals): raise ValueError("all rectifier metrics must be positive")
    vin,f,c,r,avg,pp,id=vals; il=avg/r; ce=abs(id-il)/max(id,il); est=il/(f*c); re=abs(pp-est)/est
    checks={"dc_below_input_peak":avg<vin,"load_and_diode_average_close":ce<=.02,"capacitor_discharge_ripple_order_correct":re<=.35,"output_remains_positive":avg-pp/2>0}
    return {"policy":"half_wave_rectifier_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[k for k,v in checks.items() if not v],"metrics":{"load_average_a":il,"current_relative_error":ce,"first_order_ripple_vpp":est,"ripple_relative_error":re}}
