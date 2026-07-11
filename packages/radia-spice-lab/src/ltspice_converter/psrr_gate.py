"""Solver-neutral transient power-supply rejection gate."""
import math


def transient_psrr_gate(summary: dict, minimum_psrr_db: float = 40.0) -> dict:
    names=("input_ripple_vpp","output_ripple_vpp","measured_psrr_db","raw_input_ripple_vpp","raw_output_ripple_vpp","frequency_hz","window_s")
    if not isinstance(summary,dict) or any(name not in summary for name in names): raise ValueError("complete PSRR summary required")
    vin=float(summary["input_ripple_vpp"]); vout=float(summary["output_ripple_vpp"]); rvin=float(summary["raw_input_ripple_vpp"]); rvout=float(summary["raw_output_ripple_vpp"]); measured=float(summary["measured_psrr_db"]); freq=float(summary["frequency_hz"]); window=float(summary["window_s"])
    if any(not math.isfinite(x) or x<=0 for x in (vin,vout,rvin,rvout,freq,window)): raise ValueError("ripple, frequency, and window values must be positive and finite")
    derived=20*math.log10(vin/vout); raw=20*math.log10(rvin/rvout)
    checks={"output_ripple_attenuated":vout<vin,"measure_formula_closes":abs(measured-derived)<=.02,"raw_replay_matches_measures":abs(raw-measured)<=.02 and abs(rvout-vout)/vout<=.002,"window_covers_many_cycles":window*freq>=100,"psrr_above_minimum":measured>=float(minimum_psrr_db)}
    return {"policy":"transient_psrr_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[k for k,v in checks.items() if not v],"metrics":{"derived_psrr_db":derived,"raw_psrr_db":raw,"cycles_in_window":window*freq}}
