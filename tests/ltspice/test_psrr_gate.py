from radia.ltspice.mcp_server import transient_psrr_gate
def good(): return {"input_ripple_vpp":.2,"output_ripple_vpp":35.3e-6,"measured_psrr_db":75.07,"raw_input_ripple_vpp":.2,"raw_output_ripple_vpp":35.3e-6,"frequency_hz":1e6,"window_s":5e-4}
def test_accepts_transient_psrr(): assert transient_psrr_gate(good())["status"]=="ok"
def test_rejects_stale_raw_and_short_window():
 r=good();r["raw_output_ripple_vpp"]*=2;r["window_s"]=1e-6;out=transient_psrr_gate(r);assert out["status"]=="needs_attention"
