from radia.ltspice.rectifier_gate import half_wave_rectifier_gate
def test_accepts(): assert half_wave_rectifier_gate(20,1e5,200e-9,100,15.0133,6.0698,.150133)["status"]=="ok"
def test_rejects(): assert half_wave_rectifier_gate(20,1e5,200e-9,100,15,1,.3)["status"]=="needs_attention"
