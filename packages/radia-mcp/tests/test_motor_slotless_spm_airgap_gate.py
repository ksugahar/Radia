"""Guards: motor_slotless_spm_airgap_field_gate forwards to the Radia closed form and fails loudly.

Assumes: the closed form itself is guarded by tests/analytical_formulas in radia.
"""
import json

import numpy as np

from radia_mcp.motor.server import motor_slotless_spm_airgap_field_gate

GEO = dict(pole_pairs=2, remanence_t=1.2, recoil_permeability=1.05, rotor_radius_m=0.025,
           magnet_radius_m=0.03, stator_radius_m=0.032, pole_arc_fraction=0.8, radius_m=0.031)
TH = list(np.linspace(0.0, np.pi / 2, 19))


def _gate(br, bt, **kw):
    return json.loads(motor_slotless_spm_airgap_field_gate(**{**GEO, **kw}, theta_rad=TH,
                                                           computed_b_r_t=br, computed_b_theta_t=bt))


def test_exact_passes_shifted_and_wrong_sign_fail():
    from radia.analytical_formulas.spm_slotless_airgap import slotless_spm_airgap_field
    g = {k: v for k, v in GEO.items() if k != "radius_m"}
    br, bt = slotless_spm_airgap_field(GEO["radius_m"], TH, g["pole_pairs"], g["remanence_t"],
                                       g["recoil_permeability"], g["rotor_radius_m"],
                                       g["magnet_radius_m"], g["stator_radius_m"], g["pole_arc_fraction"])
    ok = _gate(br.tolist(), bt.tolist())
    assert ok["status"] == "ok" and ok["relative_difference"] == 0.0, ok
    assert _gate((1.03 * br).tolist(), bt.tolist())["status"] == "needs_attention"
    assert _gate(br.tolist(), (-bt).tolist())["status"] == "needs_attention"


def test_invalid_inputs_are_reported():
    assert _gate([0.0] * 3, [0.0] * 19)["status"] == "invalid_input"
    assert _gate([0.0] * 19, [0.0] * 19, pole_pairs=1)["status"] == "invalid_input"
    assert _gate([0.0] * 19, [0.0] * 19, radius_m=0.02)["status"] == "invalid_input"
