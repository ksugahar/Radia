"""Guards: maglev_cuboid_magnet_force_gate forwards to the Radia closed form and fails loudly.

Assumes: radia native module (ObjRecMag, FldFrc) for the independent Maxwell-stress route.
"""
import json

import numpy as np
import pytest

from radia_mcp.maglev.server import maglev_cuboid_magnet_force_gate

MU0 = 4e-7 * np.pi
S1, S2, D = [0.01, 0.01, 0.01], [0.01, 0.02, 0.005], [0.004, 0.003, 0.012]


def _gate(force, tol=5e-2, d=D):
    return json.loads(maglev_cuboid_magnet_force_gate(1.0, S1, 1.0, S2, d, force, tol))


def test_radia_maxwell_stress_route_passes_against_closed_form():
    rad = pytest.importorskip("radia")
    rad.UtiDelAll()
    m1 = rad.ObjRecMag([0, 0, 0], S1, [0, 0, 1 / MU0])
    m2 = rad.ObjRecMag(D, S2, [0, 0, 1 / MU0])
    plane = rad.FldFrcShpRtg([0, 0, 0.007], [1.0, 1.0])   # between the magnets, normal +z
    force_on_2 = -np.asarray(rad.FldFrc(rad.ObjCnt([m1, m2]), plane))
    rad.UtiDelAll()
    res = _gate(force_on_2.tolist(), tol=2e-2)
    assert res["status"] == "ok", res
    assert res["relative_difference"] < 2e-2, (res["computed_force_n"], res["reference_force_n"])


def test_wrong_force_needs_attention_and_bad_input_is_invalid():
    ref = np.asarray(_gate([0.0, 0.0, 1.0])["reference_force_n"])
    assert _gate((-ref).tolist())["status"] == "needs_attention"
    assert _gate((ref * 1.1).tolist())["status"] == "needs_attention"
    assert _gate(ref.tolist())["relative_difference"] == 0.0
    assert _gate(ref.tolist(), d=[0.002, 0.0, 0.0])["status"] == "invalid_input"
    assert _gate([1.0, float("nan"), 0.0])["status"] == "invalid_input"
