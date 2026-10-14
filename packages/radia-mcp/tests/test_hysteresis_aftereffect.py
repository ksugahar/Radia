# -*- coding: utf-8 -*-
r"""Magnetic aftereffect / thermal viscosity for the B-input Play model -- physics-gated.

The rate/temperature-dependent side of HysterSoft (Dimian-Andrei, "Noise-Driven Phenomena in Hysteretic
Systems").  With B held, thermal activation relaxes the play states toward the held field via the
FLUCTUATION FIELD eta_f (= kT / activation-volume): the irreversible lag |p_k - B| shrinks by
eta_f*ln(t/t0), so H relaxes LOGARITHMICALLY, H(t) = H0 - S ln(t/t0), with magnetic viscosity
S = eta_f * sum_{still-relaxing cells} a_k -- proportional to eta_f and to the active differential
susceptibility (the link that ties the aftereffect to FORC, both HysterSoft features).
"""
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.hysteresis import PlayHysteresis


def _prepared(Bsat, Bhold):
    """Saturate +Bsat then descend to Bhold; return (model, px, py, Bhold)."""
    m = PlayHysteresis(eta=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], a=[150.0, 120.0, 100.0, 80.0, 60.0, 40.0])
    px = np.zeros(m.K); py = np.zeros(m.K)
    _, _, px, py = m.step(Bsat, 0.0, px, py)
    _, _, px, py = m.step(Bhold, 0.0, px, py)
    return m, px, py


def test_aftereffect_is_logarithmic_in_time():
    Bsat, Bhold = 1.4, 0.2
    m, px, py = _prepared(Bsat, Bhold)
    ln_t = np.linspace(0.0, 4.0, 40)                      # all dissipative cells still relaxing here
    for eta_f in (0.02, 0.04):
        Hx, Hy = m.aftereffect(Bhold, 0.0, px, py, eta_f, ln_t)
        p = np.polyfit(ln_t, Hx, 1)
        resid = Hx - np.polyval(p, ln_t)
        R2 = 1.0 - np.sum(resid ** 2) / np.sum((Hx - Hx.mean()) ** 2)
        S, S_ana = -p[0], m.magnetic_viscosity(eta_f)
        print(f"eta_f={eta_f}: R^2(H vs ln t)={R2:.5f}, viscosity S={S:.2f} vs eta_f*sum(a_k)={S_ana:.2f}")
        assert R2 > 0.999, "the aftereffect H(t) must be logarithmic in time"
        assert abs(S - S_ana) / S_ana < 0.03, "S = eta_f * sum_{eta_k>0} a_k (active susceptibility)"
        assert np.allclose(Hy, 0.0, atol=1e-9), "axial relaxation has no transverse H"


def test_viscosity_proportional_to_fluctuation_field():
    Bsat, Bhold = 1.4, 0.2
    m, px, py = _prepared(Bsat, Bhold)
    ln_t = np.linspace(0.0, 3.0, 40)                      # all cells still active for eta_f up to 0.06
    S = {}
    for eta_f in (0.02, 0.04, 0.06):
        Hx, _ = m.aftereffect(Bhold, 0.0, px, py, eta_f, ln_t)
        S[eta_f] = -np.polyfit(ln_t, Hx, 1)[0]
    assert abs(S[0.04] / S[0.02] - 2.0) < 0.03 and abs(S[0.06] / S[0.02] - 3.0) < 0.03, \
        "magnetic viscosity must be proportional to the fluctuation field eta_f"


def test_aftereffect_relaxes_toward_held_field():
    """Over long times every cell fully relaxes (lag -> 0): H -> the anhysteretic value sum a_k * B."""
    Bsat, Bhold = 1.4, 0.2
    m, px, py = _prepared(Bsat, Bhold)
    Hx, _ = m.aftereffect(Bhold, 0.0, px, py, eta_f=0.05, ln_times=np.array([0.0, 50.0]))
    H_anhyst = float(np.sum(m.a) * Bhold)
    assert abs(Hx[-1] - H_anhyst) < 1e-6, "fully relaxed H = sum(a_k) * B (anhysteretic)"
    assert Hx[0] != Hx[-1], "the state must actually relax"


def main():
    test_aftereffect_is_logarithmic_in_time()
    test_viscosity_proportional_to_fluctuation_field()
    test_aftereffect_relaxes_toward_held_field()
    print("[OK] magnetic aftereffect / thermal viscosity: H(t)=H0 - S ln t (logarithmic), "
          "S = eta_f * sum(active a_k) proportional to the fluctuation field, relaxing to the "
          "anhysteretic sum(a_k)*B -- HysterSoft's noise-driven side on the Play model.")


if __name__ == "__main__":
    main()
