r"""Magnetizing (air-gap) inductance (#69) -- test.

L_mu = (2 mu0/pi)(kw1 N_ph)^2 D L_stk/(p^2 g_eff), L_m = (m/2) L_mu. Validated against the explicit
MMF->B->flux->linkage derivation chain, the 1/g_eff / 1/p^2 / (kw N)^2 scalings, the m/2 synchronous
factor, and the synthesis with winding factor (#51) + Carter g_eff (#57). Tool-independent."""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (winding_factor, effective_air_gap, MU0,
                                           magnetizing_inductance_per_phase,
                                           synchronous_magnetizing_inductance)

D, L_STK, P, N_PH = 0.10, 0.10, 2, 120
KW1 = winding_factor(1, 3, 360.0 * P / 36, pitch_fraction=5.0 / 6.0)
G_EFF = effective_air_gap(math.pi * D / 36, 0.8e-3, 2.5e-3)


def _chain(D, L, g_eff, p, kw1, N):
    tau_p = math.pi * D / (2 * p)
    B1 = MU0 * (2.0 / math.pi) * kw1 * N / p / g_eff
    return kw1 * N * (2.0 / math.pi) * B1 * tau_p * L


def test_matches_derivation_chain():
    for (Dx, Lx, p, N) in [(0.10, 0.10, 2, 120), (0.20, 0.15, 3, 80), (0.05, 0.04, 1, 200)]:
        L_mu = magnetizing_inductance_per_phase(Dx, Lx, G_EFF, p, KW1, N)
        assert math.isclose(L_mu, _chain(Dx, Lx, G_EFF, p, KW1, N), rel_tol=1e-12)


def test_scalings():
    base = magnetizing_inductance_per_phase(D, L_STK, G_EFF, P, KW1, N_PH)
    # 1/g_eff
    assert math.isclose(base / magnetizing_inductance_per_phase(D, L_STK, 2 * G_EFF, P, KW1, N_PH),
                        2.0, rel_tol=1e-12)
    # 1/p^2
    assert math.isclose(base / magnetizing_inductance_per_phase(D, L_STK, G_EFF, 2 * P, KW1, N_PH),
                        4.0, rel_tol=1e-12)
    # (kw1 N)^2
    assert math.isclose(base / magnetizing_inductance_per_phase(D, L_STK, G_EFF, P, KW1, N_PH // 2),
                        4.0, rel_tol=1e-12)
    # linear in D and in L_stk
    assert math.isclose(magnetizing_inductance_per_phase(2 * D, L_STK, G_EFF, P, KW1, N_PH), 2 * base,
                        rel_tol=1e-12)
    assert math.isclose(magnetizing_inductance_per_phase(D, 3 * L_STK, G_EFF, P, KW1, N_PH), 3 * base,
                        rel_tol=1e-12)


def test_synchronous_factor_and_carter():
    L_mu = magnetizing_inductance_per_phase(D, L_STK, G_EFF, P, KW1, N_PH)
    # synchronous = (m/2) * per-phase self; 3-phase -> 3/2
    assert math.isclose(synchronous_magnetizing_inductance(D, L_STK, G_EFF, P, KW1, N_PH),
                        1.5 * L_mu, rel_tol=1e-12)
    assert math.isclose(synchronous_magnetizing_inductance(D, L_STK, G_EFF, P, KW1, N_PH, phases=5),
                        2.5 * L_mu, rel_tol=1e-12)
    # Carter slotting lowers L_mu vs a smooth gap; kw1 is a sane fundamental factor
    L_smooth = magnetizing_inductance_per_phase(D, L_STK, 0.8e-3, P, KW1, N_PH)
    assert L_mu < L_smooth
    assert 0.9 < KW1 < 0.96
    # realistic magnitude: tens of mH for this small machine
    assert 5e-3 < L_mu < 100e-3


if __name__ == "__main__":
    test_matches_derivation_chain()
    test_scalings()
    test_synchronous_factor_and_carter()
    print("[OK] magnetizing inductance: L_mu == chain; 1/g_eff, 1/p^2, (kwN)^2 scalings; L_m=(3/2)L_mu.")
