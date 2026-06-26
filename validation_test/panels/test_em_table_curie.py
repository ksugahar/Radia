"""Unit tests for the opt-in Curie mu(T) collapse helpers in
calc_em_table (curie_factor / curie_scaled_bh).  Pure logic, no ESIM
solve -- fast.  The end-to-end Z_s(H,T) collapse is exercised by the
em_table golden.
"""
from __future__ import annotations

import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PANELS = os.path.join(os.path.dirname(os.path.dirname(HERE)),
                      "src", "radia", "panels")
if PANELS not in sys.path:
    sys.path.insert(0, PANELS)

import calc_em_table as em  # noqa: E402

MU0 = 4.0e-7 * math.pi


def test_curie_factor_endpoints_clamp_monotone():
    Tc = 770.0
    assert em.curie_factor(0.0, Tc) == pytest.approx(1.0)     # f(0) = 1
    assert em.curie_factor(Tc, Tc) == 0.0                     # f(Tc) = 0
    assert em.curie_factor(2 * Tc, Tc) == 0.0                 # above Tc -> 0
    # below 0 C: x = 1 - T/Tc > 1 -> clamped to 1
    assert em.curie_factor(-100.0, Tc) == pytest.approx(1.0)
    # disabled (Tc <= 0) -> always 1 (no mu(T))
    assert em.curie_factor(500.0, 0.0) == 1.0
    fs = [em.curie_factor(float(T), Tc) for T in range(0, 850, 50)]
    assert all(0.0 <= f <= 1.0 for f in fs)                   # in [0,1]
    assert all(fs[i + 1] <= fs[i] + 1e-12
               for i in range(len(fs) - 1))                   # non-increasing


def test_curie_scaled_bh_blend():
    bh = [[0.0, 0.0], [1000.0, 1.5], [1e4, 2.0], [1e5, 2.4]]
    # f = 1 -> unchanged
    for (H0, B0), (H1, B1) in zip(bh, em.curie_scaled_bh(bh, 1.0)):
        assert H1 == pytest.approx(H0) and B1 == pytest.approx(B0)
    # f = 0 -> paramagnetic line B = mu_0 * H (mu_r = 1)
    for H, B in em.curie_scaled_bh(bh, 0.0):
        assert B == pytest.approx(MU0 * H, abs=1e-12)
    # f = 0.5 -> convex blend, monotone non-decreasing in H
    sh = em.curie_scaled_bh(bh, 0.5)
    Bs = [B for _, B in sh]
    assert all(Bs[i + 1] >= Bs[i] - 1e-12 for i in range(len(Bs) - 1))
    # blend value at H = 1e4: mu_0*H + 0.5*(2.0 - mu_0*H)
    H, B = sh[2]
    assert B == pytest.approx(MU0 * 1e4 + 0.5 * (2.0 - MU0 * 1e4), rel=1e-9)
