# -*- coding: utf-8 -*-
"""Fast AGE-related algebra checks that belong in pytest."""

import math
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import integral_slot_winding_factor


def _winding_factor(p, q, n):
    q_slots = 2 * p * 3 * q
    gamma = 2 * math.pi * p / q_slots
    slots = []
    for i in range(q_slots):
        edeg = round(360 * p * i / q_slots) % 360
        if 0 <= edeg < 60:
            slots.append((2 * math.pi * i / q_slots, +1.0))
        elif 180 <= edeg < 240:
            slots.append((2 * math.pi * i / q_slots, -1.0))
    phase_a_slots = len(slots)
    kw = abs(sum(sign * np.exp(-1j * n * p * theta) for theta, sign in slots)) / phase_a_slots
    kd = abs(math.sin(n * q * gamma / 2) / (q * math.sin(n * gamma / 2)))
    return kw, kd, phase_a_slots


def test_distributed_winding_factor_matches_distribution_factor():
    worst = 0.0
    for p, q in [(1, 2), (2, 2), (3, 3)]:
        for n in [1, 3, 5, 7]:
            kw, kd, phase_a_slots = _winding_factor(p, q, n)
            summary = integral_slot_winding_factor(slots=2 * p * 3 * q, poles=2 * p,
                                                   harmonic=n, phases=3)
            assert abs(abs(summary["winding_factor"]) - kw) < 1e-12
            assert abs(abs(summary["distribution_factor"]) - kd) < 1e-12
            assert summary["slots_per_phase"] == phase_a_slots
            if n == 1:
                assert phase_a_slots == 2 * p * q
            worst = max(worst, abs(kw - kd) / max(kd, 1e-12))
    assert worst < 1e-9
