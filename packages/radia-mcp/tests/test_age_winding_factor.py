# -*- coding: utf-8 -*-
"""Fast AGE-related algebra checks that belong in pytest."""

import math

import numpy as np


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
            if n == 1:
                assert phase_a_slots == 2 * p * q
            worst = max(worst, abs(kw - kd) / max(kd, 1e-12))
    assert worst < 1e-9
