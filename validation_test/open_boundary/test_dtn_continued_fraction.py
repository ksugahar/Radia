# -*- coding: utf-8 -*-
"""Golden tests for the continued-fraction form of the eddy DtN
(radia.open_boundary.dtn_continued_fraction); the exact DtN symbols are tested
in test_dtn_exact.py.
"""
import numpy as np
import pytest

import radia.open_boundary as ob

OMEGA = np.logspace(-1, 2, 60)
MODES = (1, 2, 3, 4, 5, 6)


def _nrmse(a, b):
    return float(np.sqrt(np.mean(np.abs(a - b) ** 2)) / np.sqrt(np.mean(np.abs(b) ** 2)))


@pytest.mark.parametrize("n", MODES)
def test_continued_fraction_exact_at_n_plus_1_stages(n):
    """The continued fraction in q terminates after n+1 quotients, well-conditioned."""
    stages = ob.continued_fraction_stages(n)
    assert len(stages) == n + 1, f"n={n}: expected n+1 stages, got {len(stages)}"
    Zc = np.array([ob.eval_continued_fraction(stages, 1j * w) for w in OMEGA])
    Zr = np.array([ob.eddy_dtn(n, 1j * w) for w in OMEGA])
    assert _nrmse(Zc, Zr) < 1e-10, f"n={n}: continued fraction not exact vs symbol"
    allc = np.abs(np.concatenate([np.asarray(s, float) for s in stages]))
    spread = float(np.max(allc) / np.min(allc[allc > 0]))
    assert spread < 1e3, f"n={n}: continued fraction ill-conditioned (spread {spread:.1e})"


@pytest.mark.parametrize("R0,mu_sigma", [(1.0, 1.0), (0.1, 0.5), (0.03, 4.0e-7 * np.pi * 5.8e7)])
def test_continued_fraction_matches_symbol_general_units(R0, mu_sigma):
    """eval_continued_fraction reproduces eddy_dtn for non-unit (R0, mu_sigma)."""
    for n in (1, 2, 3):
        stages = ob.continued_fraction_stages(n)
        Zc = np.array([ob.eval_continued_fraction(stages, 1j * w, R0, mu_sigma) for w in OMEGA])
        Zr = np.array([ob.eddy_dtn(n, 1j * w, R0, mu_sigma) for w in OMEGA])
        assert _nrmse(Zc, Zr) < 1e-9, f"n={n} R0={R0} mu_sigma={mu_sigma}: continued fraction != symbol"
