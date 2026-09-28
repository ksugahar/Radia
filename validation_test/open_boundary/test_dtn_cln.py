# -*- coding: utf-8 -*-
"""Golden tests for the Cauer realisation of the eddy DtN (radia.open_boundary.dtn_cln).

Scheduled for removal with the repository's Cauer-ladder code; the exact DtN
symbols are tested in test_dtn_exact.py.
"""
import numpy as np
import pytest

import radia.open_boundary as ob

OMEGA = np.logspace(-1, 2, 60)
MODES = (1, 2, 3, 4, 5, 6)


def _nrmse(a, b):
    return float(np.sqrt(np.mean(np.abs(a - b) ** 2)) / np.sqrt(np.mean(np.abs(b) ** 2)))


@pytest.mark.parametrize("n", MODES)
def test_cauer_exact_at_n_plus_1_stages(n):
    """The Cauer-in-q ladder is EXACT at exactly n+1 stages, well-conditioned."""
    stages = ob.cauer_ladder(n)
    assert len(stages) == n + 1, f"n={n}: expected n+1 stages, got {len(stages)}"
    Zc = np.array([ob.eval_ladder(stages, 1j * w) for w in OMEGA])
    Zr = np.array([ob.eddy_dtn(n, 1j * w) for w in OMEGA])
    assert _nrmse(Zc, Zr) < 1e-10, f"n={n}: Cauer ladder not exact vs symbol"
    allc = np.abs(np.concatenate([np.asarray(s, float) for s in stages]))
    spread = float(np.max(allc) / np.min(allc[allc > 0]))
    assert spread < 1e3, f"n={n}: ladder ill-conditioned (spread {spread:.1e})"


@pytest.mark.parametrize("R0,mu_sigma", [(1.0, 1.0), (0.1, 0.5), (0.03, 4.0e-7 * np.pi * 5.8e7)])
def test_ladder_matches_symbol_general_units(R0, mu_sigma):
    """eval_ladder reproduces eddy_dtn for non-unit (R0, mu_sigma)."""
    for n in (1, 2, 3):
        stages = ob.cauer_ladder(n)
        Zc = np.array([ob.eval_ladder(stages, 1j * w, R0, mu_sigma) for w in OMEGA])
        Zr = np.array([ob.eddy_dtn(n, 1j * w, R0, mu_sigma) for w in OMEGA])
        assert _nrmse(Zc, Zr) < 1e-9, f"n={n} R0={R0} mu_sigma={mu_sigma}: ladder != symbol"
