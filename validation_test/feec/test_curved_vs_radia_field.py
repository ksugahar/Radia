"""Golden test: (B) head-to-head -- HDiv-VIM (CURVED) external field beats the SHIPPED Radia solver
(FLAT tets) on accuracy-per-resolution, vs the ANALYTIC dipole.

Locks the quantitative curved accuracy-per-DOF win against the production code: at the same mesh h the
curved field is >10x more accurate, and curved at a COARSE mesh beats shipped-Radia-flat at a FINER one.
(Accuracy-per-resolution, geometry-driven; wall-clock is the C++ productionization, not claimed.  Radia
is the accessible FLAT stand-in for the also-flat six-face surface-charge; the reference is the analytic dipole because
Radia cannot referee curved geometry -- it facets.)
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
pytest.importorskip("radia")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "vim_legacy"))
import compare_curved_vs_radia_field as cmp  # noqa: E402


def test_curved_beats_radia_flat_per_resolution():
    """At the same h, curved field is >10x more accurate than shipped Radia (flat); and curved at a
    COARSE mesh beats Radia-flat at a FINER mesh -- the accuracy-per-resolution win vs the production
    solver."""
    ne_c, flat_coarse = cmp.radia_flat(0.6)
    nb_c, curved_coarse = cmp.hdiv_curved(0.6)
    ne_f, flat_fine = cmp.radia_flat(0.3)
    # same coarse mesh: curved >10x more accurate than shipped Radia flat
    assert curved_coarse < 0.1 * flat_coarse, \
        f"curved {100*curved_coarse:.3f}% should be >10x better than Radia-flat {100*flat_coarse:.2f}% at h=0.6"
    # curved at COARSE beats shipped Radia flat at FINER (the accuracy-per-resolution win)
    assert curved_coarse < flat_fine, \
        f"curved@h0.6 {100*curved_coarse:.3f}% should beat Radia-flat@h0.3 {100*flat_fine:.2f}%"
    # and the shipped Radia flat error is genuinely large at coarse (sanity: the win premise holds)
    assert flat_coarse > 0.05, f"Radia-flat@h0.6 {100*flat_coarse:.2f}% should be >5% (faceting)"
