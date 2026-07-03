"""Golden test: (A) the FIELD output of a curved x nonlinear body -- where the curved win is LARGE.

Unlike the magnetization (curved win modest ~0.3%, test_hdiv_vim_curved_nonlinear.py), the external
FIELD of a nonlinear soft-iron sphere inherits the ~9% volume faceting error directly (the dipole moment
m = M V), so the flat field is ~9% wrong at every external point while curved + order-3 is <0.5% -- a
~20x win, validated vs the ANALYTIC dipole (Radia cannot referee curved geometry -- it facets).

This locks where curved x nonlinear genuinely pays off: the field (the engineering deliverable).
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "vim_legacy"))
import hdiv_curved_nonlinear_field as cnf  # noqa: E402


@pytest.fixture(scope="module")
def cases():
    res = cnf.run(h=0.6, H_ext=1e4, curve_order=3)
    flat = next(c for c in res["cases"] if not c["curved"])
    curv = next(c for c in res["cases"] if c["curved"])
    return flat, curv


def test_curved_field_matches_dipole(cases):
    """curved + order-3 external field is <0.5% vs the analytic dipole at every point."""
    _, curv = cases
    assert curv["max_err"] < 5e-3, f"curved field max err {100*curv['max_err']:.2f}% not <0.5%"


def test_flat_field_is_wrong_by_volume_error(cases):
    """flat external field is ~9% off (it inherits the volume faceting error via the dipole moment)."""
    flat, _ = cases
    assert flat["max_err"] > 0.05, f"flat field max err {100*flat['max_err']:.2f}% should be >5%"


def test_curved_field_win_is_large(cases):
    """The headline: the curved field win is >10x (unlike the modest ~0.3% magnetization win)."""
    flat, curv = cases
    assert curv["max_err"] < 0.1 * flat["max_err"], \
        f"curved field should be >10x better than flat: flat {100*flat['max_err']:.2f}% vs curved {100*curv['max_err']:.2f}%"
