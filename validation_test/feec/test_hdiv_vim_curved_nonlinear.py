"""Golden test + HONEST FINDING: curved x nonlinear on a spheroid (verify-first, 2026-06-08).

A spheroid in a uniform field magnetizes UNIFORMLY for ANY M-H law (linear or nonlinear) -- so the
nonlinear magnetization is the scalar fixed point  M = M(H_ext - N_par * M)  with N_par the (curved)
demag factor.  Two things are locked here:

  1. CAPABILITY (the curved nonlinear solve is CORRECT): with N_par from the curved + order-2 single-
     layer, the nonlinear M matches the analytic-demag fixed point to <0.05% across weak->saturation.

  2. THE HONEST MAGNITUDE (the curved win on M is MODEST, ~0.3%, NOT dramatic): the demag factor is a
     VOLUME-NORMALIZED RATIO, so it cancels the ~10% volume faceting error -> the flat nonlinear M is
     only ~0.3% off (and this does NOT grow with severe curvature -- thin disks facet *better* on the
     demag).  The large curved win (~10%) lives in the FIELD output (volume error, not cancelled; see
     hdiv_demag_curved.py external-field probe), NOT in the magnetization.

So this test deliberately asserts the flat error is SMALL -- documenting that curved x nonlinear is a
real six-face surface-charge-impossible *capability* whose accuracy benefit on the *magnetization* is modest.  (Radia
is not a clean reference for curved geometry -- its ObjHex/Tet facet the body -- so the validatable
curved nonlinear case is the analytic spheroid; the genuinely non-uniform volume-charge case has no
clean curved reference.)
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("ngsolve.bem")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "vim_legacy"))
import hdiv_demag_bem_singlelayer as bem  # noqa: E402
from radia.vim import _nonlinear as nl     # noqa: E402

_CHI0, _MSAT = 5000.0, 1.6e6              # soft-iron-like saturating law
def _Mof(H):
    return _CHI0 * H / (1.0 + _CHI0 * abs(H) / _MSAT)


def test_curved_nonlinear_spheroid_matches_analytic():
    """curved + order-2 demag -> nonlinear M matches the analytic-demag fixed point to <0.05%."""
    c = 2.0
    N_an = bem.spheroid_Nz_analytic(c)
    N_cv, _ = bem.demag_spheroid(c, 0.5, 3, 2, axis=2)
    for H_ext in (1e3, 1e4, 1e5):
        M_an = nl._scalar_fixed_point(_Mof, N_an, H_ext)
        M_cv = nl._scalar_fixed_point(_Mof, N_cv, H_ext)
        assert abs(M_cv / M_an - 1) < 5e-4, \
            f"curved nonlinear M {M_cv:.1f} vs analytic {M_an:.1f} at H={H_ext:.0e}"


def test_curved_nonlinear_M_win_is_modest_not_dramatic():
    """The HONEST magnitude: the flat nonlinear M is only MODESTLY off (the demag-ratio cancels the
    volume faceting error) -- a real but small win, NOT a dramatic differentiator on the magnetization.
    Locked so the curved x nonlinear win is not overclaimed."""
    c = 2.0
    N_an = bem.spheroid_Nz_analytic(c)
    N_cv, _ = bem.demag_spheroid(c, 0.5, 3, 2, axis=2)
    N_fl, _ = bem.demag_spheroid(c, 0.5, 0, 2, axis=2)
    M_an = nl._scalar_fixed_point(_Mof, N_an, 1e4)
    curved_err = abs(nl._scalar_fixed_point(_Mof, N_cv, 1e4) / M_an - 1)
    flat_err = abs(nl._scalar_fixed_point(_Mof, N_fl, 1e4) / M_an - 1)
    assert curved_err < 5e-4, "curved must be exact on the nonlinear M"
    # the win is real (curved beats flat) but MODEST (flat off by < 1%, not a dramatic differentiator)
    assert curved_err < flat_err, "curved should beat flat"
    assert flat_err < 0.01, f"flat nonlinear M error {100*flat_err:.3f}% -- documents the MODEST win (demag-ratio robust)"
