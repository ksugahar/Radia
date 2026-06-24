"""Golden test: the CURVED-MESH geometry win for the HDiv-type VIM, measured vs ANALYTIC truth.

HDiv (RT0) lives natively on curved (isoparametric) meshes via mesh.Curve(p); six-face surface-charge flat elements
cannot.  The win is measured against the EXACT uniform-sphere dipole / volume (NOT vs Radia):

  - the external field of a uniform-M sphere is the exact point dipole; a FLAT coarse mesh gets it ~10%
    WRONG (volume-inherited), mesh.Curve(3) at the SAME ndof is <1% -> a >10x accuracy win vs truth;
  - the geometry (volume) is ~10% low flat, <0.5% curved;
  - with THIS elementary method's crude sub-point Gram the demag FACTOR does NOT cleanly discriminate
    (its ~2% quadrature bias masks the ~0.25% geometry signal) -- a limitation of the crude Gram, NOT
    of the demag factor: the proper ngsolve.bem single-layer (test_hdiv_vim_bem_demag.py) DOES
    discriminate + p-converges.  Locked so the crude offset is not mistaken for a curved-mesh bug.

See examples/vim/hdiv_demag_curved.py for the full derivation + table.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
import hdiv_demag_curved as cv  # noqa: E402


@pytest.fixture(scope="module")
def cases():
    res = cv.run(hs=(0.8,), curve_order=3)
    flat = next(c for c in res["cases"] if not c["curved"])
    curved = next(c for c in res["cases"] if c["curved"])
    return flat, curved


def test_curved_external_field_matches_dipole_truth(cases):
    """The headline win: curved external field is <1% vs the exact dipole; flat is >5% off."""
    flat, curved = cases
    assert abs(curved["ext_phi_err"]) < 0.01, \
        f"curved ext field {100*curved['ext_phi_err']:.2f}% not within 1% of analytic dipole"
    assert abs(flat["ext_phi_err"]) > 0.05, \
        f"flat ext field {100*flat['ext_phi_err']:.2f}% unexpectedly accurate (win premise broken)"
    assert abs(curved["ext_phi_err"]) < 0.2 * abs(flat["ext_phi_err"]), \
        "curved external field should be >5x more accurate than flat at the same ndof"
    assert flat["n_bnd"] == curved["n_bnd"], "flat/curved must be the SAME mesh (same ndof)"


def test_curved_geometry_exact(cases):
    """Volume + area: ~10% low flat, <0.5% curved (the isoparametric geometry is exact to order)."""
    flat, curved = cases
    assert abs(curved["V_err"]) < 5e-3, f"curved V err {100*curved['V_err']:.2f}% not <0.5%"
    assert flat["V_err"] < -0.05, f"flat V err {100*flat['V_err']:.2f}% should be <-5% (faceting)"
    assert abs(curved["area_err"]) < 0.01, f"curved area err {100*curved['area_err']:.2f}% not <1%"


def test_demag_factor_does_not_discriminate(cases):
    """The crude sub-point Gram's ~2% quadrature bias masks the ~0.25% geometry signal, so with THIS
    elementary method both flat and curved land within a few % of 1/3 and the demag factor is NOT a
    clean discriminator (use ext_potential, or the proper ngsolve.bem single-layer in
    test_hdiv_vim_bem_demag.py which DOES discriminate + p-converges).  Locked so the crude offset is
    not mistaken for a curved-mesh defect -- it is a Gram-quadrature limitation, not near-isotropy."""
    flat, curved = cases
    assert abs(flat["demag_err"]) < 0.05, "flat demag factor should already be near 1/3 (near-isotropic)"
    assert abs(curved["demag_err"]) < 0.05, "curved demag factor also near 1/3 (same near-isotropy)"
    # the curved external field win is >10x; the demag-factor 'win' is not -- assert it is NOT large
    assert abs(flat["demag_err"]) < 5 * abs(curved["ext_phi_err"]) + 0.05, \
        "demag factor must not be treated as the curved discriminator"
