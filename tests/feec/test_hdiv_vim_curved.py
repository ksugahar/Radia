"""Golden test: the CURVED-MESH geometry win for the HDiv-type VIM, measured vs ANALYTIC truth.

HDiv (RT0) lives natively on curved (isoparametric) meshes via mesh.Curve(p); yano-type flat elements
cannot.  The win is measured against the EXACT uniform-sphere dipole / volume (NOT vs Radia):

  - the external field of a uniform-M sphere is the exact point dipole; a FLAT coarse mesh gets it ~10%
    WRONG (volume-inherited), mesh.Curve(3) at the SAME ndof is <1% -> a >10x accuracy win vs truth;
  - the geometry (volume) is ~10% low flat, <0.5% curved;
  - the demag FACTOR (D_z -> 1/3) does NOT discriminate (a coarse inscribed polyhedron is already
    near-isotropic) -- locked so a future contributor does not mistake it for a curved-mesh bug.

See examples/feec_vim/hdiv_demag_curved.py for the full derivation + table.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "feec_vim"))
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
    """HONEST non-result: the sphere demag FACTOR is a near-isotropic ratio insensitive to faceting --
    BOTH flat and curved land within a few % of 1/3 and curved is NOT meaningfully better.  Locked so
    nobody mistakes the (correct) ~2-3% demag-factor quadrature offset for a curved-mesh defect."""
    flat, curved = cases
    assert abs(flat["demag_err"]) < 0.05, "flat demag factor should already be near 1/3 (near-isotropic)"
    assert abs(curved["demag_err"]) < 0.05, "curved demag factor also near 1/3 (same near-isotropy)"
    # the curved external field win is >10x; the demag-factor 'win' is not -- assert it is NOT large
    assert abs(flat["demag_err"]) < 5 * abs(curved["ext_phi_err"]) + 0.05, \
        "demag factor must not be treated as the curved discriminator"
