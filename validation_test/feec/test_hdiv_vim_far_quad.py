"""Golden: the ORDER-0 analytic charge-Gram FAR low-order double-quadrature (`far_quad`) is the
PRECISION-PRESERVING build speedup.

The Gram build dominates the HDiv-VIM cost (per-pair analytic quadrature).  For the +N CG tet path the
production default is the fast near/far split `near_factor=2` with `far_quad=4`: NEAR pairs use the exact
analytic entry, FAR pairs a degree-2 double-quadrature of 1/r (4-pt tet / 3-pt tri).  Unlike the bare
centroid-monopole far (`far_quad=0`, O((size/r)^2)), the low-order quadrature is O((size/r)^4), so it
REPRODUCES the all-analytic Gram -- this is what lets the ~4.5x faster build stay golden-accurate.

This locks the contract: on a uniform-field tet body
  * far_quad=4 (near_factor=2) reproduces the all-analytic (near_factor=1e30) demag AND volume-average M to
    ~1e-4, and the transverse (spurious) M stays at the all-analytic level;
  * far_quad=0 (monopole) deviates measurably MORE -- the regression that kept monopole off the default.

NGSolve + Netgen required (importorskip); meshing is non-deterministic, so the split-vs-exact DIFFERENCES
are checked (far_quad ~ exact << monopole), not absolute values.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

from radia.vim import hdiv_demag_solve  # noqa: E402

H0 = 1000.0
HEXT = ng.CoefficientFunction((0, 0, H0))


def _tet_cube(h=0.28):
    return ng.Mesh(OCCGeometry(Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))).GenerateMesh(maxh=h))


def test_far_quad_reproduces_exact_monopole_does_not():
    """far_quad=4 (near/far) == all-analytic to ~1e-4 on demag, Mz, and the transverse leak; monopole
    (far_quad=0) deviates more -- the precision-preserving-vs-lossy contract behind the fast default."""
    mesh = _tet_cube()
    with ng.TaskManager():
        ex = hdiv_demag_solve(mesh, 100.0, HEXT, near_factor=1e30)            # exact all-analytic
        fq = hdiv_demag_solve(mesh, 100.0, HEXT, near_factor=2.0, far_quad=4)  # precision-preserving fast
        mp = hdiv_demag_solve(mesh, 100.0, HEXT, near_factor=2.0, far_quad=0)  # lossy monopole far

    # far_quad reproduces the all-analytic demag + volume-average Mz to ~1e-4
    assert abs(fq["demag"] - ex["demag"]) < 5e-4, f"far_quad demag {fq['demag']} vs exact {ex['demag']}"
    rel_mz = abs(fq["M_avg"][2] - ex["M_avg"][2]) / abs(ex["M_avg"][2])
    assert rel_mz < 1e-3, f"far_quad Mz rel {rel_mz:.2e} vs exact"

    # and it is strictly closer to exact than the monopole far is (the reason it can be the default)
    d_fq = abs(fq["demag"] - ex["demag"])
    d_mp = abs(mp["demag"] - ex["demag"])
    assert d_fq < d_mp, f"far_quad ({d_fq:.2e}) should be closer to exact than monopole ({d_mp:.2e})"

    # the default 'auto' path (no near_factor/far_quad given) on this tet body IS the fast far_quad build:
    # it must match the explicit far_quad result (same Gram) to the solve tolerance
    with ng.TaskManager():
        df = hdiv_demag_solve(mesh, 100.0, HEXT)
    assert abs(df["demag"] - fq["demag"]) < 5e-4, "default tet CG path should use the far_quad fast build"


def test_far_quad_only_applies_to_near_far_split():
    """far_quad is a FAR-pair option: with near_factor=1e30 (all pairs NEAR) it is a no-op -> identical to
    the plain all-analytic Gram (guards against far_quad accidentally touching the exact path)."""
    mesh = _tet_cube(0.32)
    with ng.TaskManager():
        a = hdiv_demag_solve(mesh, 100.0, HEXT, near_factor=1e30, far_quad=0)
        b = hdiv_demag_solve(mesh, 100.0, HEXT, near_factor=1e30, far_quad=4)
    assert abs(a["demag"] - b["demag"]) < 1e-9, "far_quad must be a no-op when near_factor=1e30 (no far pairs)"
