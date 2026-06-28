"""Golden: the HIGH-ORDER near/far adaptive-quadrature split is ACCURACY-PRESERVING and actually fires.

The order-p charge-Gram has an OPT-IN build speedup: well-separated charge pairs use a cheap LOW-quad plain
double-Gauss (QuadDotFar) instead of the full high-quad singularity-subtraction (QuadDot).  This is sound
ONLY when (a) "well-separated" is judged by each host's BOUNDING RADIUS (so touching high-aspect-ratio
needle/sliver pairs are never misclassified FAR -- the bug fixed in rad_hacapk_hdiv.cpp: isotropic
cbrt(vol)/sqrt(area) underestimated the extent), and (b) far_quad resolves the smooth far kernel.

This test builds the SAME Gram two ways via build_charge_gram and asserts they agree:
  * ho_far_factor=inf  -> the split is DISABLED, every pair uses the exact high-quad path (the default,
    golden-locked, provably-exact Gram);
  * ho_far_factor=2.0, far_quad=3 -> the split is ON.

On an ELONGATED box the end charges are genuinely well-separated, so QuadDotFar is guaranteed to fire
(reldiff > 0 proves it ran; reldiff < 1e-3 proves it stayed accurate).  On an anisotropic NEEDLE-CHAIN bar
(long thin elements stacked end-to-end, the case the old isotropic cbrt(vol) size actually misclassified)
the bounding-radius FAR test must keep the touching pairs NEAR, so split == exact there too -- the
regression guard for the TEAM-13 thin-steel misclassification.

NGSolve + Netgen required (importorskip); the C++ charge-Gram H-matrix backs the build.  The CALLER wraps
in TaskManager per the lab policy.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

from radia.vim import build_charge_gram, DemagOperator  # noqa: E402


def _box(lx, ly, lz, h):
    return ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(lx, ly, lz))).GenerateMesh(maxh=h))


def _gram_matvec_reldiff(fes, c):
    """||G_split c - G_exact c|| / ||G_exact c||; G_split = the opt-in far split, G_exact = all high-quad.
    B is independent of ho_far_factor, so comparing G alone isolates the QuadDotFar entries."""
    _, G_exact, _ = build_charge_gram(fes, ho_far_factor=float("inf"))
    _, G_split, _ = build_charge_gram(fes, ho_far_factor=2.0, far_quad=3)
    y_exact = np.asarray(G_exact.matvec(c.tolist()))
    y_split = np.asarray(G_split.matvec(c.tolist()))
    return np.linalg.norm(y_split - y_exact) / np.linalg.norm(y_exact)


@pytest.mark.parametrize("p", [1])
def test_far_split_fires_and_matches_exact_elongated(p):
    """ELONGATED box -> guaranteed FAR pairs -> QuadDotFar fires (reldiff > 0) AND is accuracy-preserving
    (reldiff < 1e-3) vs the exact all-high-quad Gram."""
    mesh = _box(6.0, 1.0, 1.0, 1.0)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=p)
        n_charge_probe, G_ex, _ = build_charge_gram(fes, ho_far_factor=float("inf"))
        n_charge = n_charge_probe.shape[0]
        rng = np.random.default_rng(0)
        c = rng.standard_normal(n_charge)
        reldiff = _gram_matvec_reldiff(fes, c)
    assert reldiff > 0.0, (
        f"p={p}: far split produced an IDENTICAL Gram -> QuadDotFar never fired (no FAR pairs on this "
        f"mesh); the test is vacuous, use a more elongated / finer mesh")
    assert reldiff < 1e-3, f"p={p}: far split deviates from exact by {reldiff:.2e} (>1e-3) -- not accuracy-preserving"


@pytest.mark.parametrize("p", [1])
def test_far_split_matches_exact_anisotropic_bar(p):
    """ANISOTROPIC NEEDLE CHAIN (aspect ~10): a long thin bar meshes into needle tets stacked end-to-end.
    Adjacent (TOUCHING) needles have centroid separation ~ their LENGTH, which far exceeds the isotropic
    cbrt(vol) -- so the OLD cbrt-size FAR test (r > ho_far_factor*(cbrt_a+cbrt_b)) misclassified these
    touching pairs as FAR and routed their near-SINGULAR integral to the subtraction-free QuadDotFar (the
    TEAM-13 thin-steel bug).  The bounding-radius size keeps them NEAR, so the far split must reproduce the
    exact build here; genuinely separated needles (opposite ends of the bar) still go FAR, exercising
    QuadDotFar.  (Unlike a 1-cell-thick slab, whose in-plane neighbours have centroid distance ~ their size
    and so are NOT misclassified even by the old size -- a needle chain is the geometry that actually fails
    on the buggy code, making this a real regression guard.)"""
    mesh = _box(4.0, 0.1, 0.1, 1.0)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=p)
        Bm, _, _ = build_charge_gram(fes, ho_far_factor=float("inf"))
        rng = np.random.default_rng(1)
        c = rng.standard_normal(Bm.shape[0])
        reldiff = _gram_matvec_reldiff(fes, c)
    assert reldiff < 1e-3, (
        f"p={p}: on the needle-chain bar the far split deviates by {reldiff:.2e} -- a touching needle pair "
        f"was misclassified FAR (isotropic-m_size regression: bounding-radius size not in effect)")


def test_demagfactor_split_matches_exact_cube():
    """End-to-end through the public DemagOperator: the (now-default) near/far split reproduces the exact
    ~1/3 demag.  The EXACT reference is pinned explicitly (ho_far_factor=inf) so this stays a real
    split-vs-exact check independent of the default."""
    mesh = _box(1.0, 1.0, 1.0, 0.4)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        D_exact = DemagOperator(fes, ho_far_factor=float("inf")).DemagFactor(ng.CF((0, 0, 1)))
        D_split = DemagOperator(fes, ho_far_factor=2.0, far_quad=3).DemagFactor(ng.CF((0, 0, 1)))
    assert 0.30 < D_exact < 0.36, f"exact DemagFactor {D_exact:.5f} not ~1/3"
    assert abs(D_split - D_exact) < 2e-3, f"split {D_split:.5f} vs exact {D_exact:.5f} disagree >2e-3"
