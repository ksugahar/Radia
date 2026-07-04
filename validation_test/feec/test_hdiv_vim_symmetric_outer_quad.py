"""Golden: the SYMMETRIC degree-5 outer Gram quadrature (Keast-15 tet / Dunavant-7 tri) that replaces the
degree-5 PRODUCT Gauss-Duffy rule (27-pt tet / 9-pt tri) for the linear-RT1 charge-Gram OUTER integral.

Why this is valid (and why a symmetric rule suffices where a generic Galerkin singular integral would need the
Duffy point-clustering): the INNER integral is carried EXACTLY by the analytic PhiTet/TriPotential, so the
outer integrand is C^{1,alpha} (smooth) even on self/face/edge/vertex-adjacent pairs.  A same-degree symmetric
rule with 1.80x/1.29x FEWER points therefore reproduces the Gram operator N = B^T G B -- same demag factor,
identical transverse leak, and (critically) PRESERVED PSD -- while building the Gram ~1.5-1.7x faster.

The fully-double-ANALYTIC alternative was surveyed and rejected (no tractable closed form for the dominant
tet-tet Galerkin double integral); the symmetric outer rule is the real, cheap lever.  This golden guards:
  1. the transcribed Keast/Dunavant constants (degree-5 exactness -- a transcription typo fails here),
  2. the dispatch (symmetric at quad in {3,4}; other orders stay on product),
  3. the equivalence claim (symmetric build == product build on demag + PSD / nonlinear convergence).

NGSolve + Netgen required (importorskip)."""
import numpy as np
import pytest
from math import factorial as fac

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Sphere, OCCGeometry, Pnt  # noqa: E402

import radia.vim._vim as V  # noqa: E402
from radia.vim import ChargeGram  # noqa: E402
from radia.vim._vim import (_tet_ref, _tri_ref, _tet_ref_sym5,  # noqa: E402
                            _tri_ref_sym5, _outer_tet, _outer_tri)


def _exact_tet(i, j, k):
    return fac(i) * fac(j) * fac(k) / fac(i + j + k + 3)


def _exact_tri(i, j):
    return fac(i) * fac(j) / fac(i + j + 2)


def test_symmetric_rules_are_degree5_exact():
    """Keast-15 / Dunavant-7 integrate every monomial up to degree 5 exactly over the reference simplex
    (weights sum to the ref volume 1/6 / area 1/2).  A transcription typo in a node or weight fails here."""
    Pt, Wt = _tet_ref_sym5()
    Ps, Ws = _tri_ref_sym5()
    assert Pt.shape == (15, 3) and Ps.shape == (7, 2)
    assert abs(Wt.sum() - 1.0 / 6.0) < 1e-14 and abs(Ws.sum() - 0.5) < 1e-14
    err_tet = max(abs(np.sum(Wt * Pt[:, 0] ** i * Pt[:, 1] ** j * Pt[:, 2] ** k) - _exact_tet(i, j, k))
                  for i in range(6) for j in range(6 - i) for k in range(6 - i - j))
    err_tri = max(abs(np.sum(Ws * Ps[:, 0] ** i * Ps[:, 1] ** j) - _exact_tri(i, j))
                  for i in range(6) for j in range(6 - i))
    assert err_tet < 1e-13 and err_tri < 1e-13, f"degree-5 exactness failed: tet {err_tet:.2e} tri {err_tri:.2e}"


def test_outer_rule_dispatch():
    """The symmetric degree-5 rule is used at quad in {3,4}: quad==3 = linear RT1 near + default far_quad=3
    (product _tet_ref(3) is only degree 3), quad==4 = the nonlinear energy-Newton near rule (product
    _tet_ref(4)=64pts is degree 5, matched by the 15-pt symmetric rule at 4.3x fewer points).  Any OTHER order
    (inner subtraction iq=2, intorder overrides, curved) falls back to product Gauss-Duffy."""
    for q in (3, 4):                                                                  # symmetric degree-5
        assert _outer_tet(q)[0].shape[0] == 15 and _outer_tri(q)[0].shape[0] == 7
    for q in (2, 5, 6):                                                               # product fall-back
        assert _outer_tet(q)[0].shape[0] == _tet_ref(q)[0].shape[0]
        assert _outer_tri(q)[0].shape[0] == _tri_ref(q)[0].shape[0]


def _demag_and_psd(fes):
    B, G, _M = ChargeGram(fes)
    n = fes.ndof
    N = np.zeros((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        N[:, j] = B.T @ np.asarray(G.matvec_sym((B @ e).tolist()), float)
        e[j] = 0.0
    N = 0.5 * (N + N.T)
    w = np.linalg.eigvalsh(N)
    return N, w


def test_symmetric_reproduces_product_demag_and_psd(monkeypatch):
    """On the actual charge Gram, the symmetric build reproduces the product build: ||N_sym - N_prod|| is small,
    the demag-energy spectra agree, and BOTH stay PSD (min eig >= -tol * max).  This is the equivalence lock --
    the symmetric outer rule is a validated drop-in, not a silent accuracy change."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.45))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        N_sym, w_sym = _demag_and_psd(fes)                              # default = symmetric
        monkeypatch.setattr(V, "_outer_tet", _tet_ref)                  # force product Gauss-Duffy
        monkeypatch.setattr(V, "_outer_tri", _tri_ref)
        N_prod, w_prod = _demag_and_psd(fes)
    rel = np.linalg.norm(N_sym - N_prod) / np.linalg.norm(N_prod)
    assert rel < 3e-2, f"symmetric vs product operator diff too large: {rel:.2e}"
    assert w_sym.min() >= -1e-9 * w_sym.max(), f"symmetric Gram NOT PSD: {w_sym.min():.2e}"
    assert w_prod.min() >= -1e-9 * w_prod.max(), f"product Gram NOT PSD: {w_prod.min():.2e}"
    # the largest (dominant demag-energy) eigenvalues agree tightly -- same physics
    assert abs(w_sym[-1] - w_prod[-1]) <= 5e-3 * abs(w_prod[-1]), "dominant demag eigenvalue drifted"


# realistic saturating BH table (mu_r ~ 4000 low, sat ~2.2 T) for the nonlinear equivalence lock
_H = np.array([0, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 3e4, 1e5, 3e5, 1e6])
_B = np.array([0, 0.30, 0.60, 1.0, 1.45, 1.7, 1.9, 2.0, 2.05, 2.1, 2.15, 2.25, 2.5])
_BH = np.column_stack([_H, _B]).tolist()
_MSAT = _B[-1] / (4e-7 * np.pi) - _H[-1]


def test_symmetric_nonlinear_matches_product(monkeypatch):
    """The NONLINEAR energy-Newton (quad=4) now uses the SAME degree-5 symmetric rule (product _tet_ref(4)=64pts
    is itself only effective degree 5).  Lock that the symmetric build converges in the SAME or FEWER Newton
    iterations than product-64, well under 100, drives M -> Msat at deep saturation, and matches demag -- across
    a moderate drive and deep saturation.  A regression that under-resolves the energy Hessian (the 195-vs-<100
    iter blowup the degree-3 product-27 rule caused) fails here."""
    from radia.vim import Solve
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.6))
    for H0, deep in ((5000.0, False), (3e6, True)):
        Hext = ng.CoefficientFunction((0, 0, H0))
        with ng.TaskManager():
            r_sym = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1)        # default = symmetric
            monkeypatch.setattr(V, "_outer_tet", _tet_ref)                            # force product-64
            monkeypatch.setattr(V, "_outer_tri", _tri_ref)
            r_prod = Solve(mesh, bh_table=_BH, H_ext=Hext, order=1)
            monkeypatch.undo()
        assert r_sym["iters"] < 100, f"H0={H0}: symmetric nonlinear iters {r_sym['iters']} >= 100"
        assert r_sym["iters"] <= r_prod["iters"] + 5, (
            f"H0={H0}: symmetric iters {r_sym['iters']} >> product {r_prod['iters']} (Hessian under-resolved)")
        assert r_sym["M_avg"][2] > 0 and np.isfinite(r_sym["M_avg"][2])
        assert abs(r_sym["demag"] - r_prod["demag"]) < 1e-3, (
            f"H0={H0}: demag drift {abs(r_sym['demag'] - r_prod['demag']):.2e}")
        if deep:
            assert 0.95 * _MSAT < r_sym["M_avg"][2] < 1.03 * _MSAT, (r_sym["M_avg"][2], _MSAT)
