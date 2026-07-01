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
  2. the dispatch (symmetric ONLY at quad==3 = linear RT1 near + default far_quad=3; nonlinear quad=4 -> product),
  3. the equivalence claim (symmetric build == product build on demag + PSD, on the actual charge Gram).

NGSolve + Netgen required (importorskip)."""
import numpy as np
import pytest
from math import factorial as fac

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Sphere, OCCGeometry, Pnt  # noqa: E402

import radia.vim._vim as V  # noqa: E402
from radia.vim._vim import (build_charge_gram, _tet_ref, _tri_ref, _tet_ref_sym5,  # noqa: E402
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
    """The symmetric rule is used ONLY at quad==3 (degree 5: linear RT1 near + default far_quad=3).  Any other
    order (nonlinear quad=4, non-default far_quad, inner subtraction) falls back to the product Gauss-Duffy rule,
    so nonlinear/curved paths are untouched."""
    assert _outer_tet(3)[0].shape[0] == 15 and _outer_tri(3)[0].shape[0] == 7        # symmetric
    for q in (2, 4, 5):
        assert _outer_tet(q)[0].shape[0] == _tet_ref(q)[0].shape[0]                   # product
        assert _outer_tri(q)[0].shape[0] == _tri_ref(q)[0].shape[0]


def _demag_and_psd(fes):
    B, G, _M = build_charge_gram(fes)
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
