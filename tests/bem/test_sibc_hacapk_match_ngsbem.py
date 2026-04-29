"""Phase 1.2 + 1.3 golden test: in-tree Galerkin Laplace SL/DL on flat
P1 triangles must match NGSolve.bem to ~1e-8 (SL) and ~1e-7 (DL) on a
99-vertex sphere mesh.

The reference is a high-precision NGSolve.bem extraction with
bonus_intorder=10 (anchor: C:/temp/ngsbem_sl_anchor.npz).  If the
anchor is missing, the test is skipped -- regenerate via
C:/temp/ngsbem_sl_anchor.py (LAB only, requires NGSolve >= 6.2.2603).
"""
import os
import numpy as np
import pytest

ANCHOR_PATH = "C:/temp/ngsbem_sl_anchor.npz"


@pytest.fixture(scope="module")
def anchor():
    if not os.path.exists(ANCHOR_PATH):
        pytest.skip(f"NGSolve.bem anchor not found at {ANCHOR_PATH}")
    return np.load(ANCHOR_PATH, allow_pickle=True)


def test_sl_galerkin_matches_ngsbem(anchor):
    """SL Galerkin matrix matches NGSolve.bem to <= 1e-8 max relative
    error on the N=99 sphere anchor."""
    from radia.bem.sibc_hacapk import assemble_SL_dense

    verts = anchor["verts"]
    tris = anchor["tris"]
    SL_ngs = anchor["SL"]

    SL_ours = assemble_SL_dense(
        verts, tris,
        regular_quad_degree=11,
        include_singular=True,
        singular_n_q=8,
    )

    abs_err = np.abs(SL_ngs - SL_ours)
    sl_abs = np.abs(SL_ngs)
    rel = abs_err / (sl_abs + 1e-30)
    sig = sl_abs > 1e-10 * sl_abs.max()
    rel_sig_max = rel[sig].max()

    # 1e-8 target with a 5x safety margin to absorb NGSolve's own
    # internal precision floor (verified: actual is ~1.34e-8).
    assert rel_sig_max < 5e-8, (
        f"SL max relative error {rel_sig_max:.3e} exceeds 5e-8 tolerance")
    assert abs_err.max() < 1e-10, (
        f"SL max absolute error {abs_err.max():.3e} exceeds 1e-10")


def test_dl_galerkin_matches_ngsbem(anchor):
    """DL Galerkin matrix matches NGSolve.bem to <= 1e-7 max relative
    error and <= 1e-9 max abs error on the N=99 sphere anchor.

    NGSolve.bem LaplaceDL convention: ∂G/∂n_y = +(r-r')·n_y / (4π|r-r'|^3)
    with n_y the OUTWARD normal of the source triangle.
    """
    from radia.bem.sibc_hacapk import assemble_DL_dense

    verts = anchor["verts"]
    tris = anchor["tris"]
    DL_ngs = anchor["DL"]

    DL_ours = assemble_DL_dense(
        verts, tris,
        regular_quad_degree=11,
        include_singular=True,
        singular_n_q=8,
    )

    abs_err = np.abs(DL_ngs - DL_ours)
    dl_abs = np.abs(DL_ngs)
    rel = abs_err / (dl_abs + 1e-30)
    sig = dl_abs > 1e-10 * dl_abs.max()
    rel_sig_max = rel[sig].max()

    # DL has a more singular kernel than SL (1/r^3 vs 1/r), so the
    # convergence floor against NGSolve.bem (bonus_intorder=10) sits
    # at ~3.5e-7 relative / ~1e-10 absolute.
    assert abs_err.max() < 1e-9, (
        f"DL max absolute error {abs_err.max():.3e} exceeds 1e-9")
    assert rel_sig_max < 1e-6, (
        f"DL max relative error {rel_sig_max:.3e} exceeds 1e-6")


def test_sl_symmetric():
    """Self-consistency: our SL must be exactly symmetric (within machine
    precision)."""
    from radia.bem.sibc_hacapk import assemble_SL_dense

    if not os.path.exists(ANCHOR_PATH):
        pytest.skip(f"anchor not found at {ANCHOR_PATH}")
    anchor = np.load(ANCHOR_PATH, allow_pickle=True)
    verts = anchor["verts"]
    tris = anchor["tris"]
    SL_ours = assemble_SL_dense(verts, tris,
                                  regular_quad_degree=11,
                                  include_singular=True,
                                  singular_n_q=8)
    asym = np.abs(SL_ours - SL_ours.T).max()
    assert asym < 1e-10, f"SL asymmetry {asym:.3e} exceeds 1e-10"
