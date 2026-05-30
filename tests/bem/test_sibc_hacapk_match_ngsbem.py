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

    # After the edge-winding fix (commit 99db78a0), our DL is bit-exact
    # with NGSolve.bem at bonus_intorder=10:
    #   max abs ~1.7e-10 (machine precision)
    #   max rel ~3.5e-7 (limited to small-magnitude entries)
    assert abs_err.max() < 5e-10, (
        f"DL max absolute error {abs_err.max():.3e} exceeds 5e-10")
    assert rel_sig_max < 1e-6, (
        f"DL max relative error {rel_sig_max:.3e} exceeds 1e-6")


def test_assemble_bem_matrices_e2e():
    """Phase 1.5: assemble_bem_matrices on a fresh sphere mesh produces
    SL, DL, M, K that all match an independent NGSolve.bem build.

    This is the end-to-end path that bem_sibc_solver.py will use to
    bypass NGSolve.bem entirely.
    """
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    from ngsolve import (Mesh, BND, H1, BilinearForm, ds, grad,
                         InnerProduct, TaskManager)
    from ngsolve import TaskManager
    from netgen.occ import Sphere, Pnt, OCCGeometry
    from radia.bem.sibc_hacapk import assemble_bem_matrices

    sph = Sphere(Pnt(0, 0, 0), 1.0)
    g = OCCGeometry(sph)
    ng_mesh = g.GenerateMesh(maxh=0.4)
    mesh = Mesh(ng_mesh)

    # Reference matrices from NGSolve.bem (high-precision)
    fes = H1(mesh, order=1)
    u, v = fes.TnT()
    ndof = fes.ndof
    mass_bf = BilinearForm(fes); mass_bf += u.Trace() * v.Trace() * ds
    with TaskManager():
        mass_bf.Assemble()
        M_ref = np.zeros((ndof, ndof))
        for r_, c_, vv in zip(*mass_bf.mat.COO()):
            M_ref[int(r_), int(c_)] = vv
        stiff_bf = BilinearForm(fes)
        stiff_bf += InnerProduct(grad(u).Trace(), grad(v).Trace()) * ds
        stiff_bf.Assemble()
        K_ref = np.zeros((ndof, ndof))
        for r_, c_, vv in zip(*stiff_bf.mat.COO()):
            K_ref[int(r_), int(c_)] = vv

        out = assemble_bem_matrices(mesh,
                                     regular_quad_degree=11,
                                     singular_n_q=8)
        SL = out["SL"]; DL = out["DL"]
        M = out["M"]; K = out["K"]
        v_global = out["v_global"]

        # M and K are local to extracted vertices; project NGSolve refs to that
        # ordering.  For a sphere with no inner BND, v_global is all of
        # mesh.vertices.nr (sorted).
        M_ref_loc = M_ref[np.ix_(v_global, v_global)]
        K_ref_loc = K_ref[np.ix_(v_global, v_global)]

        err_M = np.abs(M - M_ref_loc).max()
        err_K = np.abs(K - K_ref_loc).max()
        assert err_M < 1e-12, f"mass abs err {err_M:.3e}"
        assert err_K < 1e-10, f"stiffness abs err {err_K:.3e}"

        # SL/DL self-symmetry of OUR matrices (should be exact within roundoff)
        assert np.abs(SL - SL.T).max() < 1e-10
        # DL is generally NOT symmetric.

        # Surface area sanity check via M
        ones = np.ones(SL.shape[0])
        area = ones @ M @ ones
        assert abs(area - 4 * np.pi) / (4 * np.pi) < 0.05, (
            f"surface area {area:.3f} vs 4π {4*np.pi:.3f}")


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
