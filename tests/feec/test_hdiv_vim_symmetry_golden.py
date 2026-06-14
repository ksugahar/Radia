"""Golden test for the symmetric HDiv-type VIM demag operator (Radia C++ core, rad_hdiv_vim).

Locks the structural foundation of the HDiv-type VIM (the symmetric alternative to the collocation
MSC kernel): on a structured nx*ny*nz hex grid (RT0 faces), the demag operator N = B^T G B is
  (1) SYMMETRIC (Galerkin energy form) -- ||N - N^T||/||N|| ~ machine eps, and
  (2) loops are FIELD-NULL BY CONSTRUCTION (loops = ker B; B.loop = 0 => N.loop = 0).
Golden values come from the NGSolve prototype (examples/vim/hdiv_demag_quad_self.json):
regular 3x3x3 -> ndof=108, n_loop=28, asym~1e-16, loop_res~1e-16.  The C++ hand-enumerated topology
(rad_hdiv_vim, no NGSolve) reproduces these exactly (verified standalone before integration).
"""
import numpy as np
import pytest

import radia._radia_pybind as _rp


def _assemble(nx, ny, nz):
    d = _rp._hdiv_vim_assemble(nx, ny, nz)
    nf, n_charge = d["nf"], d["n_charge"]
    N = np.asarray(d["N"], float).reshape(nf, nf)
    B = np.asarray(d["B"], float).reshape(n_charge, nf)
    return d, N, B


@pytest.mark.parametrize("nx,ny,nz,golden_ndof,golden_nloop", [
    (1, 1, 1, 6, 0),
    (2, 2, 2, 36, 5),
    (3, 3, 3, 108, 28),   # the NGSolve prototype golden (hdiv_demag_quad_self.json)
    (4, 4, 4, 240, 81),   # scaling table (hdiv_demag_scaling.json)
    (5, 5, 5, 450, 176),
])
def test_hdiv_vim_symmetry_and_loop_nullity(nx, ny, nz, golden_ndof, golden_nloop):
    d, N, B = _assemble(nx, ny, nz)
    nf = d["nf"]
    assert nf == golden_ndof, f"ndof {nf} != golden {golden_ndof}"

    # (1) symmetry (Galerkin energy form)
    nN = np.linalg.norm(N, 2)
    asym = np.linalg.norm(N - N.T) / nN
    assert asym < 1e-12, f"N not symmetric: ||N-N^T||/||N|| = {asym:.2e}"

    # (2) loops = ker(B); count + field-nullity
    s = np.linalg.svd(B, compute_uv=False)
    rankB = int(np.sum(s > 1e-9 * s.max()))
    n_loop = nf - rankB
    assert n_loop == golden_nloop, f"n_loop {n_loop} != golden {golden_nloop}"

    if n_loop:
        _, _, Vt = np.linalg.svd(B)
        loops = Vt[rankB:, :]                      # ker(B) basis
        loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / nN
        assert loop_res < 1e-10, f"loops not field-null: max||N.loop||/||N|| = {loop_res:.2e}"


def test_hdiv_vim_charge_conservation_gauss():
    """rank(B) = n_charge - 1 (Gauss: the total-charge mode is the single redundancy)."""
    d, N, B = _assemble(3, 3, 3)
    s = np.linalg.svd(B, compute_uv=False)
    rankB = int(np.sum(s > 1e-9 * s.max()))
    assert rankB == d["n_charge"] - 1, f"rank(B)={rankB} != n_charge-1={d['n_charge']-1}"


def test_hdiv_vim_demag_factors_physical():
    """PHYSICS lock (Phase 4 accurate Gram, nsub>=1): the generalized eigenvalues of (N, M_mass) are
    the demag factors and must be PHYSICAL -- positive semi-definite (>= -1e-6; a negative demag
    factor is unphysical) and all <= 1.  The accurate sub-point Gram makes N PSD and drives the
    uniform mode to the true cube 1/3 (uniform M_z demag -> 0.32 -> 1/3 with nsub, validated in
    examples/vim/hdiv_vim_accurate_g.py).  Also catches the surface-charge SIGN bug (M.n must
    use the OUTWARD normal; the global-normal bug gave demag factors up to 19.5) -- which the
    symmetry/loop-nullity structural tests do NOT catch (a row sign-flip in B preserves ker B +
    N's symmetry)."""
    import scipy.linalg as sla
    d = _rp._hdiv_vim_assemble(3, 3, 3, 4)         # nsub=4 -> accurate, positive-semidefinite Gram
    nf = d["nf"]
    N = np.asarray(d["N"], float).reshape(nf, nf)
    M = np.asarray(d["M_mass"], float).reshape(nf, nf)
    assert np.allclose(M, M.T) and np.all(np.linalg.eigvalsh(M) > 0), "M_mass not SPD"
    asym = np.linalg.norm(N - N.T) / np.linalg.norm(N, 2)
    assert asym < 1e-12, f"accurate N not symmetric: {asym:.2e}"
    df = np.sort(sla.eigh(N, M, eigvals_only=True))
    assert df[0] > -1e-6, f"accurate Gram not PSD: min demag factor {df[0]:.3e} (should be >= 0)"
    assert df[-1] < 1.0 + 1e-6, f"unphysical demag factor > 1: max={df[-1]:.4f} (surface-charge sign bug?)"
    assert 0.80 < df[-1] < 0.95, f"top demag factor {df[-1]:.4f} not ~0.874"


def test_hdiv_vim_distorted_trilinear_psd():
    """Phase 3: on a DISTORTED grid the accurate Gram uses TRILINEAR (hex) / BILINEAR (quad)
    sub-points that follow the actual cell geometry, keeping N positive semi-definite (demag factors
    in [0,1]).  The regular cube-grid sub-points do NOT follow the distorted cell -> unphysical demag
    factors (negative AND > 1).  Symmetry + loops-field-null still hold on the distorted grid."""
    import scipy.linalg as sla
    nf_ref = None
    d = _rp._hdiv_vim_assemble(3, 3, 3, 4, 0.18)        # nsub=4 accurate (trilinear), distort=0.18
    nf = d["nf"]
    N = np.asarray(d["N"], float).reshape(nf, nf)
    M = np.asarray(d["M_mass"], float).reshape(nf, nf)
    B = np.asarray(d["B"], float).reshape(d["n_charge"], nf)
    assert np.linalg.norm(N - N.T) / np.linalg.norm(N, 2) < 1e-12, "not symmetric on distorted"
    s = np.linalg.svd(B, compute_uv=False); rankB = int(np.sum(s > 1e-9 * s.max()))
    Vt = np.linalg.svd(B)[2]; loops = Vt[rankB:]
    lres = max(np.linalg.norm(N @ loops[k]) for k in range(nf - rankB)) / np.linalg.norm(N, 2)
    assert lres < 1e-10, f"loops not field-null on distorted: {lres:.1e}"
    df = np.sort(sla.eigh(N, M, eigvals_only=True))
    assert df[0] > -1e-6, f"accurate (trilinear) Gram not PSD on distorted: min demag {df[0]:.2e}"
    assert df[-1] < 1.0 + 1e-6, f"unphysical demag factor > 1 on distorted: {df[-1]:.4f}"
    # contrast: the cube-grid sub-points (nsub=0 centroid-monopole) are non-PD on the distorted grid
    d0 = _rp._hdiv_vim_assemble(3, 3, 3, 0, 0.18)
    N0 = np.asarray(d0["N"], float).reshape(nf, nf)
    df0 = np.sort(sla.eigh(N0, M, eigvals_only=True))
    assert df0[0] < -1e-3, "expected the cube-grid Gram to be non-PD on the distorted grid (contrast)"


def test_hdiv_vim_centroid_gram_fast_baseline():
    """The fast centroid-monopole Gram (nsub=0, default) stays symmetric + loop-null (G-agnostic
    properties) and is ~7% off / slightly non-PD -- a documented crudeness used for the structural
    tests; the accurate Gram (nsub>=1) is the physical one (test above)."""
    d, N, _ = _assemble(3, 3, 3)                    # nsub=0 default
    nf = d["nf"]
    assert np.linalg.norm(N - N.T) / np.linalg.norm(N, 2) < 1e-12
