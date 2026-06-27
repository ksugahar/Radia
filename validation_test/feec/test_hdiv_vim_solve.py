"""Golden test: end-to-end SCALABLE SOLVE of the HDiv-type VIM material system using the HACApK
H-matrix operator (rad_hacapk_hdiv, _HDivVimHMatrix).  This is Layer A.5 -- turning the H-matrix
operator into an actual solver.

The material system is A x = b with A = inv_chi * M_mass - N (symmetric INDEFINITE; its generalized
eigenvalues vs M_mass are inv_chi minus the demag factors, which lie in [0,1]).  We solve it with MINRES
(the principled symmetric-indefinite Krylov method) driven entirely by the build-once H-matrix
operator's apply_system (= inv_chi M_mass x - N x; N x is the O(N log N) H-matvec) and Jacobi
preconditioner (diag_system).  Nothing dense enters the solve -- the H-matvec is the only cost.

Four locks:
  (1) apply_system equals the dense A x within the ACA tolerance (operator correctness);
  (2) when the grid is all-dense (exact N) the MINRES solution equals the dense direct solve to the
      MINRES tolerance (SOLVER correctness, decoupled from ACA);
  (3) on larger grids the MINRES solution matches the dense direct solve within ~ACA tol, on
      regular AND distorted grids, nsub=0 and nsub>=1;
  (4) the preconditioned-MINRES iteration count stays BOUNDED as mu_r ranges 10 -> 1e5
      (mu_r-INDEPENDENCE -- the property that makes the loop-field-null HDiv operator strong).
"""
import inspect

import numpy as np
import pytest

import radia._radia_pybind as _rp

minres = pytest.importorskip("scipy.sparse.linalg").minres
from scipy.sparse.linalg import LinearOperator  # noqa: E402

# scipy renamed the MINRES tolerance kwarg tol -> rtol around 1.12; support both.
_TOLKW = "rtol" if "rtol" in inspect.signature(minres).parameters else "tol"


def _dense_A(nx, ny, nz, nsub, distort, inv_chi):
    d = _rp._hdiv_vim_assemble(nx, ny, nz, nsub, distort)
    nf = d["nf"]
    N = np.asarray(d["N"], float).reshape(nf, nf)
    M = np.asarray(d["M_mass"], float).reshape(nf, nf)
    return inv_chi * M - N, nf


def _op(H, inv_chi):
    nf = H.ndof()
    A = LinearOperator((nf, nf),
                       matvec=lambda v: np.asarray(H.apply_system(np.asarray(v).tolist(), inv_chi), float))
    dvec = np.maximum(np.abs(np.asarray(H.diag_system(inv_chi), float)), 1e-30)
    Minv = LinearOperator((nf, nf), matvec=lambda v: np.asarray(v, float) / dvec)
    return A, Minv


class _Counter:
    def __init__(self):
        self.n = 0

    def __call__(self, *a):
        self.n += 1


def test_apply_system_matches_dense():
    """apply_system(x, inv_chi) reproduces the dense A x = (inv_chi M_mass - N) x within the ACA tol."""
    nx = ny = nz = 4
    inv_chi = 1e-3
    A_dense, nf = _dense_A(nx, ny, nz, 0, 0.0, inv_chi)
    H = _rp._HDivVimHMatrix(nx, ny, nz, 0, 0.0, 1e-6, 32, 2.0)
    rng = np.random.default_rng(0)
    x = rng.standard_normal(nf)
    y_h = np.asarray(H.apply_system(x.tolist(), inv_chi), float)
    y_d = A_dense @ x
    rel = np.linalg.norm(y_h - y_d) / np.linalg.norm(y_d)
    assert rel < 1e-4, f"apply_system vs dense A: rel err {rel:.2e}"


def test_minres_exact_when_all_dense():
    """2x2x2 leaf=64 -> exact (all-dense) N -> the MINRES solution equals the dense direct solve to
    the MINRES tolerance.  Locks the SOLVER independent of any ACA approximation."""
    inv_chi = 1e-3
    A_dense, nf = _dense_A(2, 2, 2, 0, 0.0, inv_chi)
    H = _rp._HDivVimHMatrix(2, 2, 2, 0, 0.0, 1e-4, 64, 2.0)
    assert H.stats()["n_lowrank"] == 0
    rng = np.random.default_rng(1)
    b = rng.standard_normal(nf)
    x_dense = np.linalg.solve(A_dense, b)
    A, Minv = _op(H, inv_chi)
    x_h, info = minres(A, b, M=Minv, maxiter=4000, **{_TOLKW: 1e-10})
    rel = np.linalg.norm(x_h - x_dense) / np.linalg.norm(x_dense)
    assert rel < 1e-7, f"all-dense MINRES vs direct: rel err {rel:.2e} (info={info})"


@pytest.mark.parametrize("nx,ny,nz,nsub,distort", [
    (3, 3, 3, 0, 0.0),
    (4, 4, 4, 0, 0.0),
    (4, 4, 4, 2, 0.0),
    (4, 4, 4, 0, 0.18),
])
def test_minres_matches_dense_within_aca_tol(nx, ny, nz, nsub, distort):
    """The H-matrix MINRES solution matches the dense direct solve within ~ACA tol (the H-matrix
    approximates N to eps; that error propagates through A^{-1})."""
    eps = 1e-6
    inv_chi = 1e-3
    A_dense, nf = _dense_A(nx, ny, nz, nsub, distort, inv_chi)
    H = _rp._HDivVimHMatrix(nx, ny, nz, nsub, distort, eps, 32, 2.0)
    rng = np.random.default_rng(2)
    b = rng.standard_normal(nf)
    x_dense = np.linalg.solve(A_dense, b)
    A, Minv = _op(H, inv_chi)
    x_h, info = minres(A, b, M=Minv, maxiter=5000, **{_TOLKW: 1e-9})
    rel = np.linalg.norm(x_h - x_dense) / np.linalg.norm(x_dense)
    assert rel < 1e-3, f"MINRES vs dense rel err {rel:.2e} (nx={nx} nsub={nsub} distort={distort} info={info})"


def test_minres_iters_bounded_vs_mu_r():
    """Preconditioned-MINRES iteration count stays BOUNDED as mu_r ranges 10 -> 1e5: the HDiv
    operator has the loops field-null by construction, so high mu_r does not create a growing
    near-null space that wrecks convergence (the property that motivates the HDiv formulation)."""
    nx = ny = nz = 4
    H = _rp._HDivVimHMatrix(nx, ny, nz, 0, 0.0, 1e-6, 32, 2.0)
    nf = H.ndof()
    rng = np.random.default_rng(3)
    b = rng.standard_normal(nf)
    iters = {}
    for mu_r in (1e1, 1e3, 1e5):
        inv_chi = 1.0 / (mu_r - 1.0)
        A, Minv = _op(H, inv_chi)
        c = _Counter()
        _x, _info = minres(A, b, M=Minv, callback=c, maxiter=5000, **{_TOLKW: 1e-8})
        iters[mu_r] = c.n
    print(f"  MINRES iters vs mu_r: {iters}")
    lo, hi = min(iters.values()), max(iters.values())
    assert hi < 200, f"MINRES iteration count too high: {iters}"
    assert hi <= 4 * lo + 30, f"iteration count not mu_r-bounded (blows up with mu_r): {iters}"


def test_minres_iters_bounded_vs_mu_r_distorted():
    """The mu_r-INDEPENDENCE SURVIVES MESH DISTORTION -- the practical payoff of the HDiv formulation.

    On a DISTORTED grid (distort=0.18, accurate trilinear Gram) the preconditioned-MINRES iteration
    count stays BOUNDED as mu_r ranges 10 -> 1e4 (all real soft iron).  This works because the loops are
    field-null BY CONSTRUCTION on ANY mesh (N = B^T G B, so any ker(B) loop gives N.loop = 0 exactly,
    affine OR distorted) -- distortion does NOT create a growing near-null space.  This is the de-Rham
    HDiv advantage over the six-face surface-charge combinatorial +/-1 loops, which are field-null ONLY on affine
    hexes and carry field on distorted ones (the de-Rham defect the hand-crafted surface-charge elements /
    the shipped MSC's installCycle retrofit have to work around).  Here the iters even DECREASE with mu_r
    (the loops, which would otherwise blow up, sit exactly at field-null)."""
    nx = ny = nz = 4
    H = _rp._HDivVimHMatrix(nx, ny, nz, 4, 0.18, 1e-6, 32, 2.0)   # nsub=4 accurate Gram, distort=0.18
    nf = H.ndof()
    rng = np.random.default_rng(3)
    b = rng.standard_normal(nf)
    iters = {}
    for mu_r in (1e1, 1e2, 1e3, 1e4):
        inv_chi = 1.0 / (mu_r - 1.0)
        A, Minv = _op(H, inv_chi)
        c = _Counter()
        minres(A, b, M=Minv, callback=c, maxiter=5000, **{_TOLKW: 1e-8})
        iters[mu_r] = c.n
    print(f"  distorted MINRES iters vs mu_r: {iters}")
    assert max(iters.values()) < 200, f"distorted MINRES iters too high: {iters}"
    # high mu_r must not blow up the count relative to low mu_r (the mu_r-independence statement)
    assert iters[1e4] <= 2.0 * iters[1e1], f"distorted iters blow up with mu_r: {iters}"
