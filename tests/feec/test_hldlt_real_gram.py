"""Golden test: production #2 (integration) -- the rk-aware symmetric H-LDL^T factoring a REAL
HACApK-built Radia H-matrix end-to-end, via the cHACApK_hldlt_factor_leafmtxp driver exposed as
_ChargeGramHMatrix.factor_solve_hldlt.

The charge-charge Coulomb Gram G (the unstructured HDiv-VIM Gram, #1b) is symmetric, so it is the
natural target for the symmetric factorization.  Built on a DETERMINISTIC structured charge cloud
(no Netgen variability) so the HACApK cluster tree shape is reproducible.

Two locks (the HONEST state of #2):
  - the driver factors a real multi-leaf HACApK Gram H-matrix and solves G x = b to MACHINE precision
    (single dense leaf at leaf>=N, AND a genuine 7-leaf diagonal-refined tree at leaf=N//2);
  - on a DEEPER tree (small leaf -> the natural HACApK build creates INTERNAL off-diagonal blocks),
    the driver returns CHACAPK_HARITH_ERR_NEED_RECURSIVE (-5) -- it fails LOUD (no fallback), recording
    the boundary: the lower-only H-LDL^T does not yet recurse internal off-diagonal blocks (the symmetric
    analog of the H-LU's recursive arith is the remaining work for fully-compressed deep trees).
The rk-aware off-diagonal-LEAF path is verified to machine precision by test_hldlt_factorization.py
(genuine rank-rk off-diagonal leaves); this compact Gram is near-field-heavy so its small-size tree
has dense (not rk) leaves, but the driver + permutation + solve are identical.
"""
import numpy as np

import radia._radia_pybind as _rp

PI = np.pi
INV4PI = 1.0 / (4.0 * PI)
C_CUBE = 1.88231
NEED_RECURSIVE = -5


def _charge_cloud(n=5):
    """Deterministic n^3 structured charge cloud (cell centroids of a unit grid) + dense Gram G
    matching the C++ entry: G[a!=b] = meas_a meas_b/(4pi r), G[a][a] = cube self-energy."""
    h = 1.0 / n
    pts = np.array([[(i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h]
                    for i in range(n) for j in range(n) for k in range(n)])
    N = len(pts)
    meas = np.full(N, h ** 3)
    se = C_CUBE * meas ** (5.0 / 3.0) * INV4PI
    D = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    np.fill_diagonal(D, np.inf)
    G = np.outer(meas, meas) * INV4PI / D
    np.fill_diagonal(G, se)
    return pts, meas, se, G, N


def test_hldlt_factors_real_gram_hmatrix():
    """H-LDL^T factors a real HACApK Gram H-matrix (single leaf AND a 7-leaf diagonal-refined tree)
    and solves G x = b to machine precision."""
    pts, meas, se, G, N = _charge_cloud(5)
    cent, meL, seL = list(pts.ravel()), list(meas), list(se)
    rng = np.random.default_rng(0)
    b = rng.standard_normal(N)
    x_dense = np.linalg.solve(G, b)
    for leaf in (N, N // 2):
        H = _rp._ChargeGramHMatrix(cent, meL, seL, 1e-6, leaf, 2.0)
        res = H.factor_solve_hldlt(b.tolist())
        assert res["ok"], f"H-LDL^T factor failed at leaf={leaf}: factor_rc={res['factor_rc']}"
        x_h = np.asarray(res["x"], float)
        rel = np.linalg.norm(x_h - x_dense) / np.linalg.norm(x_dense)
        assert rel < 1e-11, f"H-LDL^T real-Gram solve rel err {rel:.2e} at leaf={leaf}"


def test_hldlt_deep_tree_fails_loud_need_recursive():
    """A DEEP tree (small leaf -> internal off-diagonal blocks) returns NEED_RECURSIVE cleanly --
    fail-loud, no fallback.  Locks the documented boundary (internal-off-diagonal symmetric recursion
    is future work); the day it lands, this assertion flips and the test above extends to small leaf."""
    pts, meas, se, G, N = _charge_cloud(5)
    cent, meL, seL = list(pts.ravel()), list(meas), list(se)
    b = list(np.ones(N))
    H = _rp._ChargeGramHMatrix(cent, meL, seL, 1e-6, 16, 2.0)
    res = H.factor_solve_hldlt(b)
    assert not res["ok"], "expected deep-tree factor to fail (internal off-diagonal not yet recursed)"
    assert res["factor_rc"] == NEED_RECURSIVE, \
        f"expected NEED_RECURSIVE ({NEED_RECURSIVE}), got factor_rc={res['factor_rc']}"
