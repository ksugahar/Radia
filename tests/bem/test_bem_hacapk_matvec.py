"""BEM HACApK compressed-MatVec accuracy regression.

The 2026-06-26 HACApK-subclass audit found that `RadHACApKBEMManager`
(the in-tree BEM H-matrix wrapper, pybind `HACApKBEMManager`) had NO test
of its COMPRESSED MatVec: the existing BEM tests cover the dense Galerkin
SL/DL assembly (`tests/bem/test_sibc_hacapk_match_ngsbem.py`) and an
end-to-end heating solve, but neither asserts that the ACA-compressed
H-matrix MatVec reproduces the dense `A @ x` -- so a subtle ACA-compression
corruption could land silently.

`HACApKBEMManager` is geometry/kernel-agnostic: it takes (coords, dense)
and wraps the dense matrix as an ACA H-matrix. The at-risk code is the ACA
compression + MatVec, NOT the SL assembly (separately tested). So this test
feeds it the genuine BEM single-layer Laplace kernel `1/(4 pi r)` over a
PLANAR point patch -- a flat patch separates clusters so admissible far
blocks become genuinely low-rank and the ACA path is actually exercised
(a compact sphere at test scale stays all-dense and would not engage ACA).
Fast: pure numpy + radia, no NGSolve, no O(N^2) singular quadrature.
"""
import numpy as np
import pytest

_pb = pytest.importorskip("radia._radia_pybind")
if not hasattr(_pb, "HACApKBEMManager"):
    pytest.skip("radia built without HACApKBEMManager", allow_module_level=True)

INV4PI = 1.0 / (4.0 * np.pi)


def _sl_kernel_matrix(coords, h):
    """Symmetric BEM single-layer Laplace matrix A[i,j] = 1/(4 pi |xi-xj|),
    self term regularized by the panel size h."""
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    np.fill_diagonal(d, h)
    return np.ascontiguousarray(INV4PI / d)


def _planar_patch(n_side, L=1.0):
    g = np.linspace(0.0, L, n_side)
    xx, yy = np.meshgrid(g, g)
    pts = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(xx.size)])
    return np.ascontiguousarray(pts, dtype=np.float64), L / (n_side - 1)


def _build(coords, A, aca_eps, leaf=16, eta=2.0):
    mgr = _pb.HACApKBEMManager(coords, A)
    ok = mgr.BuildHMatrix(aca_eps=aca_eps, leaf_size=leaf, eta=eta,
                          max_rank=-1, print_level=0)
    assert ok and mgr.IsValid(), "HACApKBEMManager.BuildHMatrix failed"
    return mgr


def test_bem_hacapk_engages_low_rank():
    """Guard the test's own validity: on a planar patch the ACA path MUST
    produce low-rank blocks (else the matvec test below is vacuously dense)."""
    coords, h = _planar_patch(16)          # N = 256
    A = _sl_kernel_matrix(coords, h)
    mgr = _build(coords, A, aca_eps=1e-6)
    st = mgr.GetStats()
    assert mgr.GetNDOF() == coords.shape[0]
    assert st["n_lowrank"] > 0, (
        f"ACA produced no low-rank blocks (n_lowrank={st['n_lowrank']}); "
        "the compressed-matvec path is not being exercised")


@pytest.mark.parametrize("aca_eps", [1e-4, 1e-6, 1e-8])
def test_bem_hacapk_matvec_matches_dense(aca_eps):
    """Compressed H-matrix MatVec reproduces the dense A @ x within a few x
    aca_eps. Catches ACA-compression corruption that the SL-assembly and
    e2e heating tests cannot see."""
    coords, h = _planar_patch(16)          # N = 256
    A = _sl_kernel_matrix(coords, h)
    mgr = _build(coords, A, aca_eps=aca_eps)

    rng = np.random.default_rng(0)
    rel_max = 0.0
    for _ in range(4):
        x = np.ascontiguousarray(rng.standard_normal(coords.shape[0]))
        y_h = mgr.MatVec(x)
        y_d = A @ x
        assert np.all(np.isfinite(y_h))
        rel_max = max(rel_max, np.linalg.norm(y_h - y_d) / np.linalg.norm(y_d))

    tol = 50.0 * aca_eps + 1e-9
    assert rel_max < tol, (
        f"compressed matvec rel err {rel_max:.3e} exceeds {tol:.3e} "
        f"(aca_eps={aca_eps:.0e})")


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
