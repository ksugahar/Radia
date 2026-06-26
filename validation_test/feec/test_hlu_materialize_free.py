"""Golden lock: the HACApK H-LU is MATERIALIZE-FREE on the real Coulomb-kernel
collocation system matrix, so it scales sub-cubically instead of degenerating to
the dense materialize-and-redo fallback.

Background: the block-recursive H-LU (cHACApK_harith.c) used to densify whole
internal subtrees in two places -- the addmul "A internal x B internal -> rk C"
Schur update, and the htrsm "dense-leaf L x internal X" triangular solve. On deep
(small-leaf) or unbalanced (large-N) trees this drove the factor to O(N^3.7).
The fix replaces both with materialize-FREE paths: a randomized R(A*B) range-finder
(addmul) and a column/row-son recursion (htrsm). The HLUMaterializeStats n_calls
counter must therefore stay ZERO on these structured RT0 hex grids, and the H-LU
round-trip must stay accurate.

These build a structured n^3 hex cube of MatLin soft iron -> a hex MSC (6 DOF/cell)
collocation interaction matrix, then HLUTestOnHACApK factors A = -N + 1/chi and
round-trips A^-1(A x) = x. Each rebuilds its own HACApK H-matrix at the configured
ACA eps; HLUSetTruncTol sets the H-LU recompression tolerance.
"""
import sys, os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia"))
import radia as rad


def _hex_cube(n, L=0.1, mu_r=50.0):
    rad.UtiDelAll()
    h = L / n
    mat = rad.MatLin(mu_r)
    objs = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x0, y0, z0 = i * h, j * h, k * h
                verts = [[x0, y0, z0], [x0 + h, y0, z0], [x0 + h, y0 + h, z0], [x0, y0 + h, z0],
                         [x0, y0, z0 + h], [x0 + h, y0, z0 + h], [x0 + h, y0 + h, z0 + h], [x0, y0 + h, z0 + h]]
                e = rad.ObjHexahedron(verts, [0.0, 0.0, 0.0])
                rad.MatApl(e, mat)
                objs.append(e)
    return rad.ObjCnt(objs), 6 * n * n * n


def _run_hlu(n, leaf, tol=1e-4, eps=1e-5):
    rad.SolverConfig(hacapk_eps=eps, hacapk_leaf=leaf, hacapk_eta=2.0)
    cnt, ndof = _hex_cube(n)
    handle = rad.BuildMatrix(cnt)
    rad.HLUSetTruncTol(tol)
    err = rad.HLUTestOnHACApK(handle)
    mat = rad.HLUMaterializeStats()
    rad.UtiDelAll()
    # n_internal = internal-subtree densifications (the cubic driver); n_leaf =
    # benign dense-leaf copies for leaf-level dgemms (always present, O(leaf^2)).
    return ndof, err, int(mat["n_internal"])


@pytest.mark.parametrize("n,leaf", [
    (8, 80),    # moderate leaf, shallow tree
    (10, 80),   # 6000 DOF
    (12, 80),   # 10368 DOF -- pre-fix this REappeared (240 materialize calls)
    (8, 16),    # SMALL leaf = deep tree -- pre-fix the heavy-materialize regime
])
def test_hlu_materialize_free(n, leaf):
    """H-LU on the real Coulomb kernel must do ZERO internal-subtree materializations
    (the cubic driver) and stay accurate. Benign dense-leaf copies (n_leaf) are fine."""
    ndof, err, n_internal = _run_hlu(n, leaf)
    assert n_internal == 0, f"H-LU densified {n_internal} internal subtrees (ndof={ndof}, leaf={leaf}) -- the cubic-driving materialize is back"
    assert err < 5e-3, f"H-LU round-trip inaccurate: err={err:.3e} (ndof={ndof}, leaf={leaf})"


def test_hlu_exact_at_tight_tol():
    """At a tight trunc tol the leaf=80 H-LU is a near-exact direct solve (machine
    precision), confirming the materialize-free paths preserve correctness."""
    rad.SolverConfig(hacapk_eps=1e-5, hacapk_leaf=80, hacapk_eta=2.0)
    cnt, ndof = _hex_cube(8)
    handle = rad.BuildMatrix(cnt)
    rad.HLUSetTruncTol(1e-12)
    err = rad.HLUTestOnHACApK(handle)
    mat = rad.HLUMaterializeStats()
    rad.UtiDelAll()
    assert int(mat["n_internal"]) == 0, f"densified {mat['n_internal']} internal subtrees at leaf=80 n=8"
    assert err < 1e-8, f"tight-tol H-LU not near-exact: err={err:.3e}"
