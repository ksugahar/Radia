"""Charge-basis normalization contract (the roundoff-amplification fix).

The ChargeGram H-matrix stores the normalized Gram Ghat = S^-1 G S^-1
(sigma_p = sqrt(raw G_pp)) and wraps S back inside every public apply, so
ALL external semantics remain the physical Gram.  These tests lock that
contract on a small sphere:

* ``charge_sigma`` is populated by the build, strictly positive, and
  genuinely non-uniform (it IS the per-charge scale);
* ``entry(i, j)`` keeps physical values: the diagonal reproduces
  ``sigma_p**2`` and the matrix stays symmetric;
* ``matvec_sym`` equals the dense physical Gram assembled from
  ``entry`` -- i.e. the internal normalization is invisible outside;
* the assembled demag operator keeps its physical spectrum band.

Why this exists: on extreme element-size-ratio meshes (marching-cubes
micro-facet blobs, 3700x volume ratio) the unnormalized basis let the
congruence B^T G B amplify G's absolute rounding band by ||B||^2 ~ 1e16
into O(1) NEGATIVE demag eigenvalues (bug pattern
facet-tet-charge-gram-indefinite-cg-stall).
"""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

from netgen.occ import OCCGeometry, Sphere, Pnt  # noqa: E402

from radia import vim  # noqa: E402


@pytest.fixture(scope="module")
def sphere_gram():
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(
            Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.55))
        fes = ng.HDiv(mesh, order=1)
        operator = vim.DemagOperator(fes, eps=1e-10)
    return mesh, fes, operator


def test_charge_sigma_is_populated_and_positive(sphere_gram):
    _, _, operator = sphere_gram
    sigma = np.asarray(operator._G.charge_sigma(), dtype=float)
    n_charge = operator._B.shape[0]
    assert sigma.shape == (n_charge,)
    assert np.all(sigma > 0.0) and np.all(np.isfinite(sigma))
    # the scale is real, not a constant no-op
    assert sigma.max() / sigma.min() > 1.5


def test_entry_stays_physical_and_matches_sigma_diagonal(sphere_gram):
    _, _, operator = sphere_gram
    G = operator._G
    sigma = np.asarray(G.charge_sigma(), dtype=float)
    n = sigma.size
    for p in (0, n // 3, n - 1):
        assert G.entry(p, p) == pytest.approx(sigma[p] ** 2, rel=1e-12)
    rng = np.random.default_rng(11)
    for _ in range(12):
        i, j = rng.integers(0, n, size=2)
        assert G.entry(int(i), int(j)) == pytest.approx(
            G.entry(int(j), int(i)), rel=1e-11, abs=1e-300)


def test_matvec_sym_equals_dense_physical_gram(sphere_gram):
    _, _, operator = sphere_gram
    G = operator._G
    sigma = np.asarray(G.charge_sigma(), dtype=float)
    n = sigma.size
    dense = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            dense[i, j] = G.entry(i, j)
    rng = np.random.default_rng(7)
    for _ in range(3):
        x = rng.standard_normal(n)
        y = np.asarray(G.matvec_sym(x), dtype=float)
        # eps=1e-10 ACA build vs the exact entry oracle
        assert np.linalg.norm(y - dense @ x) <= 1e-8 * np.linalg.norm(
            dense @ x) + 1e-16


def test_demag_operator_keeps_physical_spectrum_band(sphere_gram):
    mesh, fes, operator = sphere_gram
    import scipy.linalg as sla
    import scipy.sparse as sp

    with ng.TaskManager():
        u, v = fes.TnT()
        mass = ng.BilinearForm(fes, symmetric=True)
        mass += u * v * ng.dx
        mass.Assemble()
        n = fes.ndof
        x = ng.GridFunction(fes)
        y = x.vec.CreateVector()
        Nd = np.empty((n, n))
        for i in range(n):
            x.vec[:] = 0.0
            x.vec[i] = 1.0
            operator.mat.Mult(x.vec, y)
            Nd[:, i] = y.FV().NumPy()
        rows, cols, vals = mass.mat.COO()
        Md = sp.coo_matrix(
            (np.array(vals), (np.array(rows), np.array(cols))),
            shape=(n, n)).toarray()
    w = sla.eigh(0.5 * (Nd + Nd.T), Md, eigvals_only=True)
    # physical magnetization-operator band [0, 1] with discretization slack;
    # a uniform sphere sits near 1/3
    assert w[0] > -1e-8
    assert w[-1] < 1.5
