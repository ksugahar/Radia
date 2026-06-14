"""Golden test: the SCALABLE (skip_dense_gram) build_demag path is genuinely sparse -- no dense N^2.

The "scalable" HDiv-VIM Newton used to call build_demag(mesh) with DEFAULT flags, which built the dense
O(N^2) Gram G, the dense demag N, AND the O(N^2) loop SVD, then ignored almost all of it (it drives the
demag apply from the C++ analytic charge-Gram H-matvec instead).  That hidden dense N^2 would OOM at
C-type scale (a 150k-charge dense M_mass is ~180 GB).  This locks the fix:

  build_demag(mesh, skip_dense_gram=True) returns
    - M_mass, B, B_csr  as scipy SPARSE (the RT0 HDiv mass + the L2/SurfaceL2 charge map are LOCAL -> sparse)
    - N, G, loops, n_loop  as None (the dense N^2 objects are never formed)

and the sparse matrices are NUMERICALLY IDENTICAL to the dense reference path (analytic_gram=True), so
nothing downstream changes except the memory/time footprint.  The dense REFERENCE path is unchanged
(still dense, for the small-N demag_factor + dense Newton).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import scipy.sparse as sp  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt  # noqa: E402
from radia.vim import _core as tet  # noqa: E402


def _sphere(h=0.6):
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0)); return ng.Mesh(g.GenerateMesh(maxh=h))


def _cube(h=0.5):
    g = CSGeometry(); g.Add(OrthoBrick(Pnt(-.5, -.5, -.5), Pnt(.5, .5, .5))); return ng.Mesh(g.GenerateMesh(maxh=h))


@pytest.mark.parametrize("mk", [_sphere, _cube])
def test_skip_dense_gram_is_sparse_and_dense_free(mk):
    with ng.TaskManager():
        d = tet.build_demag(mk(), skip_dense_gram=True)
    # the charge map + HDiv mass come back SPARSE
    assert sp.issparse(d["M_mass"]), "skip_dense_gram M_mass must be sparse (would be ~180 GB dense at C-type scale)"
    assert sp.issparse(d["B"]), "skip_dense_gram B must be sparse"
    assert sp.issparse(d["B_csr"]), "B_csr must always be sparse"
    # the dense O(N^2) objects are NEVER formed on the scalable path
    assert d["N"] is None and d["G"] is None and d["loops"] is None and d["n_loop"] is None, \
        "skip_dense_gram must not build the dense N / G / loop-SVD objects"


@pytest.mark.parametrize("mk", [_sphere, _cube])
def test_sparse_build_matches_dense_reference(mk):
    """The sparse scalable build is bit-for-bit the same matrices as the dense reference -- only the
    storage differs.  (Guards against a sparse-assembly bug silently changing the operator.)"""
    with ng.TaskManager():
        dd = tet.build_demag(mk(), analytic_gram=True)        # dense reference path
        ds = tet.build_demag(mk(), skip_dense_gram=True)       # scalable sparse path
    assert isinstance(dd["M_mass"], np.ndarray), "the dense reference path must stay dense"
    assert np.allclose(ds["M_mass"].toarray(), dd["M_mass"], atol=1e-12), "sparse M_mass != dense M_mass"
    assert np.allclose(ds["B"].toarray(), dd["B"], atol=1e-12), "sparse B != dense B"
    assert np.allclose(ds["B_csr"].toarray(), dd["B_csr"].toarray(), atol=1e-12), "B_csr mismatch"
    # the Rayleigh denominator (the only M_mass use shared by both paths) is identical
    mu = dd["m_unit"]
    denom_dense = float(mu @ dd["M_mass"] @ mu)
    denom_sparse = float(mu @ (ds["M_mass"] @ mu))
    assert abs(denom_dense - denom_sparse) <= 1e-12 * abs(denom_dense), "Rayleigh denom drifted"
