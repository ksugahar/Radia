"""Golden test: build_demag is SPARSE-ONLY -- no dense N^2.

The dense Python charge-Gram path (the O(N^2) Gram G, the dense demag N, the O(N^2) loop SVD) was
REMOVED; the C++ analytic _ChargeGramHMatrix is the demag operator.  A hidden dense N^2 would OOM at
C-type scale (a 150k-charge dense M_mass is ~180 GB).  This locks that build_demag(mesh) returns
  - M_mass, B, B_csr  as scipy SPARSE (the RT0 HDiv mass + the L2/SurfaceL2 charge map are LOCAL -> sparse),
  - NO dense N / G / loops / n_loop keys (those objects are never formed).
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
def test_build_demag_is_sparse_and_dense_free(mk):
    with ng.TaskManager():
        d = tet.build_demag(mk())
    # the charge map + HDiv mass come back SPARSE
    assert sp.issparse(d["M_mass"]), "build_demag M_mass must be sparse (would be ~180 GB dense at C-type scale)"
    assert sp.issparse(d["B"]), "build_demag B must be sparse"
    assert sp.issparse(d["B_csr"]), "B_csr must always be sparse"
    # the dense O(N^2) objects are NEVER formed -- the dense Python Gram path was removed
    for k in ("N", "G", "loops", "n_loop"):
        assert k not in d, f"build_demag must not return the dense '{k}' object (removed dense Gram path)"
    # the C++ charge-Gram geometry IS present (cell/face verts for tet, poly for hex/wedge)
    assert "cell_verts" in d and "face_verts" in d and "poly" in d and "n_el" in d
