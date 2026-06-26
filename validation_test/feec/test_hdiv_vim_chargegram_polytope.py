"""Golden (Item-1): the C++ POLYTOPE charge-Gram H-matrix (_ChargeGramHMatrix polytope ctor) on HEX and
WEDGE meshes -- the hex/wedge analog of the tet/triangle charge Gram.  The C++ entry
G[a][b] = 0.5*(QuadDot(a,b)+QuadDot(b,a)), QuadDot(t,s) = (1/4pi) sum_p w_p Phi_s(p) over t's outer
quadrature (centroid-fan / Dunavant), Phi_s = the divergence-theorem polytope potential (cell) /
sum-of-sub-triangle flat-triangle potential (face).  The dense Python Gram path was removed, so this
locks the PHYSICS the C++ polytope Gram must deliver: the hex / wedge cube demag factor through the
H-matrix (N = B^T G_h B Rayleigh quotient) is ~1/3 (isotropic uniform-M body).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")

import radia._radia_pybind as _rp  # noqa: E402
from radia.vim import _core as tet  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

L = 0.02


def _cube(hexes, n):
    mp = lambda x, y, z: (L * x, L * y, L * z)  # noqa: E731
    with ng.TaskManager():
        return MakeStructured3DMesh(hexes=hexes, prism=(not hexes), nx=n, ny=n, nz=n, mapping=mp)


def _poly_gram(d, eps, leaf):
    p = d["poly"]
    return _rp._ChargeGramHMatrix(
        cell_tris=list(p["cell_tris"]), cell_troff=list(p["cell_troff"]),
        cell_cent=list(p["cell_cent"]), cell_meas=list(p["cell_meas"]),
        face_tris=list(p["face_tris"]), face_troff=list(p["face_troff"]),
        face_cent=list(p["face_cent"]), face_meas=list(p["face_meas"]),
        n_el=int(d["n_el"]), eps=eps, leaf=leaf, eta=2.0)


@pytest.mark.parametrize("hexes", [True, False], ids=["hex", "wedge"])
def test_polytope_chargegram_demag(hexes):
    """The cube demag factor through the C++ polytope charge-Gram H-matrix (N = B^T G_h B Rayleigh
    quotient) is ~1/3 -- the hex / wedge physics gate."""
    with ng.TaskManager():
        d = tet.build_demag(_cube(hexes, 3))
    H = _poly_gram(d, eps=1e-7, leaf=32)
    B = d["B_csr"]
    m = d["m_unit"]
    Nm = B.T @ np.asarray(H.matvec((B @ m).tolist()), float)
    Dz = float((m @ Nm) / (m @ (d["M_mass"] @ m)))
    assert 0.30 < Dz < 0.36, f"{'hex' if hexes else 'wedge'} polytope demag {Dz:.4f} not ~1/3"
