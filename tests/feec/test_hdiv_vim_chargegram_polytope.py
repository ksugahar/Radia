"""Golden (Item-1): the C++ POLYTOPE charge-Gram H-matrix (_ChargeGramHMatrix polytope ctor) reproduces
the dense Python analytic_charge_gram POLYTOPE path on HEX and WEDGE meshes -- the hex/wedge analog of
test_hdiv_vim_chargegram_analytic (tet/triangle).  The C++ entry G[a][b] = 0.5*(QuadDot(a,b)+QuadDot(b,a)),
QuadDot(t,s) = (1/4pi) sum_p w_p Phi_s(p) over t's outer quadrature (centroid-fan / Dunavant), Phi_s =
the divergence-theorem polytope potential (cell) / sum-of-sub-triangle Wilton potential (face) -- the SAME
construction (same convex-hull triangulation, same built-in GL4 / Dunavant rules) as
radia.vim._core.analytic_charge_gram, so it must match the dense G entry-by-entry.

Locks (against the dense polytope reference on the SAME hex / wedge mesh):
  (1) all-dense (leaf >= n_charge, no ACA): the polytope H-matvec equals the dense G q to near-MACHINE
      precision -> the C++ polytope entry reproduces the dense analytic_charge_gram polytope entry;
  (2) compressed (ACA): the H-matvec matches the dense G q within the ACA tolerance;
  (3) the hex / wedge cube demag factor through the H-matrix (N = B^T G_h B Rayleigh quotient) is ~1/3.
"""
import math

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
def test_polytope_chargegram_alldense_exact(hexes):
    """leaf >= n_charge -> a single dense leaf, no ACA -> the polytope H-matvec must reproduce the dense
    analytic_charge_gram polytope G to near-machine precision (proves the C++ entry == the dense entry)."""
    with ng.TaskManager():
        d = tet.build_demag(_cube(hexes, 2), analytic_gram=True)
    n = d["n_charge"]
    H = _poly_gram(d, eps=1e-8, leaf=n + 10)
    assert H.ndof() == n
    rng = np.random.default_rng(0)
    q = rng.standard_normal(n)
    y_h = np.asarray(H.matvec(q.tolist()), float)
    y_d = d["G"] @ q
    rel = np.linalg.norm(y_h - y_d) / np.linalg.norm(y_d)
    assert rel < 1e-9, f"{'hex' if hexes else 'wedge'} all-dense polytope H-matvec vs dense G: rel {rel:.2e}"


@pytest.mark.parametrize("hexes", [True, False], ids=["hex", "wedge"])
def test_polytope_chargegram_compressed(hexes):
    """With ACA compression the polytope H-matvec matches the dense polytope G q within the tolerance."""
    with ng.TaskManager():
        d = tet.build_demag(_cube(hexes, 3), analytic_gram=True)
    n = d["n_charge"]
    H = _poly_gram(d, eps=1e-5, leaf=32)
    rng = np.random.default_rng(1)
    q = rng.standard_normal(n)
    y_h = np.asarray(H.matvec(q.tolist()), float)
    y_d = d["G"] @ q
    rel = np.linalg.norm(y_h - y_d) / np.linalg.norm(y_d)
    assert rel < 2e-3, f"{'hex' if hexes else 'wedge'} compressed polytope H-matvec vs dense G: rel {rel:.2e}"


@pytest.mark.parametrize("hexes", [True, False], ids=["hex", "wedge"])
def test_polytope_chargegram_demag(hexes):
    """The cube demag factor through the polytope H-matrix (N = B^T G_h B Rayleigh quotient) is ~1/3 and
    matches the dense polytope demag factor -> the physics survives the C++ polytope H-matrix path."""
    with ng.TaskManager():
        d = tet.build_demag(_cube(hexes, 3), analytic_gram=True)
    H = _poly_gram(d, eps=1e-7, leaf=32)
    B = d["B_csr"]
    m = d["m_unit"]
    Nm = B.T @ np.asarray(H.matvec((B @ m).tolist()), float)
    Dz = float((m @ Nm) / (m @ d["M_mass"] @ m))
    Dz_dense = tet.demag_factor(d)
    assert abs(Dz - Dz_dense) < 5e-3 * abs(Dz_dense) + 1e-3, \
        f"polytope H demag {Dz:.4f} vs dense {Dz_dense:.4f}"
    assert 0.30 < Dz < 0.36, f"{'hex' if hexes else 'wedge'} polytope demag {Dz:.4f} not ~1/3"
