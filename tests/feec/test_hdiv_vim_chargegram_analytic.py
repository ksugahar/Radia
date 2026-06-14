"""Golden test (M2b): the C++ ANALYTIC charge-Gram H-matrix (_ChargeGramHMatrix analytic mode)
reproduces the dense Python build_demag(analytic_gram=True) Gram -- i.e. the accurate Wilton/phi_tet
charge Gram now lives in the C++ SCALABLE path, not just the dense Python prototype. The C++ entry is
G[a][b] = 0.5*(QuadDot(a,b)+QuadDot(b,a)), QuadDot(t,s) = (1/4pi) sum_p w_p Phi_s(p) over t's outer
quadrature (tet barycentric sub-points / Dunavant-5), Phi_s = rad_hdiv::PhiTet / TriPotential -- the
same construction as radia.vim._core.analytic_charge_gram, so it must match entry-by-entry.

Locks (against the dense analytic reference on the SAME tet mesh):
  (1) all-dense (leaf >= n_charge, no ACA): the analytic H-matvec equals the dense analytic G q to
      MACHINE precision -> the C++ analytic entry reproduces analytic_charge_gram exactly;
  (2) compressed (ACA): the analytic H-matvec matches the dense analytic G q within the ACA tolerance;
  (3) the sphere demag factor through the analytic H-matrix is ~1/3 and matches the dense analytic one.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

import radia._radia_pybind as _rp  # noqa: E402
from radia.vim import _core as tet  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402


def _sphere(h):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def _analytic_gram(d, eps, leaf):
    return _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]),
                                  face_verts=list(d["face_verts"]),
                                  n_el=int(d["n_el"]), eps=eps, leaf=leaf, eta=2.0)


def test_analytic_chargegram_alldense_exact():
    """leaf >= n_charge -> a single dense leaf, no ACA -> the analytic H-matvec must reproduce the dense
    analytic_charge_gram G to MACHINE precision (proves the C++ entry == the Python analytic entry)."""
    d = tet.build_demag(_sphere(0.8), analytic_gram=True)
    n = d["n_charge"]
    H = _analytic_gram(d, eps=1e-6, leaf=n + 10)
    assert H.ndof() == n
    rng = np.random.default_rng(0)
    q = rng.standard_normal(n)
    y_h = np.asarray(H.matvec(q.tolist()), float)
    y_d = d["G"] @ q
    rel = np.linalg.norm(y_h - y_d) / np.linalg.norm(y_d)
    assert rel < 1e-9, f"all-dense analytic H-matvec vs dense analytic G: rel {rel:.2e}"


def test_analytic_chargegram_compressed():
    """With ACA compression the analytic H-matvec matches the dense analytic G q within the tolerance."""
    d = tet.build_demag(_sphere(0.5), analytic_gram=True)
    n = d["n_charge"]
    H = _analytic_gram(d, eps=1e-5, leaf=32)
    rng = np.random.default_rng(1)
    q = rng.standard_normal(n)
    y_h = np.asarray(H.matvec(q.tolist()), float)
    y_d = d["G"] @ q
    rel = np.linalg.norm(y_h - y_d) / np.linalg.norm(y_d)
    assert rel < 1e-3, f"compressed analytic H-matvec vs dense analytic G: rel {rel:.2e}"


def test_analytic_chargegram_demag_sphere():
    """The sphere demag factor through the analytic H-matrix (N = B^T G_h B Rayleigh quotient) is ~1/3
    and matches the dense analytic demag factor -> the physics survives the analytic H-matrix path."""
    d = tet.build_demag(_sphere(0.5), analytic_gram=True)
    H = _analytic_gram(d, eps=1e-6, leaf=32)
    B = d["B_csr"]
    m = d["m_unit"]
    Nm = B.T @ np.asarray(H.matvec((B @ m).tolist()), float)
    Dz = float((m @ Nm) / (m @ d["M_mass"] @ m))
    Dz_dense = tet.demag_factor(d)
    assert abs(Dz - Dz_dense) < 5e-3 * abs(Dz_dense) + 1e-3, \
        f"analytic H demag {Dz:.4f} vs dense analytic {Dz_dense:.4f}"
    assert 0.30 < Dz < 0.36, f"sphere analytic demag {Dz:.4f} not ~1/3"
