"""Golden test: production step #1b -- the C++ HACApK charge-Gram H-matrix (_ChargeGramHMatrix) on
UNSTRUCTURED TET meshes.  This is the SCALABLE half of #1: the dense O(n_charge^2) Coulomb Gram G
of the tet ingest (#1a) is replaced by a HACApK H-matrix (O(N log N) matvec), and the demag operator
N = B^T G B is applied as B^T (G-Hmatvec (B m)) with B the sparse charge map from NGSolve HDiv(0).

Locks (against the dense #1a reference on the SAME tet mesh):
  - the Gram H-matvec reproduces the dense G q within the ACA tolerance;
  - the demag factor via the H-matrix (N_h = B^T G_h B Rayleigh quotient) matches the dense demag
    factor -> the PHYSICS survives the H-matrix compression on real tets;
  - a MINRES solve of A = (1/chi) M_mass - N driven by the Gram H-matvec matches the dense direct
    solve -> the scalable tet solver is correct.
"""
import os
import sys
import inspect

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
minres = pytest.importorskip("scipy.sparse.linalg").minres
from scipy.sparse.linalg import LinearOperator  # noqa: E402

import radia._radia_pybind as _rp  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
from radia.vim import _core as tet  # noqa: E402

import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

_TOLKW = "rtol" if "rtol" in inspect.signature(minres).parameters else "tol"


def _sphere(h=0.45):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def _gram(d, eps=1e-5):
    return _rp._ChargeGramHMatrix(list(d["cent"].ravel()), list(d["meas"]),
                                  list(d["self_energy"]), eps, 32, 2.0)


def test_charge_gram_hmatvec_matches_dense():
    """The charge-Gram H-matvec reproduces the dense Coulomb Gram G q within the ACA tolerance."""
    d = tet.build_demag(_sphere(), nsub=4)
    H = _gram(d, eps=1e-5)
    assert H.ndof() == d["n_charge"]
    rng = np.random.default_rng(0)
    q = rng.standard_normal(d["n_charge"])
    y_h = np.asarray(H.matvec(q.tolist()), float)
    y_d = d["G"] @ q
    rel = np.linalg.norm(y_h - y_d) / np.linalg.norm(y_d)
    assert rel < 1e-3, f"charge-Gram H-matvec rel err {rel:.2e}"


def test_tet_demag_factor_via_hmatrix_matches_dense():
    """The demag factor computed through the H-matrix (N_h = B^T G_h B) matches the dense one --
    the tet physics (sphere demag -> 1/3) survives the H-matrix compression."""
    d = tet.build_demag(_sphere(), nsub=4)
    H = _gram(d, eps=1e-6)
    B = d["B_csr"]
    m = d["m_unit"]
    Nm_h = B.T @ np.asarray(H.matvec((B @ m).tolist()), float)
    Dz_h = float((m @ Nm_h) / (m @ d["M_mass"] @ m))
    Dz_dense = tet.demag_factor(d)
    assert abs(Dz_h - Dz_dense) < 5e-3 * abs(Dz_dense) + 1e-3, \
        f"H-matrix demag {Dz_h:.4f} vs dense {Dz_dense:.4f}"
    assert 0.28 < Dz_h < 0.34, f"H-matrix sphere demag {Dz_h:.4f} not ~1/3"


def test_tet_minres_solve_via_hmatrix_matches_dense():
    """A MINRES solve of A = (1/chi) M_mass - N driven by the Gram H-matvec (N m = B^T G_h(B m))
    matches the dense direct solve -> the scalable unstructured-tet solver is correct."""
    d = tet.build_demag(_sphere(), nsub=4)
    H = _gram(d, eps=1e-7)   # tight Gram: the solve error is ~ ACA_eps * cond(A); cond ~ 1/inv_chi
    B = d["B_csr"]
    ndof = d["ndof"]
    inv_chi = 1e-3
    A_dense = inv_chi * d["M_mass"] - d["N"]
    rng = np.random.default_rng(1)
    b = rng.standard_normal(ndof)
    x_dense = np.linalg.solve(A_dense, b)

    def apply(v):
        v = np.asarray(v, float)
        Nv = B.T @ np.asarray(H.matvec((B @ v).tolist()), float)
        return inv_chi * (d["M_mass"] @ v) - Nv

    A = LinearOperator((ndof, ndof), matvec=apply)
    diag = np.maximum(np.abs(inv_chi * np.diag(d["M_mass"]) - np.diag(d["N"])), 1e-30)
    Minv = LinearOperator((ndof, ndof), matvec=lambda v: np.asarray(v, float) / diag)
    x_h, info = minres(A, b, M=Minv, maxiter=5000, **{_TOLKW: 1e-9})
    rel = np.linalg.norm(x_h - x_dense) / np.linalg.norm(x_dense)
    # the residual ||A_dense x_h - b|| confirms x_h solves the (slightly H-approximated) system; the
    # solution-vector rel err is ACA_eps amplified by cond(A) ~ 1/inv_chi (loops sit at eig=inv_chi).
    assert rel < 5e-3, f"tet H-matrix MINRES vs dense rel err {rel:.2e} (info={info})"
