"""Golden test (M3, the NONLINEAR solve in C++): RadHACApKChargeGram::SolveNonlinearPicard runs the
scalar-chi Picard for the isotropic nonlinear demag M = Mof(H0 - Dscal*M) ENTIRELY in C++ -- each
Picard step is a SolveLinearMaterial solve of ((1/chi)M_mass + B^T G B) m = H0*(M_mass mu) with G the
analytic charge-Gram H-matvec, then chi <- 0.5 chi + 0.5 chi_sec(|H|) with the closed-form saturating
curve.  No NGSolve per iteration -- the full nonlinear physics for an isotropic body in C++.

Two locks:
  (1) IMPLEMENTATION: the C++ Picard reproduces a Python scalar-chi Picard driven by the IDENTICAL
      analytic-N operator (same H-matvec, same closed-form chi update) -> same converged M_avg / Dscal;
  (2) PHYSICS: on a sphere the converged M_avg lands on the analytic scalar fixed point
      M = Mof(H0 - Dscal*M) (the same fixed point the Python dense/scalar solvers hit).
"""
import inspect

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import scipy.sparse as sp  # noqa: E402
from scipy.sparse.linalg import LinearOperator, minres  # noqa: E402

import radia._radia_pybind as _rp  # noqa: E402
from radia.vim import _core as tet, _nonlinear as nl  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

_TOLKW = "rtol" if "rtol" in inspect.signature(minres).parameters else "tol"


def _sphere(h=0.6):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_cpp_nonlinear_picard():
    chi0, Msat, H0 = 1000.0, 1.0e6, 3.0e5        # near saturation onset (~Dscal*Msat): genuinely nonlinear
    d = tet.build_demag(_sphere())
    n_face = int(d["ndof"])
    Mm = d["M_mass"]
    mu = d["m_unit"]
    denom = float(mu @ np.asarray(Mm @ mu).ravel())
    H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                               n_el=int(d["n_el"]), eps=1e-8, leaf=32, eta=2.0)
    B = d["B_csr"]
    mco = sp.coo_matrix(Mm)
    Mm_dense = Mm.toarray() if sp.issparse(Mm) else np.asarray(Mm)
    Mmass_diag = np.diag(Mm_dense)
    # the demag N = B^T G B via the C++ analytic charge-Gram H-matvec; its diagonal feeds the Jacobi
    # preconditioner (no dense Python Gram).
    def _Napply0(v):
        v = np.asarray(v, float)
        return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    N_diag = np.array([_Napply0(np.eye(n_face)[:, k])[k] for k in range(n_face)])

    res = H.solve_nonlinear_picard(
        list(map(int, B.indptr)), list(map(int, B.indices)), list(B.data.astype(float)), n_face,
        list(map(int, mco.row)), list(map(int, mco.col)), list(mco.data.astype(float)),
        list(Mmass_diag.astype(float)), list(N_diag.astype(float)), list(mu.astype(float)),
        denom, chi0, Msat, H0, 400, 1e-11, 5000)
    Mavg_cpp, Dscal_cpp = res["Mavg"], res["Dscal"]
    # the 0.5-damped scalar Picard is robust but slow (~230 iters to the 1e-10 break at this drive)
    assert res["iters"] < 400, f"C++ Picard did not converge ({res['iters']} iters)"
    assert 0.0 < Mavg_cpp < Msat, f"C++ Picard M_avg out of range: {Mavg_cpp}"

    # (1) IMPLEMENTATION: a Python scalar-chi Picard on the IDENTICAL analytic-N operator.
    def Napply(v):
        v = np.asarray(v, float)
        return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    Dscal_py = float(mu @ Napply(mu) / denom)
    b0 = Mm @ mu
    chi, Mprev, Mavg_py = chi0, 0.0, 0.0
    Mof = nl._bh_curve(chi0, Msat)
    for it in range(400):
        Aop = LinearOperator((n_face, n_face),
                             matvec=lambda v, c=chi: (1.0 / c) * (Mm @ np.asarray(v, float)) + Napply(v))
        prec = (1.0 / chi) * Mmass_diag + N_diag
        Minv = LinearOperator((n_face, n_face), matvec=lambda v: np.asarray(v, float) / prec)
        mvec, _ = minres(Aop, H0 * b0, M=Minv, maxiter=5000, **{_TOLKW: 1e-11})
        Mavg_py = float(mu @ (Mm @ mvec) / denom)
        Hi = H0 - Dscal_py * Mavg_py
        chi = 0.5 * chi + 0.5 * (chi0 / (1.0 + chi0 * abs(Hi) / Msat))
        if it > 0 and abs(Mavg_py - Mprev) < 1e-10 * (abs(Mavg_py) + 1e-30):
            break
        Mprev = Mavg_py
    assert abs(Dscal_cpp - Dscal_py) < 1e-9 * abs(Dscal_py) + 1e-12, \
        f"Dscal C++ {Dscal_cpp} vs Python {Dscal_py}"
    assert abs(Mavg_cpp - Mavg_py) < 1e-5 * abs(Mavg_py) + 1e-9, \
        f"C++ Picard M_avg {Mavg_cpp} vs Python {Mavg_py}"

    # (2) PHYSICS: sphere -> the analytic scalar fixed point M = Mof(H0 - Dscal*M).
    M_fp = nl._scalar_fixed_point(Mof, Dscal_cpp, H0)
    assert abs(Mavg_cpp - M_fp) < 1e-2 * abs(M_fp) + 1e-6, \
        f"C++ Picard M_avg {Mavg_cpp} vs scalar fixed point {M_fp}"
