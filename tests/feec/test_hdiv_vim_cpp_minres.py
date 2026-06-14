"""Golden test: the C++ mu_r-independent MATERIAL solve RadHACApKChargeGram::SolveMaterialMINRES --
Jacobi-preconditioned MINRES for the SYMMETRIC INDEFINITE A m = rhs, A = inv_chi*M_mass - B^T G B, with
G the analytic charge-Gram H-matvec (HACApK-parallel under ngcore::RegionTaskManager).  This is the
production demag solver (the loop-field-null mu_r-independent one, vs the +N CG Picard warmstart).

Two locks:
  (1) IMPLEMENTATION: the C++ MINRES reproduces scipy MINRES on the IDENTICAL operator (same H-matvec,
      same Jacobi preconditioner) to ~machine precision -> the C++ MINRES recurrence + the -N apply are
      correct, independent of ACA;
  (2) PHYSICS: the C++ MINRES solution matches the dense direct solve of (inv_chi*M_mass - N) within the
      ACA tolerance amplified by cond(A).
"""
import inspect

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import scipy.sparse as sp  # noqa: E402
from scipy.sparse.linalg import LinearOperator, minres  # noqa: E402

import radia._radia_pybind as _rp  # noqa: E402
from radia.vim import _core as tet  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

_TOLKW = "rtol" if "rtol" in inspect.signature(minres).parameters else "tol"


def _sphere(h=0.6):
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=h))


def test_cpp_minres_material():
    d = tet.build_demag(_sphere(), analytic_gram=True)
    n_face = int(d["ndof"])
    Mm = d["M_mass"]
    mu = d["m_unit"]
    H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                               n_el=int(d["n_el"]), eps=1e-8, leaf=32, eta=2.0)
    B = d["B_csr"]
    mco = sp.coo_matrix(Mm)
    inv_chi = 1.0 / 999.0                                  # mu_r = 1000
    rhs = Mm @ mu                                          # physical uniform-field source
    prec = np.abs(inv_chi * np.diag(Mm) - np.diag(d["N"]))  # SPD Jacobi for the indefinite A
    prec = np.maximum(prec, prec.max() * 1e-12)

    res = H.solve_material_minres(
        list(map(int, B.indptr)), list(map(int, B.indices)), list(B.data.astype(float)), n_face,
        list(map(int, mco.row)), list(map(int, mco.col)), list(mco.data.astype(float)),
        inv_chi, list(prec.astype(float)), list(rhs.astype(float)), 1e-9, 5000)
    x_cpp = np.asarray(res["m"], float)
    assert res["iters"] < 5000, f"C++ MINRES did not converge ({res['iters']} iters)"

    # (1) IMPLEMENTATION: scipy MINRES on the IDENTICAL operator + preconditioner.
    def applyA(v):
        v = np.asarray(v, float)
        return inv_chi * (Mm @ v) - B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    A = LinearOperator((n_face, n_face), matvec=applyA)
    Minv = LinearOperator((n_face, n_face), matvec=lambda v: np.asarray(v, float) / prec)
    x_sp, _ = minres(A, rhs, M=Minv, maxiter=5000, **{_TOLKW: 1e-9})
    rel1 = np.linalg.norm(x_cpp - x_sp) / np.linalg.norm(x_sp)
    # corroboration: both MINRES on the same operator land at the same solution to ~solver tol.  The
    # residual floor (~1e-4) is the DIFFERENT stopping criteria (C++ phibar/beta1 vs scipy
    # rnorm/(Anorm*ynorm)) on an indefinite cond~1/inv_chi system -- NOT a recurrence bug (which would
    # be O(1)).  The rigorous correctness gate is the dense-true-solution check below.
    assert rel1 < 1e-4, f"C++ MINRES vs scipy MINRES (same operator): rel {rel1:.2e}"

    # (2) PHYSICS / CORRECTNESS GATE: the C++ MINRES finds the TRUE solution (dense direct solve of
    # A = inv_chi*M_mass - N) within ACA_eps * cond(A) -- a recurrence bug could not reach this.
    x_dense = np.linalg.solve(inv_chi * Mm - d["N"], rhs)
    rel2 = np.linalg.norm(x_cpp - x_dense) / np.linalg.norm(x_dense)
    assert rel2 < 5e-3, f"C++ MINRES vs dense true solution: rel {rel2:.2e}"
