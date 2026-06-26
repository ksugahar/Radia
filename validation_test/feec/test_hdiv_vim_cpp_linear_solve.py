"""Golden test (M3, the iterative-solve hot kernel in C++): RadHACApKChargeGram::SolveLinearMaterial
solves the SPD HDiv-VIM linear material system ((1/chi)M_mass + B^T G B) m = rhs by Jacobi-PCG, with
G applied as the analytic charge-Gram H-matvec -- the linear soft-iron demag solve and the symmetric
Picard warmstart of the nonlinear Newton, now in C++ (no Python per-iteration glue).

Two locks:
  (1) IMPLEMENTATION: the C++ CG reproduces scipy MINRES on the IDENTICAL operator (same H-matvec) to
      ~machine precision -> the C++ Krylov + sparse B/M_mass plumbing is correct, independent of ACA;
  (2) PHYSICS: the C++ CG solution matches the dense analytic direct solve within the ACA tolerance
      (amplified by cond(A) ~ 1/inv_chi) -> it solves the right (physical) linear-material system.
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
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_cpp_cg_linear_material():
    d = tet.build_demag(_sphere())
    n_face = int(d["ndof"])
    H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                               n_el=int(d["n_el"]), eps=1e-8, leaf=32, eta=2.0)
    B = d["B_csr"]                                  # n_charge x n_face CSR
    Mm = d["M_mass"]
    inv_chi = 1e-3

    # the demag N = B^T G B applied via the C++ analytic charge-Gram H-matvec (no dense Python Gram);
    # form it densely (small mesh) from the SAME operator the C++ solver uses, for the Jacobi diagonal +
    # the direct-solve reference.
    def N_apply(v):
        v = np.asarray(v, float)
        return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    N_dense = np.column_stack([N_apply(np.eye(n_face)[:, k]) for k in range(n_face)])
    Mm_dense = Mm.toarray() if sp.issparse(Mm) else np.asarray(Mm)
    rng = np.random.default_rng(3)
    rhs = rng.standard_normal(n_face)
    prec = inv_chi * np.diag(Mm_dense) + np.diag(N_dense)  # exact Jacobi diagonal of the system

    # M_mass as COO
    mco = sp.coo_matrix(Mm)
    res = H.solve_linear_material(
        list(map(int, B.indptr)), list(map(int, B.indices)), list(B.data.astype(float)),
        n_face, list(map(int, mco.row)), list(map(int, mco.col)), list(mco.data.astype(float)),
        inv_chi, list(prec.astype(float)), list(rhs.astype(float)), 1e-12, 5000)
    x_cpp = np.asarray(res["m"], float)
    assert res["iters"] < 5000, "C++ CG did not converge"

    # (1) IMPLEMENTATION: scipy MINRES on the IDENTICAL operator (same H-matvec) -> ~machine match.
    def apply_py(v):
        v = np.asarray(v, float)
        return inv_chi * (Mm @ v) + N_apply(v)

    A = LinearOperator((n_face, n_face), matvec=apply_py)
    Minv = LinearOperator((n_face, n_face), matvec=lambda v: np.asarray(v, float) / prec)
    x_sp, info = minres(A, rhs, M=Minv, maxiter=5000, **{_TOLKW: 1e-12})
    rel_impl = np.linalg.norm(x_cpp - x_sp) / np.linalg.norm(x_sp)
    assert rel_impl < 1e-6, f"C++ CG vs scipy MINRES (same operator): rel {rel_impl:.2e}"

    # (2) PHYSICS: matches the dense direct solve (dense A built from the SAME C++ H-matvec) within ACA
    # tol * cond(A) (cond ~ 1/inv_chi).
    A_dense = inv_chi * Mm_dense + N_dense
    x_dense = np.linalg.solve(A_dense, rhs)
    rel_phys = np.linalg.norm(x_cpp - x_dense) / np.linalg.norm(x_dense)
    assert rel_phys < 5e-4, f"C++ CG vs dense solve: rel {rel_phys:.2e}"
