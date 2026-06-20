# -*- coding: utf-8 -*-
"""Golden tests for the C++ static INFINITE-ELEMENT DtN surface operator
(src/core/rad_infinite_element.cpp), the port of the Python prototypes
examples/kelvin_transformation/DtN_spectrum/act7_32 (assembly + Steklov spectrum) and
act7_33 (coupled BVP).  On a sphere the IE == Kelvin (act7_28); these lock the C++ kernel.

Two tiers:
  * PURE (needs only _radia_pybind): the C++ radial decay operators (R1,R0,g) reproduce the
    Python prototype radial_RR to ~1e-12, the trace is clean (g=e_1), and the per-mode DtN
    R1+n(n+1)R0 is EXACT = -(n+1) for P>=n+1 (the IE spectral exactness, act7_25/28).
  * NGSOLVE (importorskip): the condensed DtN surface stiffness S_Gamma on a CURVED sphere has
    discrete Steklov spectrum eig(S, M^S) == the analytic ladder (n+1) with multiplicity (2n+1)
    (the act7_32 result, now produced by the C++ kernel), and == the Python prototype.
"""
import numpy as np
import pytest

import radia._radia_pybind as _rp
from radia import infinite_element as ie


# ----------------------------- Python prototype reference (act7_32 radial_RR) -----------------------------
def _gauss01(nq):
    x, w = np.polynomial.legendre.leggauss(nq)
    return 0.5 * (x + 1.0), 0.5 * w


def _legval(j, xi):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(xi, c)


def _radial_eval(P, t):
    t = np.asarray(t, float)
    N = np.zeros((P, t.size)); Np = np.zeros((P, t.size))
    N[0] = t; Np[0] = np.ones_like(t)
    xi = 2.0 * t - 1.0
    for k in range(2, P + 1):
        N[k - 1] = (_legval(k, xi) - _legval(k - 2, xi)) / (2.0 * k - 1.0)
        Np[k - 1] = _legval(k - 1, xi) * 2.0
    return N, Np


def _radial_RR_py(P, a=1.0, nq=160):
    t, w = _gauss01(nq)
    N, Np = _radial_eval(P, t)
    return a * (Np * w) @ Np.T, a * (N / t ** 2 * w) @ N.T


# ----------------------------- PURE tier -----------------------------
@pytest.mark.parametrize("a", (1.0, 2.5))
def test_cpp_radial_operators_match_python(a):
    P = 6
    R1c, R0c, g = ie.radial_operators(P, a=a)
    R1p, R0p = _radial_RR_py(P, a=a)
    assert np.max(np.abs(R1c - R1p)) < 1e-11, "C++ R1 != Python radial_RR"
    assert np.max(np.abs(R0c - R0p)) < 1e-11, "C++ R0 != Python radial_RR"
    assert abs(g[0] - 1.0) < 1e-13 and np.max(np.abs(g[1:])) < 1e-13, "trace g must be e_1 (clean trace)"


@pytest.mark.parametrize("n", (0, 1, 2, 3, 4))
def test_cpp_radial_per_mode_dtn_exact(n):
    """Per-mode DtN from the C++ radial operators is exact -(n+1) once P>=n+1 (IE spectral)."""
    P = 6
    a = 1.0
    R1, R0, g = ie.radial_operators(P, a=a)
    E = R1 + n * (n + 1) * R0
    dtn = -1.0 / (g @ np.linalg.solve(E, g))
    assert abs(dtn - (-(n + 1)) / a) < 1e-7, f"n={n}: C++ IE DtN {dtn:.6f} != {-(n+1)/a}"


def test_cpp_radial_nodal_well_conditioned():
    """The nodal radial basis is well-conditioned (the act7_28 lesson; monomial would be ~1e6)."""
    P = 6
    R1, R0, _ = ie.radial_operators(P, a=1.0)
    cond = np.linalg.cond(R1 + 1 * 2 * R0)  # n=1 energy
    assert cond < 1e3, f"nodal radial energy cond {cond:.1e} -- should be well-conditioned"


# ----------------------------- NGSOLVE tier -----------------------------
def _curved_sphere_surface_matrices(p=3, h=0.5, a=1.0):
    import ngsolve as ng
    from netgen.occ import Sphere, Pnt, OCCGeometry
    import scipy.sparse as sp
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=h))
    with ng.TaskManager():
        mesh.Curve(p)
        fes = ng.H1(mesh, order=p)
        u, v = fes.TnT()
        n = ng.specialcf.normal(3)
        bm = ng.BilinearForm(fes, symmetric=True, check_unused=False); bm += u * v * ng.ds; bm.Assemble()
        gu, gv = ng.grad(u).Trace(), ng.grad(v).Trace()
        gut = gu - (gu * n) * n; gvt = gv - (gv * n) * n
        bk = ng.BilinearForm(fes, symmetric=True, check_unused=False); bk += (gut * gvt) * ng.ds; bk.Assemble()

    def _csr(m):
        r, c, val = m.COO()
        return sp.csr_matrix((np.asarray(val), (np.asarray(r), np.asarray(c))), shape=(m.height, m.height))

    bnd = np.array([i for i in range(fes.ndof) if fes.GetDofs(mesh.Boundaries(".*"))[i]], dtype=int)
    MS = _csr(bm.mat)[np.ix_(bnd, bnd)].toarray(); MS = 0.5 * (MS + MS.T)
    KS = _csr(bk.mat)[np.ix_(bnd, bnd)].toarray(); KS = 0.5 * (KS + KS.T)
    return MS, KS


def test_cpp_dtn_steklov_spectrum_matches_analytic_ladder():
    """The C++ condensed DtN surface operator reproduces the analytic Steklov ladder (n+1)
    with the correct (2n+1) multiplicities -- the act7_32 result via the C++ kernel."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import scipy.linalg as sla
    P, a = 6, 1.0
    MS, KS = _curved_sphere_surface_matrices(p=3, h=0.5, a=a)   # a=1 -> Mtil = M^S
    S = ie.dtn_surface_operator(MS, KS, P, a=a)
    assert S.shape == MS.shape
    w = np.sort(sla.eigh(S, MS, eigvals_only=True))
    ladder = []
    for n in range(4):
        ladder += [n + 1] * (2 * n + 1)
    ladder = np.array(ladder, float)
    relerr = np.max(np.abs(w[:len(ladder)] - ladder) / ladder)
    assert relerr < 5e-3, f"C++ Steklov spectrum != analytic ladder (max relerr {relerr:.1e})"
    assert w[0] > 0, "DtN surface operator must be SPD"


def test_cpp_matches_python_prototype_spectrum():
    """C++ S_Gamma spectrum == the Python prototype ie_surface_operator spectrum (same kernel)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import scipy.linalg as sla
    P, a = 5, 1.0
    MS, KS = _curved_sphere_surface_matrices(p=2, h=0.6, a=a)

    # Python prototype condensation (act7_32 schur_blocks)
    R1, R0 = _radial_RR_py(P, a)
    block = lambda k, l: R1[k, l] * MS + R0[k, l] * KS
    A11 = block(0, 0)
    b = list(range(1, P))
    A1b = np.hstack([block(0, l) for l in b])
    Abb = np.vstack([np.hstack([block(k, l) for l in b]) for k in b])
    S_py = A11 - A1b @ np.linalg.solve(Abb, A1b.T)
    S_py = 0.5 * (S_py + S_py.T)

    S_cpp = ie.dtn_surface_operator(MS, KS, P, a=a)
    w_py = np.sort(sla.eigh(S_py, MS, eigvals_only=True))
    w_cpp = np.sort(sla.eigh(S_cpp, MS, eigvals_only=True))
    assert np.max(np.abs(w_py - w_cpp)) < 1e-8, "C++ != Python prototype condensation"


def test_surface_dtn_from_mesh_runs():
    """The convenience closure surface_dtn_from_mesh returns an SPD S of the right size."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import ngsolve as ng
    import scipy.linalg as sla
    from netgen.occ import Sphere, Pnt, OCCGeometry
    a = 1.0
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=0.6))
    with ng.TaskManager():
        mesh.Curve(2)
        S, bnd = ie.surface_dtn_from_mesh(mesh, P=5, a=a, order=2)
    assert S.shape == (len(bnd), len(bnd))
    assert np.all(np.isfinite(S))
    assert np.min(sla.eigvalsh(S)) > -1e-8, "S must be (numerically) SPD"
