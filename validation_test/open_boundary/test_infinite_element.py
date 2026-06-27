# -*- coding: utf-8 -*-
"""Golden tests for radia.infinite_element -- the NGSolve-native static infinite-element open boundary
(the port of the archived Python prototypes act7_32 = assembly +
Steklov spectrum, act7_33 = coupled BVP).  On a sphere the IE == Kelvin (act7_28).

Tiers:
  * PURE numpy: the radial decay operators (well-conditioned vertex+integrated-Legendre basis) give
    the EXACT per-mode DtN -(n+1) for P>=n+1 and are well-conditioned (the act7_28 lesson).
  * NGSOLVE: (a) the condensed DtN matrix has discrete Steklov spectrum == the analytic ladder (n+1)
    with (2n+1) multiplicity; (b) the NGSolve-NATIVE coupled solve (compound space + boundary IE
    terms, monolithic SPARSE, no condensation) reproduces the analytic permeable-sphere demag
    3/(mu_r+2), the exterior dipole, and the (a/r)^2 exterior decay recovered EVERYWHERE via
    exterior_field (the IE is not interior-only).
"""
import numpy as np
import pytest

from radia import infinite_element as ie


# ----------------------------- PURE numpy tier -----------------------------
@pytest.mark.parametrize("n", (0, 1, 2, 3, 4))
def test_radial_per_mode_dtn_exact(n):
    P, a = 6, 1.0
    R1, R0, g = ie.radial_operators(P, a=a)
    E = R1 + n * (n + 1) * R0
    dtn = -1.0 / (g @ np.linalg.solve(E, g))
    assert abs(dtn - (-(n + 1)) / a) < 1e-7, f"n={n}: IE DtN {dtn:.6f} != {-(n+1)/a}"


@pytest.mark.parametrize("a", (1.0, 2.5))
def test_radial_trace_and_conditioning(a):
    P = 6
    R1, R0, g = ie.radial_operators(P, a=a)
    assert abs(g[0] - 1.0) < 1e-13 and np.max(np.abs(g[1:])) < 1e-13, "trace must be e_1"
    cond = np.linalg.cond(R1 + 1 * 2 * R0)   # n=1 energy; monomial basis would be ~1e6
    assert cond < 1e3, f"nodal radial energy cond {cond:.1e} -- should be well-conditioned"


# ----------------------------- NGSOLVE helpers -----------------------------
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


def _solve_permeable_sphere(mu_r, P=6, a=1.0, order=3, h=0.5):
    """NGSolve-native RSP solve: reduced scalar potential phi_red of a permeable sphere (mu_r) in a
    uniform axial field, closed by the IE.  Returns (gf, mesh)."""
    import ngsolve as ng
    from netgen.occ import Sphere, Pnt, OCCGeometry
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), a)).GenerateMesh(maxh=h))
    with ng.TaskManager():
        mesh.Curve(order)
        X = ie.ie_compound_space(mesh, P, order=order)
        trial = X.TrialFunction(); test = X.TestFunction()
        n = ng.specialcf.normal(3)
        a_bf = ng.BilinearForm(X, symmetric=True, check_unused=False)
        a_bf += mu_r * ng.grad(trial[0]) * ng.grad(test[0]) * ng.dx       # interior iron
        ie.add_exterior_ie(a_bf, X, P, a=a)                              # IE exterior closure
        a_bf.Assemble()
        f = ng.LinearForm(X)
        f += (mu_r - 1.0) * n[2] * test[0].Trace() * ng.ds              # (mu_r-1) H0 (n.z), H0=1
        f.Assemble()
        gf = ng.GridFunction(X)
        gf.vec.data = a_bf.mat.Inverse(X.FreeDofs(), inverse="sparsecholesky") * f.vec
    return gf, mesh


# ----------------------------- NGSOLVE tier -----------------------------
def test_dtn_steklov_spectrum_matches_analytic_ladder():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import scipy.linalg as sla
    P, a = 6, 1.0
    MS, KS = _curved_sphere_surface_matrices(p=3, h=0.5, a=a)   # a=1 -> Mtil = M^S
    S = ie.dtn_surface_matrix(MS, KS, P, a=a)
    w = np.sort(sla.eigh(S, MS, eigvals_only=True))
    ladder = []
    for n in range(4):
        ladder += [n + 1] * (2 * n + 1)
    ladder = np.array(ladder, float)
    relerr = np.max(np.abs(w[:len(ladder)] - ladder) / ladder)
    assert relerr < 5e-3, f"Steklov spectrum != analytic ladder (max relerr {relerr:.1e})"
    assert w[0] > 0, "DtN surface operator must be SPD"


@pytest.mark.parametrize("mu_r", (5.0, 50.0, 1000.0))
def test_ngsolve_native_permeable_sphere_interior_and_dipole(mu_r):
    """The NGSolve-native monolithic IE solve reproduces the analytic permeable-sphere demag and
    exterior dipole (== act7_33), with NO condensation (radial levels are explicit DOFs)."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import ngsolve as ng
    a, P = 1.0, 6
    gf, mesh = _solve_permeable_sphere(mu_r, P=P, a=a, order=3, h=0.5)
    grad_an = (mu_r - 1.0) / (mu_r + 2.0)
    C_an = a ** 3 * (mu_r - 1.0) / (mu_r + 2.0)
    gc = ng.grad(gf.components[0])(mesh(0, 0, 0))                       # interior reduced field
    C_fem = (ng.Integrate(gf.components[0] * ng.z, mesh.Boundaries(".*"))
             / ng.Integrate(ng.z * ng.z, mesh.Boundaries(".*")))       # exterior dipole coeff
    assert abs(gc[2] - grad_an) / grad_an < 5e-3, f"interior field {gc[2]:.5f} != {grad_an:.5f}"
    assert (abs(gc[0]) + abs(gc[1])) / grad_an < 5e-3, "transverse field must vanish (axial symmetry)"
    assert abs(C_fem - C_an) / C_an < 5e-3, f"dipole C {C_fem:.5f} != {C_an:.5f}"


def test_exterior_field_recovered_everywhere():
    """exterior_field recovers the analytic dipole C cos(theta)/r^2 at arbitrary EXTERIOR points
    (|x|>a) -- the IE defines the exterior field everywhere, not interior-only."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    a, P, mu_r = 1.0, 6, 50.0
    gf, mesh = _solve_permeable_sphere(mu_r, P=P, a=a, order=3, h=0.5)
    C = a ** 3 * (mu_r - 1.0) / (mu_r + 2.0)
    pts = np.array([(0, 0, 1.5), (0, 0, 3.0), (0, 0, 6.0),
                    (2.0, 0.0, 0.0), (1.5, 0.0, 1.5)], float)
    got = ie.exterior_field(gf, P, a, pts)
    for x, val in zip(pts, got):
        r = np.linalg.norm(x); ct = x[2] / r
        exact = C * ct / r ** 2                                         # phi_red,out = C cos(theta)/r^2
        assert abs(val - exact) < 5e-3 * (abs(exact) + 1e-3), \
            f"exterior field at {tuple(x)}: {val:.5e} vs analytic {exact:.5e}"
