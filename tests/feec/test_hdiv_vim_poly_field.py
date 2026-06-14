"""Golden test: the POLYNOMIAL-CHARGE field kernel (reconstruct_field_polynomial) -- the order>=2 field.

The centroid `reconstruct_field` is only correct for piecewise-constant M (surface charge sigma = M.n).
For a polynomial M the VOLUME charge rho = -div M is non-zero, and its omission is a 90-230% error.
This kernel adds the volume charge:

    H(r) = (1/4pi)[ INT_V (-div M)(r-r')/|r-r'|^3 dV' + INT_S (M.n)(r-r')/|r-r'|^3 dS' ]

Three locks:
  (1) uniform M sphere (div M = 0): center field == -M/3 (the surface charge alone, exact);
  (2) uniform M sphere: external field == the analytic point dipole (modulo the flat-mesh faceting);
  (3) linear M (div M = const != 0): the volume charge is ESSENTIAL -- dropping rho (== the old
      centroid/surface-only reconstruct_field) is >50% wrong vs a finer-mesh reference, while the full
      kernel is coarse->fine self-convergent.
"""
import numpy as np
from math import pi
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402
from radia.hdiv_vim import (  # noqa: E402
    reconstruct_field_polynomial,
    reconstruct_field_internal,
    flat_triangle_charge_field,
    tet_self_volume_field,
    triangle_potential_const,
    triangle_potential_moment,
    tet_newtonian_potential,
    tet_volume_field_linear,
    linear_triangle_charge_field,
    triangle_potential_moment2,
    tet_newtonian_moment,
    tet_volume_field_quadratic,
    quadratic_triangle_charge_field,
    polynomial_triangle_charge_field,
    tet_volume_field_polynomial,
)


def _sphere(h):
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=h))


def _H_dipole(r, m):
    rn = np.linalg.norm(r)
    rh = r / rn
    mv = np.array([0.0, 0.0, m])
    return (1.0 / (4 * pi)) * (3.0 * np.dot(mv, rh) * rh - mv) / rn ** 3


@pytest.fixture(scope="module")
def uniform_sphere():
    mesh = _sphere(0.4)
    Mval = 5.0e5
    fes = ng.HDiv(mesh, order=1)
    gf = ng.GridFunction(fes)
    with ng.TaskManager():
        gf.Set(ng.CoefficientFunction((0, 0, Mval)))
    return mesh, gf, Mval


def test_uniform_center_is_minus_M_over_3(uniform_sphere):
    """div M = 0 (uniform): the surface charge alone gives the exact -M/3 at the sphere center."""
    mesh, gf, Mval = uniform_sphere
    with ng.TaskManager():
        Hc = reconstruct_field_polynomial(mesh, gf, np.array([[0.0, 0.0, 0.0]]), quad=4)
    rel = abs(Hc[0, 2] + Mval / 3) / (Mval / 3)
    assert rel < 5e-3, f"center H_z {Hc[0,2]:.4e} vs -M/3 {-Mval/3:.4e} (rel {rel:.2e})"
    assert abs(Hc[0, 0]) < 1e-2 * Mval and abs(Hc[0, 1]) < 1e-2 * Mval, "transverse leak at center"


def test_uniform_external_is_dipole(uniform_sphere):
    """uniform M: external field == the analytic dipole (within the flat-mesh faceting ~6% at h=0.4)."""
    mesh, gf, Mval = uniform_sphere
    m_dip = Mval * (4 * pi / 3)
    obs = np.array([[0, 0, 2.0], [0, 0, 3.0], [2.0, 0, 0.0]], float)
    with ng.TaskManager():
        H = reconstruct_field_polynomial(mesh, gf, obs, quad=4)
    for i, r in enumerate(obs):
        Hd = _H_dipole(r, m_dip)
        rel = np.linalg.norm(H[i] - Hd) / np.linalg.norm(Hd)
        assert rel < 0.08, f"external r={r.tolist()} rel {rel:.3e} (faceting expected <8% at h=0.4)"


def test_volume_charge_is_essential_for_div_M():
    """linear M (div M = M0 != 0): the full kernel is coarse->fine convergent, but dropping rho
    (surface-only = the old centroid reconstruct_field) is grossly wrong (>50%)."""
    M0 = 3.0e5
    obs = np.array([[0, 0, 2.5], [2.5, 0, 0.0], [1.5, 1.5, 1.0]], float)

    def fld(h, drop_rho=False):
        me = _sphere(h)
        fe = ng.HDiv(me, order=1)
        g = ng.GridFunction(fe)
        with ng.TaskManager():
            g.Set(ng.CoefficientFunction((0, 0, M0 * (1 + ng.z))))
            return reconstruct_field_polynomial(me, g, obs, quad=4, include_volume=not drop_rho)

    Hcoarse = fld(0.5)
    Hfine = fld(0.3)
    Hsig = fld(0.5, drop_rho=True)
    for i in range(len(obs)):
        rel_cf = np.linalg.norm(Hcoarse[i] - Hfine[i]) / np.linalg.norm(Hfine[i])
        rel_sig = np.linalg.norm(Hsig[i] - Hfine[i]) / np.linalg.norm(Hfine[i])
        assert rel_cf < 0.12, f"full kernel coarse->fine not convergent at obs {i}: {rel_cf:.2e}"
        assert rel_sig > 0.4, f"dropping rho should be grossly wrong at obs {i}, got {rel_sig:.2e}"


# --------------------------------------------------------------------------------------------------
# HEXAHEDRAL meshes: the kernel is element-type agnostic (NGSolve GetTrafo + IntegrationRule(el.type)
# + specialcf.normal), so a hex mesh of the SAME body must give the SAME external field as the tet
# mesh -- the physics is identical, only the interior meshing differs.  This locks the hex path.
# --------------------------------------------------------------------------------------------------
def _box_field(hexes, M_cf, obs, n=3, order=1):
    me = MakeStructured3DMesh(hexes=hexes, nx=n, ny=n, nz=n, mapping=lambda x, y, z: (x, y, z))
    fe = ng.HDiv(me, order=order)
    g = ng.GridFunction(fe)
    with ng.TaskManager():
        g.Set(M_cf)
        return reconstruct_field_polynomial(me, g, obs, quad=4)


def test_hex_matches_tet_uniform_M():
    """uniform M (div M = 0, surface charge only): a HEX box and a TET box give the same external
    field to ~machine precision (the boundary surface charge M.n is identical; only the interior
    meshing differs)."""
    Mval = 5.0e5
    obs = np.array([[0.5, 0.5, 3.0], [3.0, 0.5, 0.5], [2.0, 2.0, 2.0]], float)
    Mcf = ng.CoefficientFunction((0, 0, Mval))
    Hhex = _box_field(True, Mcf, obs)
    Htet = _box_field(False, Mcf, obs)
    for i in range(len(obs)):
        rel = np.linalg.norm(Hhex[i] - Htet[i]) / np.linalg.norm(Htet[i])
        assert rel < 1e-9, f"hex vs tet uniform-M external field differ at obs {i}: {rel:.2e}"


def test_hex_matches_tet_linear_M():
    """linear M (div M = M0 != 0, the VOLUME charge path): a HEX box and a TET box of the same body
    with the same (exactly representable) linear M give the same external field -- locks the hex
    volume-charge quadrature, not just the surface charge."""
    M0 = 3.0e5
    obs = np.array([[0.5, 0.5, 3.0], [3.0, 0.5, 0.5], [2.0, 2.0, 2.0]], float)
    Mcf = ng.CoefficientFunction((0, 0, M0 * (1 + ng.z)))
    Hhex = _box_field(True, Mcf, obs)
    Htet = _box_field(False, Mcf, obs)
    for i in range(len(obs)):
        rel = np.linalg.norm(Hhex[i] - Htet[i]) / np.linalg.norm(Htet[i])
        assert rel < 1e-6, f"hex vs tet linear-M (div M != 0) external field differ at obs {i}: {rel:.2e}"


# --------------------------------------------------------------------------------------------------
# Step 2 building blocks (the singular / near-singular INTERNAL field), each validated standalone.
# --------------------------------------------------------------------------------------------------
def _tri_field_gauss(P, r, ng=40):
    """fine Duffy-Gauss reference for INT_T (r-r')/|r-r'|^3 dS' (non-singular at a moderate r)."""
    P = np.asarray(P, float)
    x, w = np.polynomial.legendre.leggauss(ng)
    s = 0.5 * (x + 1); ws = 0.5 * w
    e1, e2 = P[1] - P[0], P[2] - P[0]
    area2 = np.linalg.norm(np.cross(e1, e2))
    F = np.zeros(3)
    for u, wu in zip(s, ws):
        for v, wv in zip(s, ws):
            rp = P[0] + u * e1 + (v * (1 - u)) * e2
            d = r - rp
            F += (wu * wv * (1 - u) * area2) * d / np.linalg.norm(d) ** 3
    return F


def test_flat_triangle_charge_field_exact():
    """The analytic uniform-triangle field matches a fine Gauss reference to ~machine precision."""
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    for r in [np.array([0.3, 0.3, 0.5]), np.array([0.3, 0.3, -0.5]), np.array([1.0, 1.0, 0.3])]:
        Fa = flat_triangle_charge_field(P, r)
        Fg = _tri_field_gauss(P, r)
        rel = np.linalg.norm(Fa - Fg) / np.linalg.norm(Fg)
        assert rel < 1e-7, f"analytic triangle field vs Gauss r={r.tolist()}: rel {rel:.2e}"


def test_tet_self_volume_field_vs_phitet_gradient():
    """The spherical ray-trace self volume-charge field (constant rho) matches -(rho/4pi) grad(phi_tet),
    the gradient of the EXACT analytic Newtonian potential (radia _hdiv_phi_tet) -- the 1/r^2 singularity
    is removed analytically, so it is accurate at INTERIOR points."""
    import radia._radia_pybind as rp
    V = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]   # unit tet
    Varr = np.array(V, float).reshape(4, 3)
    rho0 = 1.0
    delta = 1e-5
    for r in [np.array([0.25, 0.25, 0.25]), np.array([0.15, 0.30, 0.20]), np.array([0.4, 0.3, 0.1])]:
        Hs = tet_self_volume_field(Varr, r, lambda p: rho0, nth=48, nph=96)
        g = np.zeros(3)
        for k in range(3):
            rp_, rm_ = list(r), list(r)
            rp_[k] += delta; rm_[k] -= delta
            g[k] = (rp._hdiv_phi_tet(V, [float(x) for x in rp_])
                    - rp._hdiv_phi_tet(V, [float(x) for x in rm_])) / (2 * delta)
        Href = -(rho0 / (4 * np.pi)) * g
        rel = np.linalg.norm(Hs - Href) / np.linalg.norm(Href)
        assert rel < 8e-3, f"spherical self volume field vs grad(phi_tet) r={r.tolist()}: rel {rel:.2e}"


def test_internal_field_assembly_uniform_sphere():
    """The assembled INTERNAL field (self-volume spherical + far-volume + analytic surface) on a uniform
    sphere: center -> -M/3 (validates the assembly factors); and near-surface the analytic surface beats
    the Step-1 plain Gauss-Duffy by a large margin (Step-1 is near-singular there)."""
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = ng.Mesh(g.GenerateMesh(maxh=0.4))
    Mval = 5.0e5
    fes = ng.HDiv(mesh, order=1)
    gf = ng.GridFunction(fes)
    with ng.TaskManager():
        gf.Set(ng.CoefficientFunction((0, 0, Mval)))
    d = np.array([0.3, 0.2, 0.9]); d = d / np.linalg.norm(d)
    pts = np.array([[0, 0, 0.0], (0.95 * d).tolist()])
    with ng.TaskManager():
        Hint = reconstruct_field_internal(mesh, gf, pts)
        Hstep1 = reconstruct_field_polynomial(mesh, gf, pts, quad=4)
    # center: -M/3
    rel_c = abs(Hint[0, 2] + Mval / 3) / (Mval / 3)
    assert rel_c < 5e-3, f"assembled center H_z {Hint[0,2]:.4e} vs -M/3: rel {rel_c:.2e}"
    # near surface (0.95R): assembled is far better than Step-1 (which is ~58% off there)
    rel_int = abs(Hint[1, 2] + Mval / 3) / (Mval / 3)
    rel_s1 = abs(Hstep1[1, 2] + Mval / 3) / (Mval / 3)
    assert rel_int < 0.05, f"assembled near-surface relZ {rel_int:.2e} not < 5%"
    assert rel_int < 0.3 * rel_s1, f"assembled ({rel_int:.2e}) should beat Step-1 ({rel_s1:.2e}) >3x near surface"


# --------------------------------------------------------------------------------------------------
# C++ field kernels (the order>=2 field accelerated in C++) vs the Python goldens -- the lab pattern
# (C++ probe validated entry-by-entry against the Python reference, like _hdiv_phi_tet vs phi_tet).
# --------------------------------------------------------------------------------------------------
def test_cpp_tri_field_matches_python():
    """C++ _hdiv_tri_field (Wilton triangle field) == Python flat_triangle_charge_field to ~machine
    precision (same closed form)."""
    import radia._radia_pybind as rp
    P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    Vflat = P.ravel().tolist()
    for r in [[0.3, 0.3, 0.5], [0.3, 0.3, -0.5], [1.0, 1.0, 0.3], [0.2, 0.2, 0.05]]:
        Fc = np.array(rp._hdiv_tri_field(Vflat, r))
        Fp = flat_triangle_charge_field(P, np.array(r))
        rel = np.linalg.norm(Fc - Fp) / np.linalg.norm(Fp)
        assert rel < 1e-12, f"C++ tri_field vs Python r={r}: rel {rel:.2e}"


def test_cpp_tet_field_matches_grad_phitet():
    """C++ _hdiv_tet_field (= -grad PhiTet, the analytic tet volume-charge field) matches a central FD
    of the analytic _hdiv_phi_tet to ~machine precision (valid near AND far), and matches the Python
    spherical ray-trace (tet_self_volume_field*4pi) to the spherical method's ~1e-3 at interior points."""
    import radia._radia_pybind as rp
    Vtet = [0.0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]
    Varr = np.array(Vtet, float).reshape(4, 3)
    d = 1e-5
    for r in [[0.25, 0.25, 0.25], [0.4, 0.3, 0.1], [2.0, 0.0, 0.0]]:
        Fc = np.array(rp._hdiv_tet_field(Vtet, r))
        g = np.zeros(3)
        for k in range(3):
            rpp, rmm = list(r), list(r)
            rpp[k] += d; rmm[k] -= d
            g[k] = (rp._hdiv_phi_tet(Vtet, rpp) - rp._hdiv_phi_tet(Vtet, rmm)) / (2 * d)
        rel_fd = np.linalg.norm(Fc + g) / np.linalg.norm(g)         # Fc == -grad(phi_tet)
        assert rel_fd < 1e-6, f"C++ tet_field vs -grad(phi_tet) FD r={r}: rel {rel_fd:.2e}"
        if max(r) < 1 and min(r) >= 0 and sum(r) < 1:               # interior: vs the spherical golden
            Fsph = 4 * np.pi * tet_self_volume_field(Varr, np.array(r), lambda p: 1.0, nth=48, nph=96)
            rel_sph = np.linalg.norm(Fc - Fsph) / np.linalg.norm(Fsph)
            assert rel_sph < 8e-3, f"C++ tet_field vs spherical r={r}: rel {rel_sph:.2e}"


# ---------------------------------------------------------------------------------------------------
# Analytic LINEAR volume-charge field (tet_volume_field_linear): EXACT closed form, the order-2 -div M
# term.  Validated 4 ways -- pure-Python building blocks vs C++ probes (machine), constant case vs the
# C++ TetField (machine), linear case vs far Gauss (machine), interior linear vs spherical (~1e-3).
# ---------------------------------------------------------------------------------------------------
def _tri_pot_gauss(P, r, nq=40):
    P = np.asarray(P, float); r = np.asarray(r, float)
    J = np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    I0 = 0.0; M1 = np.zeros(3)
    for i in range(nq):
        for j in range(nq):
            u, w = xs[i], xs[j]; jac = (1 - u)
            p = P[0] + u * (P[1] - P[0]) + w * (1 - u) * (P[2] - P[0])
            wgt = ws[i] * ws[j] * jac * J / np.linalg.norm(r - p)
            I0 += wgt; M1 += wgt * p
    return I0, M1


def _tet_field_gauss(tet, r, rho0, g, nq=24):
    t = np.asarray(tet, float); r = np.asarray(r, float)
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    V6 = abs(np.linalg.det(np.array([t[1] - t[0], t[2] - t[0], t[3] - t[0]])))
    acc = np.zeros(3)
    for i in range(nq):
        for j in range(nq):
            for k in range(nq):
                a, b, c = xs[i], xs[j], xs[k]
                l1, l2, l3 = a, b * (1 - a), c * (1 - a) * (1 - b)
                jac = (1 - a) ** 2 * (1 - b)
                p = t[0] + l1 * (t[1] - t[0]) + l2 * (t[2] - t[0]) + l3 * (t[3] - t[0])
                acc += ws[i] * ws[j] * ws[k] * jac * V6 * (rho0 + np.dot(g, p)) * (r - p) / np.linalg.norm(r - p) ** 3
    return acc


_TET = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)


def test_triangle_potential_pure_python_vs_cpp_and_gauss():
    """Pure-Python triangle potential (const I0 + first moment M1) == C++ _hdiv_tri_potential (machine)
    and == off-plane Gauss (machine)."""
    import radia._radia_pybind as rp
    P = _TET[[0, 1, 2]]
    for r in [[0.3, 0.3, 0.7], [0.3, 0.3, 0.2], [1.5, 0.5, 0.4]]:
        r = np.array(r, float)
        I0 = triangle_potential_const(P, r)
        M1 = triangle_potential_moment(P, r)
        I0c = rp._hdiv_tri_potential(P.ravel().tolist(), r.tolist())
        I0g, M1g = _tri_pot_gauss(P, r, 40)
        assert abs(I0 - I0c) < 1e-12, f"I0 vs C++ r={r}"
        assert abs(I0 - I0g) / abs(I0g) < 1e-10, f"I0 vs Gauss r={r}"
        assert np.linalg.norm(M1 - M1g) / np.linalg.norm(M1g) < 1e-10, f"M1 vs Gauss r={r}"


def test_tet_newtonian_potential_vs_cpp():
    """Pure-Python PhiTet (-1/2 SUM h_f I0_f) == C++ _hdiv_phi_tet (machine), interior AND exterior."""
    import radia._radia_pybind as rp
    for r in [[0.25, 0.25, 0.25], [0.4, 0.3, 0.1], [2.0, 0.0, 0.0], [0.1, 0.1, 0.1]]:
        r = np.array(r, float)
        a = tet_newtonian_potential(_TET, r)
        b = rp._hdiv_phi_tet(_TET.ravel().tolist(), r.tolist())
        assert abs(a - b) / abs(b) < 1e-12, f"PhiTet vs C++ r={r}: {abs(a - b) / abs(b):.2e}"


def test_tet_volume_field_linear_constant_vs_cpp_tetfield():
    """Constant-rho closed-form volume field (SUM_f n_f I0_f) == the C++ TetField (= -grad PhiTet),
    two INDEPENDENT derivations, to machine precision."""
    import radia._radia_pybind as rp
    for r in [[0.25, 0.25, 0.25], [0.4, 0.3, 0.1], [2.0, 0.0, 0.0], [0.3, 0.3, 0.7]]:
        r = np.array(r, float)
        Fc = tet_volume_field_linear(_TET, r, 1.0, np.zeros(3))
        Ft = np.array(rp._hdiv_tet_field(_TET.ravel().tolist(), r.tolist()))
        rel = np.linalg.norm(Fc - Ft) / np.linalg.norm(Ft)
        assert rel < 1e-12, f"const volume field vs C++ TetField r={r}: {rel:.2e}"


def test_tet_volume_field_linear_vs_gauss_far():
    """Linear-rho closed-form volume field == high-order tet Gauss at FAR (non-singular) points."""
    g = np.array([0.7, -0.3, 0.5]); rho0 = 0.4
    for r in [[2.0, 0.0, 0.0], [1.5, 1.5, 1.0], [-1.0, 0.5, 0.5], [0.5, 0.5, 2.5]]:
        r = np.array(r, float)
        Fc = tet_volume_field_linear(_TET, r, rho0, g)
        Fg = _tet_field_gauss(_TET, r, rho0, g, 24)
        rel = np.linalg.norm(Fc - Fg) / np.linalg.norm(Fg)
        assert rel < 1e-10, f"linear volume field vs far Gauss r={r}: {rel:.2e}"


def test_tet_volume_field_linear_interior_vs_spherical():
    """Linear-rho closed form == the spherical ray-trace (tet_self_volume_field, which integrates a
    polynomial rho) at INTERIOR points, to the spherical method's own ~1e-3 accuracy -- confirms the
    closed form IS the exact value the ~1e-3 spherical method converges to."""
    g = np.array([0.7, -0.3, 0.5]); rho0 = 0.4
    for r in [[0.25, 0.25, 0.25], [0.15, 0.30, 0.20], [0.1, 0.1, 0.1]]:
        r = np.array(r, float)
        Fc = tet_volume_field_linear(_TET, r, rho0, g)
        Fs = 4 * np.pi * tet_self_volume_field(_TET, r, lambda p: rho0 + float(np.dot(g, p)),
                                               nth=64, nph=128, ns=12)
        rel = np.linalg.norm(Fc - Fs) / np.linalg.norm(Fs)
        assert rel < 2e-3, f"linear volume field interior vs spherical r={r}: {rel:.2e}"


# ---------------------------------------------------------------------------------------------------
# Analytic LINEAR surface-charge field (linear_triangle_charge_field): EXACT closed form, the order-2
# M.n term.  H_surf = -grad_r phi_sigma with phi_sigma = sigma0 I0 + s.M1 (the degree-1 triangle
# potential).  Validated vs off-plane Gauss (machine); the s=0 case reproduces the constant field.
# ---------------------------------------------------------------------------------------------------
def _lin_tri_field_gauss(P, r, sigma0, s, nq=44):
    P = np.asarray(P, float); r = np.asarray(r, float); s = np.asarray(s, float)
    J = np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    F = np.zeros(3)
    for i in range(nq):
        for j in range(nq):
            u, v = xs[i], xs[j]; jac = (1 - u)
            p = P[0] + u * (P[1] - P[0]) + v * (1 - u) * (P[2] - P[0])
            d = r - p
            F += ws[i] * ws[j] * jac * J * (sigma0 + np.dot(s, p)) * d / np.linalg.norm(d) ** 3
    return F


def test_linear_triangle_charge_field_vs_gauss():
    """Linear-sigma triangle field == off-plane Gauss to machine precision (near AND far)."""
    P = _TET[[0, 1, 2]]
    s = np.array([0.6, -0.4, 0.3]); sigma0 = 0.5
    for r in [[0.3, 0.3, 0.7], [0.3, 0.3, 0.2], [1.5, 0.5, 0.4], [0.2, 0.2, 0.3], [-0.5, 0.4, 0.6]]:
        r = np.array(r, float)
        Fc = linear_triangle_charge_field(P, r, sigma0, s)
        Fg = _lin_tri_field_gauss(P, r, sigma0, s, 44)
        rel = np.linalg.norm(Fc - Fg) / np.linalg.norm(Fg)
        assert rel < 1e-10, f"linear-sigma tri field vs Gauss r={r}: {rel:.2e}"


def test_linear_triangle_charge_field_constant_reduces():
    """s = 0 reproduces sigma0 * flat_triangle_charge_field exactly (bit-identical)."""
    P = _TET[[0, 1, 2]]
    sigma0 = 0.5
    for r in [[0.3, 0.3, 0.7], [1.5, 0.5, 0.4], [0.2, 0.2, 0.3]]:
        r = np.array(r, float)
        a = linear_triangle_charge_field(P, r, sigma0, np.zeros(3))
        b = sigma0 * flat_triangle_charge_field(P, r)
        assert np.linalg.norm(a - b) < 1e-13, f"s=0 vs sigma0*const r={r}"


# ---------------------------------------------------------------------------------------------------
# QUADRATIC volume-charge field (tet_volume_field_quadratic) + its degree-2 moment building blocks:
# surface second moment M2 = INT_T r'(x)r'/R dS', volume first moment V1 = INT_V r'/R dV'.
# ---------------------------------------------------------------------------------------------------
def _tri_moment2_gauss(P, r, nq=44):
    P = np.asarray(P, float); r = np.asarray(r, float)
    J = np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    M = np.zeros((3, 3))
    for i in range(nq):
        for j in range(nq):
            u, v = xs[i], xs[j]; jac = (1 - u)
            p = P[0] + u * (P[1] - P[0]) + v * (1 - u) * (P[2] - P[0])
            M += ws[i] * ws[j] * jac * J * np.outer(p, p) / np.linalg.norm(r - p)
    return M


def _tet_field_quad_gauss(V, r, rho0, g, Q, nq=26):
    V = np.asarray(V, float); r = np.asarray(r, float)
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    V6 = abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]])))
    acc = np.zeros(3)
    for i in range(nq):
        for j in range(nq):
            for k in range(nq):
                a, b, c = xs[i], xs[j], xs[k]
                l1, l2, l3 = a, b * (1 - a), c * (1 - a) * (1 - b); jac = (1 - a) ** 2 * (1 - b)
                p = V[0] + l1 * (V[1] - V[0]) + l2 * (V[2] - V[0]) + l3 * (V[3] - V[0])
                rho = rho0 + np.dot(g, p) + p @ Q @ p
                acc += ws[i] * ws[j] * ws[k] * jac * V6 * rho * (r - p) / np.linalg.norm(r - p) ** 3
    return acc


_QSYM = np.array([[0.3, 0.1, -0.2], [0.1, -0.4, 0.15], [-0.2, 0.15, 0.25]])


def test_triangle_potential_moment2_vs_gauss_symmetric():
    """Surface second moment M2 == off-plane Gauss (machine) and is symmetric."""
    P = _TET[[0, 1, 2]]
    for r in [[0.3, 0.3, 0.7], [0.3, 0.3, 0.2], [1.5, 0.5, 0.4], [0.2, 0.2, 0.3]]:
        r = np.array(r, float)
        a = triangle_potential_moment2(P, r)
        b = _tri_moment2_gauss(P, r, 44)
        assert np.linalg.norm(a - b) / np.linalg.norm(b) < 1e-10, f"M2 vs Gauss r={r}"
        assert np.linalg.norm(a - a.T) / np.linalg.norm(a) < 1e-12, f"M2 not symmetric r={r}"


def test_tet_newtonian_moment_vs_gauss_far():
    """Volume first moment V1 == far tet Gauss (machine)."""
    def gV1(r, nq=24):
        xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
        V6 = abs(np.linalg.det(np.array([_TET[1] - _TET[0], _TET[2] - _TET[0], _TET[3] - _TET[0]])))
        acc = np.zeros(3)
        for i in range(nq):
            for j in range(nq):
                for k in range(nq):
                    a, b, c = xs[i], xs[j], xs[k]
                    l1, l2, l3 = a, b * (1 - a), c * (1 - a) * (1 - b); jac = (1 - a) ** 2 * (1 - b)
                    p = _TET[0] + l1 * (_TET[1] - _TET[0]) + l2 * (_TET[2] - _TET[0]) + l3 * (_TET[3] - _TET[0])
                    acc += ws[i] * ws[j] * ws[k] * jac * V6 * p / np.linalg.norm(r - p)
        return acc
    for r in [[2.0, 0, 0], [1.5, 1.5, 1.0], [-1.0, 0.5, 0.5]]:
        r = np.array(r, float)
        a = tet_newtonian_moment(_TET, r)
        b = gV1(r, 24)
        assert np.linalg.norm(a - b) / np.linalg.norm(b) < 1e-10, f"V1 vs Gauss r={r}"


def test_tet_volume_field_quadratic_vs_gauss_far():
    """Quadratic-rho closed-form volume field == far tet Gauss (machine)."""
    g = np.array([0.7, -0.3, 0.5]); rho0 = 0.4
    for r in [[2.0, 0, 0], [1.5, 1.5, 1.0], [-1.0, 0.5, 0.5], [0.5, 0.5, 2.5]]:
        r = np.array(r, float)
        Fc = tet_volume_field_quadratic(_TET, r, rho0, g, _QSYM)
        Fg = _tet_field_quad_gauss(_TET, r, rho0, g, _QSYM, 26)
        assert np.linalg.norm(Fc - Fg) / np.linalg.norm(Fg) < 1e-9, f"quad volume field vs Gauss r={r}"


def test_tet_volume_field_quadratic_reduces_to_linear():
    """Q = 0 reproduces tet_volume_field_linear bit-identically."""
    g = np.array([0.7, -0.3, 0.5]); rho0 = 0.4
    for r in [[2.0, 0, 0], [0.25, 0.25, 0.25], [0.1, 0.1, 0.1]]:
        r = np.array(r, float)
        a = tet_volume_field_quadratic(_TET, r, rho0, g, np.zeros((3, 3)))
        b = tet_volume_field_linear(_TET, r, rho0, g)
        assert np.linalg.norm(a - b) < 1e-13, f"Q=0 vs linear r={r}"


# ---------------------------------------------------------------------------------------------------
# QUADRATIC surface-charge field (quadratic_triangle_charge_field): EXACT closed form via the
# in-plane/normal split + 1/R^3 moments.  The degree-2 surface companion of tet_volume_field_quadratic.
# ---------------------------------------------------------------------------------------------------
def _quad_tri_field_gauss(P, r, sigma0, s, S, nq=46):
    P = np.asarray(P, float); r = np.asarray(r, float); s = np.asarray(s, float); S = np.asarray(S, float)
    J = np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    F = np.zeros(3)
    for i in range(nq):
        for j in range(nq):
            u, v = xs[i], xs[j]; jac = (1 - u)
            p = P[0] + u * (P[1] - P[0]) + v * (1 - u) * (P[2] - P[0])
            d = r - p
            F += ws[i] * ws[j] * jac * J * (sigma0 + np.dot(s, p) + p @ S @ p) * d / np.linalg.norm(d) ** 3
    return F


def test_quadratic_triangle_charge_field_vs_gauss():
    """Quadratic-sigma triangle field == off-plane Gauss to machine precision (near AND far)."""
    P = _TET[[0, 1, 2]]
    s = np.array([0.6, -0.4, 0.3]); sigma0 = 0.5
    for r in [[0.3, 0.3, 0.7], [0.3, 0.3, 0.2], [1.5, 0.5, 0.4], [0.2, 0.2, 0.3], [-0.5, 0.4, 0.6]]:
        r = np.array(r, float)
        Fc = quadratic_triangle_charge_field(P, r, sigma0, s, _QSYM)
        Fg = _quad_tri_field_gauss(P, r, sigma0, s, _QSYM, 46)
        assert np.linalg.norm(Fc - Fg) / np.linalg.norm(Fg) < 1e-10, f"quad-sigma tri field r={r}"


def test_quadratic_triangle_charge_field_reduces_to_linear():
    """S = 0 matches linear_triangle_charge_field (two independent derivations -- in-plane/normal vs
    -grad phi -- agree to machine precision)."""
    P = _TET[[0, 1, 2]]
    s = np.array([0.6, -0.4, 0.3]); sigma0 = 0.5
    for r in [[0.3, 0.3, 0.7], [1.5, 0.5, 0.4], [0.2, 0.2, 0.3]]:
        r = np.array(r, float)
        a = quadratic_triangle_charge_field(P, r, sigma0, s, np.zeros((3, 3)))
        b = linear_triangle_charge_field(P, r, sigma0, s)
        assert np.linalg.norm(a - b) / np.linalg.norm(b) < 1e-11, f"S=0 vs linear r={r}"


# ---------------------------------------------------------------------------------------------------
# ARBITRARY-DEGREE general assembler (polynomial_triangle_charge_field / tet_volume_field_polynomial):
# the general moment recursion.  Validated vs Gauss for cubic charge, and consistent with the degree-2
# closed forms (two independent code paths agree).
# ---------------------------------------------------------------------------------------------------
def _rand_poly(seed, degree):
    rng = np.random.RandomState(seed)
    cf = {(ax, ay, k - ax - ay): rng.randn()
          for k in range(degree + 1) for ax in range(k + 1) for ay in range(k + 1 - ax)}
    return lambda p: sum(c * p[0] ** ax * p[1] ** ay * p[2] ** az for (ax, ay, az), c in cf.items())


def test_polynomial_triangle_charge_field_cubic_vs_gauss():
    """General surface assembler: cubic sigma == off-plane Gauss to machine precision (near AND far)."""
    P = _TET[[0, 1, 2]]
    sig = _rand_poly(3, 3)
    for r in [[0.3, 0.3, 0.7], [0.3, 0.3, 0.2], [1.5, 0.5, 0.4], [-0.5, 0.4, 0.6]]:
        r = np.array(r, float)
        Fc = polynomial_triangle_charge_field(P, r, sig, 3)
        Fg = _quad_tri_field_gauss_fn(P, r, sig, 52)
        assert np.linalg.norm(Fc - Fg) / np.linalg.norm(Fg) < 1e-10, f"cubic surface r={r}"


def test_polynomial_triangle_charge_field_matches_quadratic():
    """General assembler at degree 2 == the closed-form quadratic_triangle_charge_field (two paths)."""
    P = _TET[[0, 1, 2]]
    s = np.array([0.6, -0.4, 0.3]); sigma0 = 0.5
    sig = lambda p: sigma0 + np.dot(s, p) + p @ _QSYM @ p
    for r in [[0.3, 0.3, 0.7], [0.2, 0.2, 0.3], [1.5, 0.5, 0.4]]:
        r = np.array(r, float)
        a = polynomial_triangle_charge_field(P, r, sig, 2)
        b = quadratic_triangle_charge_field(P, r, sigma0, s, _QSYM)
        assert np.linalg.norm(a - b) / np.linalg.norm(b) < 1e-11, f"deg2 general vs closed r={r}"


def test_tet_volume_field_polynomial_cubic_vs_gauss():
    """General volume assembler: cubic rho == far tet Gauss to machine precision."""
    rho = _rand_poly(4, 3)
    for r in [[2.0, 0, 0], [1.5, 1.5, 1.0], [-1.0, 0.5, 0.5], [0.5, 0.5, 2.5]]:
        r = np.array(r, float)
        Fc = tet_volume_field_polynomial(_TET, r, rho, 3)
        Fg = _vol_field_gauss_fn(_TET, r, rho, 24)
        assert np.linalg.norm(Fc - Fg) / np.linalg.norm(Fg) < 1e-9, f"cubic volume r={r}"


def test_tet_volume_field_polynomial_matches_quadratic():
    """General volume assembler at degree 2 == the closed-form tet_volume_field_quadratic."""
    g = np.array([0.7, -0.3, 0.5]); rho0 = 0.4
    rho = lambda p: rho0 + np.dot(g, p) + p @ _QSYM @ p
    for r in [[2.0, 0, 0], [0.25, 0.25, 0.25], [-1.0, 0.5, 0.5]]:
        r = np.array(r, float)
        a = tet_volume_field_polynomial(_TET, r, rho, 2)
        b = tet_volume_field_quadratic(_TET, r, rho0, g, _QSYM)
        assert np.linalg.norm(a - b) / np.linalg.norm(b) < 1e-10, f"deg2 general vs closed r={r}"


def _quad_tri_field_gauss_fn(P, r, fn, nq=50):
    P = np.asarray(P, float); r = np.asarray(r, float)
    J = np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    F = np.zeros(3)
    for i in range(nq):
        for j in range(nq):
            u, v = xs[i], xs[j]; jac = (1 - u)
            p = P[0] + u * (P[1] - P[0]) + v * (1 - u) * (P[2] - P[0]); d = r - p
            F += ws[i] * ws[j] * jac * J * fn(p) * d / np.linalg.norm(d) ** 3
    return F


def _vol_field_gauss_fn(V, r, fn, nq=24):
    V = np.asarray(V, float); r = np.asarray(r, float)
    xs, ws = np.polynomial.legendre.leggauss(nq); xs = 0.5 * (xs + 1); ws = 0.5 * ws
    V6 = abs(np.linalg.det(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]])))
    F = np.zeros(3)
    for i in range(nq):
        for j in range(nq):
            for k in range(nq):
                a, b, cc = xs[i], xs[j], xs[k]
                l1, l2, l3 = a, b * (1 - a), cc * (1 - a) * (1 - b); jac = (1 - a) ** 2 * (1 - b)
                p = V[0] + l1 * (V[1] - V[0]) + l2 * (V[2] - V[0]) + l3 * (V[3] - V[0]); d = r - p
                F += ws[i] * ws[j] * ws[k] * jac * V6 * fn(p) * d / np.linalg.norm(d) ** 3
    return F
