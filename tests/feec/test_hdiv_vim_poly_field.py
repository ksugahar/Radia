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
from radia.hdiv_vim import reconstruct_field_polynomial  # noqa: E402


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
