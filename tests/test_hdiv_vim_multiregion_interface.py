"""Multi-region conforming meshes: internal interface faces carry no charge.

Locks the 2026-07-28 finding: the charge layer used to treat EVERY surface
element as an exterior boundary face (single-sided sigma = M.n).  On a
conforming multi-region mesh (OCC Glue of two materials with a shared internal
face) the internal surface elements then formed a spurious charge sheet: the
uniform chi=100 sanity on a two-region glued ball read <Mz> = 2.559 vs the
exact sphere demag 2.9126 (-12 %), with the interface flux clamped ~23x.
`_exterior_bnd_elements` now drops facets owned by two volume elements from the
charge layer (their true single-layer charge is the jump = 0 for the
HDiv-continuous normal trace), which restores the exact single-region result
and makes conforming multi-region meshes valid DemagOperator inputs.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, Glue, OCCGeometry, Pnt, Sphere  # noqa: E402
from ngsolve import (  # noqa: E402
    BilinearForm, CoefficientFunction, GridFunction, HDiv, Integrate,
    LinearForm, Mesh, TaskManager, dx,
)
from ngsolve.krylovspace import CGSolver  # noqa: E402

from radia.vim import DemagOperator  # noqa: E402
from radia.vim._vim import _exterior_bnd_elements  # noqa: E402

CHI = 100.0
EXACT_SPHERE = CHI / (1.0 + CHI / 3.0)      # uniform-sphere demag chi_eff
ONE = CoefficientFunction(1.0)
MAXH = 0.35


def _ball():
    return Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=MAXH))


def _glued_ball():
    """Conforming two-region ball: upper 'iron', lower 'void', named internal
    interface at z=0 (netgen creates internal surface elements for it)."""
    sph = Sphere(Pnt(0, 0, 0), 1.0)
    cut = Box(Pnt(-2, -2, 0), Pnt(2, 2, 2))
    up = sph * cut
    dn = sph - cut
    up.mat("iron")
    dn.mat("void")
    for solid in (up, dn):
        for face in solid.faces:
            if abs(face.center.z) < 1e-9:
                face.name = "interface"
    return Mesh(OCCGeometry(Glue([up, dn])).GenerateMesh(maxh=MAXH))


def _solve_m(mesh, weight):
    """(M_s + N) m = P(H0):  the standard linear VIM solve, H0 = z-hat."""
    fes = HDiv(mesh, order=1)
    demag = DemagOperator(fes, eps=1e-7)
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += weight * u * v * dx
    a.Assemble()
    rhs = LinearForm(fes)
    rhs += CoefficientFunction((0.0, 0.0, 1.0)) * v * dx
    rhs.Assemble()
    pre = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    inv = CGSolver(a.mat + demag.mat, pre=pre, tol=1e-10, maxiter=3000)
    gf = GridFunction(fes)
    gf.vec.data = inv * rhs.vec
    return gf


def test_exterior_filter_is_noop_on_single_region_mesh():
    mesh = _ball()
    kept = _exterior_bnd_elements(mesh)
    assert len(kept) == mesh.GetNE(ng.BND)


def test_internal_interface_faces_are_dropped():
    mesh = _glued_ball()
    iface_nrs = {el.nr for el in mesh.Elements(ng.BND) if el.mat == "interface"}
    assert iface_nrs, "test mesh must contain internal surface elements"
    kept = _exterior_bnd_elements(mesh)
    kept_nrs = {e.nr for e in kept}
    assert len(kept) == mesh.GetNE(ng.BND) - len(iface_nrs)
    assert not (kept_nrs & iface_nrs)


def test_uniform_chi_multiregion_ball_matches_sphere_demag():
    """Pre-fix this read ~2.559 (-12 %); the measured post-fix value is 2.9190."""
    mesh = _glued_ball()
    with TaskManager():
        gf = _solve_m(mesh, CoefficientFunction(1.0 / CHI))
        mz = Integrate(gf[2], mesh) / Integrate(ONE, mesh)
    assert abs(mz - EXACT_SPHERE) / EXACT_SPHERE < 6e-3


def test_per_region_chi_void_response_is_physical():
    """Iron/void split: the void mean Mz is the physical chi_v * H_local
    (measured 1.32e-3 at chi_v = 1e-3, local H ~ 1.3), not a charge-sheet or
    leak artefact; the iron half stays in its measured band (1.9797)."""
    chi_void = 1e-3
    mesh = _glued_ball()
    with TaskManager():
        weight = mesh.MaterialCF({"iron": 1.0 / CHI, "void": 1.0 / chi_void})
        gf = _solve_m(mesh, weight)
        v_iron = Integrate(ONE, mesh, definedon=mesh.Materials("iron"))
        v_void = Integrate(ONE, mesh, definedon=mesh.Materials("void"))
        m_iron = Integrate(gf[2], mesh, definedon=mesh.Materials("iron")) / v_iron
        m_void = Integrate(gf[2], mesh, definedon=mesh.Materials("void")) / v_void
    assert 1.85 < m_iron < 2.10
    assert 0.8e-3 < m_void < 1.8e-3
