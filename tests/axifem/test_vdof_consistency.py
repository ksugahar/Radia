# -*- coding: utf-8 -*-
r"""One field convention across the native matrices, the evaluators and the traces.

Guards:
  * The native stiffness and sigma-mass integrators assemble the same bilinear
    forms as NGSolve symbolic forms built from the evaluated basis u = A_phi
    (``u``, ``grad(u)``): (2 pi / mu) INT (u_z v_z + (u_r + u/r)(v_r + v/r)) r
    and 2 pi sigma INT u v r.  A basis or gradient convention that differs from
    the integrators (psi instead of A_phi, a missing 1/r, a reference-gradient
    chain rule on a map that is not NGSolve's) breaks the equality.  Covered:
    P1 and P2 triangles, closed-form Q1/Q2 and general Q2 quads, off-axis and
    axis-touching quads, and the counterexample triangle (1,0), (2,1), (3,2.1)
    whose affine (s, z) chord map folds.
  * The boundary trace (``ds`` evaluator) equals the volume value of the
    adjacent element at the same physical point, for every DOF, on horizontal,
    vertical, oblique and axis-connected boundary edges.
  * The axis limit of the basis: A_phi -> 0 linearly in r, dA/dr finite and
    continuous, dA/dz -> 0, no NaN at r = 0.
  * The evaluated gradient is the derivative of the evaluated value: ``grad``
    equals central finite differences of the value at interior points, also on
    elements with an oblique edge ending on the axis, whose target map shifts z
    (the counterexample triangle (0,0), (1,0), (1,1) gave (-1/3, 5/3) for the
    true (0, 1) before the chain rule through that shift).
Assumes: curved (Curve(2)) edges off the axis; the curved meshes carry a curved
interface shared by two elements, whose traces must agree from both sides.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
import ngsolve as ng
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")
from netgen.csg import Pnt
from netgen.geom2d import SplineGeometry
from netgen.meshing import Element2D, FaceDescriptor, Mesh as NGMesh, MeshPoint
from radia.axifem import AxiHenrotteSigmaMassBFI, AxiHenrotteStiffnessBFI, H1Henrotte

from _vol_mesh import structured_rect_vol_mesh

MU, SIGMA = 1.3, 2.1
MATRIX_RTOL = 1.0e-6     # native 8x8 rule vs NGSolve's rule on the 1/r integrand
TRACE_LIMIT = 1.0e-10
REF_VERTICES = {
    ng.ET.TRIG: [(1.0, 0.0), (0.0, 1.0), (0.0, 0.0)],
    ng.ET.QUAD: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
}


def _mesh(points, elements):
    m = NGMesh(dim=2)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    ids = [m.Add(MeshPoint(Pnt(px, py, 0))) for px, py in points]
    for e in elements:
        m.Add(Element2D(1, [ids[i] for i in e]))
    m.SetMaterial(1, "c")
    return ng.Mesh(m)


def _off_axis_trig_mesh():
    geo = SplineGeometry()
    geo.AddRectangle((1.0, -0.5), (2.0, 0.5), bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=0.35))


def _skew_quad_mesh(r0):
    return structured_rect_vol_mesh(
        r0, r0 + 1.0, 0.0, 1.0, quads=True, nx=2, ny=2,
        mapping=lambda xi, eta: (r0 + xi + 0.3 * eta * xi, eta + 0.2 * xi),
        stem="axifem_vdof_skew")


def _curved_two_domain_mesh(quads=False):
    # A Curve(2) circle inside a rectangle: the curved interface edges are
    # shared by two volume elements, the outer edges are straight.
    geo = SplineGeometry()
    geo.AddRectangle((0.6, -0.7), (1.8, 0.7), leftdomain=2, rightdomain=0, bc="outer")
    geo.AddCircle((1.2, 0.0), 0.4, leftdomain=1, rightdomain=2, bc="interface")
    geo.SetMaterial(1, "c")
    geo.SetMaterial(2, "a")
    mesh = ng.Mesh(geo.GenerateMesh(maxh=0.25, quad_dominated=quads))
    mesh.Curve(2)
    return mesh


def _curved_quad_sector_mesh():
    # Two annular sectors (r 0.8-1.1 and 1.1-1.4) meshed with Curve(2) quads;
    # Netgen's quad-dominated meshes of a circle in a box contain nearly
    # degenerate quads, these do not.
    from netgen.occ import Glue, MoveTo, OCCGeometry, WorkPlane
    box = MoveTo(0.3, -0.45).Rectangle(1.3, 0.9).Face()
    a = (WorkPlane().Circle(1.1).Face() - WorkPlane().Circle(0.8).Face()) * box
    b = (WorkPlane().Circle(1.4).Face() - WorkPlane().Circle(1.1).Face()) * box
    a.faces.name = "c"
    b.faces.name = "a"
    mesh = ng.Mesh(OCCGeometry(Glue([a, b]), dim=2).GenerateMesh(maxh=0.3, quad_dominated=True))
    mesh.Curve(2)
    return mesh


def _curved_half_disk_mesh():
    # A half disk on the axis inside a larger half disk: curved arcs end on the
    # axis (s-matched curved axis edges) and the inner arc is an interface.
    from netgen.occ import Glue, MoveTo, OCCGeometry, WorkPlane, X
    half = MoveTo(0.0, -3.0).Rectangle(3.0, 6.0).Face()
    inner = WorkPlane().Circle(1.0).Face() * half
    outer = (WorkPlane().Circle(2.0).Face() * half) - inner
    inner.faces.name = "c"
    outer.faces.name = "a"
    mesh = ng.Mesh(OCCGeometry(Glue([inner, outer]), dim=2).GenerateMesh(maxh=0.45))
    mesh.Curve(2)
    return mesh


def _rect_quad_mesh(r0):
    return structured_rect_vol_mesh(r0, r0 + 1.0, 0.0, 1.0, quads=True, nx=2, ny=2,
                                    stem="axifem_vdof_rect")


MATRIX_CASES = {
    # name: (mesh factory, order, curvedquad)
    "p1_off_axis": (_off_axis_trig_mesh, 1, False),
    "p2_off_axis": (_off_axis_trig_mesh, 2, False),
    "p1_chord_fold_triangle": (lambda: _mesh([(1, 0), (2, 1), (3, 2.1)], [(0, 1, 2)]), 1, False),
    "p2_chord_fold_triangle": (lambda: _mesh([(1, 0), (2, 1), (3, 2.1)], [(0, 1, 2)]), 2, False),
    "q1_off_axis": (lambda: _rect_quad_mesh(0.5), 1, False),
    "q2aa_off_axis": (lambda: _rect_quad_mesh(0.5), 2, False),
    "q2_general_skew": (lambda: _skew_quad_mesh(0.5), 2, True),
    "q1_axis": (lambda: _rect_quad_mesh(0.0), 1, False),
    "q2aa_axis": (lambda: _rect_quad_mesh(0.0), 2, False),
    "q2_general_axis": (lambda: _rect_quad_mesh(0.0), 2, True),
    # oblique edges ending on the axis (target map with a z shift)
    "p1_oblique_axis_triangle": (lambda: _mesh([(0, 0), (1, 0), (1, 1)], [(0, 1, 2)]), 1, False),
    "p2_oblique_axis_triangle": (lambda: _mesh([(0, 0), (1, 0), (1, 1)], [(0, 1, 2)]), 2, False),
    "p1_axis_triangles": (lambda: _boundary_mesh_triangles(), 1, False),
    "p2_axis_triangles": (lambda: _boundary_mesh_triangles(), 2, False),
    "q2_general_skew_axis": (lambda: _skew_quad_mesh(0.0), 2, True),
    # curved (Curve(2)) edges off the axis (blended target s)
    "p1_curved": (lambda: _curved_two_domain_mesh(), 1, False),
    "p2_curved": (lambda: _curved_two_domain_mesh(), 2, False),
    "q2_curved": (_curved_quad_sector_mesh, 2, True),
    "p1_curved_axis": (_curved_half_disk_mesh, 1, False),
    "p2_curved_axis": (_curved_half_disk_mesh, 2, False),
}


def _dense(bf):
    return np.array(bf.mat.ToDense())


@pytest.mark.parametrize("name", list(MATRIX_CASES))
def test_native_matrices_equal_symbolic_forms_of_the_evaluated_basis(name):
    factory, order, curvedquad = MATRIX_CASES[name]
    mesh = factory()
    fes = H1Henrotte(mesh, order=order, curvedquad=curvedquad)
    u, v = fes.TnT()
    x = ng.x
    # curved elements make the integrands rational: bonus 12 leaves ~7e-7
    bonus = ng.dx(bonus_intorder=20)
    k_sym = ng.BilinearForm(fes, check_unused=False)
    k_sym += (2 * math.pi / MU) * (ng.grad(u)[1] * ng.grad(v)[1]
                                   + (ng.grad(u)[0] + u / x) * (ng.grad(v)[0] + v / x)) * x * bonus
    m_sym = ng.BilinearForm(fes, check_unused=False)
    m_sym += 2 * math.pi * SIGMA * u * v * x * bonus
    k_nat = ng.BilinearForm(fes, check_unused=False)
    k_nat += AxiHenrotteStiffnessBFI(ng.CF(MU))
    m_nat = ng.BilinearForm(fes, check_unused=False)
    m_nat += AxiHenrotteSigmaMassBFI(ng.CF(SIGMA))
    for bf in (k_sym, m_sym, k_nat, m_nat):
        bf.Assemble()
    for label, nat, sym in (("stiffness", k_nat, k_sym), ("sigma mass", m_nat, m_sym)):
        a, b = _dense(nat), _dense(sym)
        rel = np.linalg.norm(a - b) / np.linalg.norm(b)
        assert rel <= MATRIX_RTOL, (
            f"{name} {label}: native vs symbolic relative difference {rel:.3e}, "
            f"expected <= {MATRIX_RTOL:.0e} (norms {np.linalg.norm(a):.6e} vs {np.linalg.norm(b):.6e})")


def _boundary_mesh_triangles():
    # Axis (r = 0), a horizontal edge ending on the axis, an off-axis oblique
    # edge, an oblique edge ending on the axis and a vertical edge.
    geo = SplineGeometry()
    pts = [geo.AppendPoint(*p) for p in [(0.0, 0.0), (1.0, 0.0), (1.4, 0.8), (1.4, 1.2), (0.0, 1.6)]]
    for i in range(5):
        geo.Append(["line", pts[i], pts[(i + 1) % 5]], bc="b%d" % i)
    return ng.Mesh(geo.GenerateMesh(maxh=0.5))


BOUNDARY_CASES = [
    ("triangles-p1", _boundary_mesh_triangles, 1, False),
    ("triangles-p2", _boundary_mesh_triangles, 2, False),
    ("axis-quads-q1", lambda: _rect_quad_mesh(0.0), 1, False),
    ("axis-quads-q2aa", lambda: _rect_quad_mesh(0.0), 2, False),
    ("skew-quads-q2", lambda: _skew_quad_mesh(0.0), 2, True),
    ("curved-p1", _curved_two_domain_mesh, 1, False),
    ("curved-p2", _curved_two_domain_mesh, 2, False),
    ("curved-q2", _curved_quad_sector_mesh, 2, True),
    ("curved-axis-p1", _curved_half_disk_mesh, 1, False),
    ("curved-axis-p2", _curved_half_disk_mesh, 2, False),
]


@pytest.mark.parametrize("name,factory,order,curvedquad", BOUNDARY_CASES,
                         ids=[c[0] for c in BOUNDARY_CASES])
def test_boundary_trace_equals_volume_value(name, factory, order, curvedquad):
    mesh = factory()
    fes = H1Henrotte(mesh, order=order, curvedquad=curvedquad)
    gf = ng.GridFunction(fes)
    vol_of = {}
    for el in mesh.Elements(ng.VOL):
        for vv in el.vertices:
            vol_of.setdefault(vv.nr, []).append(el)
    ts = np.linspace(0.04, 0.96, 13)
    probes = []
    for bel in mesh.Elements(ng.BND):
        a, b = [vv.nr for vv in bel.vertices]
        vels = [el for el in vol_of[a] if b in [vv.nr for vv in el.vertices]]
        assert len(vels) in (1, 2), f"segment ({a},{b}) has {len(vels)} volume neighbours"
        tb = mesh.GetTrafo(ng.ElementId(ng.BND, bel.nr))
        for vel in vels:   # both sides of an interface segment
            verts = [vv.nr for vv in vel.vertices]
            ra = np.array(REF_VERTICES[vel.type][verts.index(a)])
            rb = np.array(REF_VERTICES[vel.type][verts.index(b)])
            tv = mesh.GetTrafo(ng.ElementId(ng.VOL, vel.nr))
            for t in ts:
                # NGSolve's segment reference: t = 1 at vertex 0 (a), t = 0 at vertex 1 (b).
                pb = tb(float(t))
                pv = tv(*(t * ra + (1 - t) * rb))
                assert np.allclose(pb.point[:2], pv.point[:2], atol=1e-13), (pb.point, pv.point)
                # keep the transformations alive: a mapped point refers to its trafo
                probes.append((pb, pv, a, b, float(t), tb, tv))
    worst = (0.0, None)
    for dof in range(fes.ndof):
        gf.vec[:] = 0.0
        gf.vec[dof] = 1.0
        for pb, pv, a, b, t, _tb, _tv in probes:
            d = abs(gf(pb) - gf(pv))
            if d > worst[0]:
                worst = (d, (dof, a, b, t))
    assert worst[0] < TRACE_LIMIT, (
        f"{name}: boundary trace differs from the volume value by {worst[0]:.3e} "
        f"(DOF, segment vertices, t = {worst[1]}), expected < {TRACE_LIMIT:.0e}")


AXIS_CASES = [
    ("p1", _boundary_mesh_triangles, 1, False),
    ("p2", _boundary_mesh_triangles, 2, False),
    ("q1", lambda: _rect_quad_mesh(0.0), 1, False),
    ("q2aa", lambda: _rect_quad_mesh(0.0), 2, False),
    ("q2_general", lambda: _skew_quad_mesh(0.0), 2, True),
    ("p2_curved_axis", _curved_half_disk_mesh, 2, False),
]


@pytest.mark.parametrize("name,factory,order,curvedquad", AXIS_CASES,
                         ids=[c[0] for c in AXIS_CASES])
def test_axis_limit_of_the_basis(name, factory, order, curvedquad):
    mesh = factory()
    fes = H1Henrotte(mesh, order=order, curvedquad=curvedquad)
    gf = ng.GridFunction(fes)
    rng = np.random.default_rng(7)
    gf.vec.FV().NumPy()[:] = rng.standard_normal(fes.ndof)
    gr = ng.grad(gf)
    extent = max(max(abs(c) for c in mesh[v].point) for v in mesh.vertices)
    for z in (0.23, 0.61):
        a0 = gf(mesh(0.0, z))
        g0 = np.array(gr(mesh(0.0, z)))
        assert np.all(np.isfinite([a0, *g0])), (name, z, a0, g0)
        assert a0 == 0.0 and g0[1] == 0.0, f"{name}: A_phi = {a0}, dA/dz = {g0[1]} on the axis"
        scale = max(1.0, abs(g0[0]))
        for eps in (1e-3, 1e-6, 1e-9):
            mp = mesh(eps, z)
            a = gf(mp)
            g = np.array(gr(mp))
            r = ng.x(mp)   # radius of the point NGSolve located
            # A = (dA/dr)(0) r + O(r^2); the gradient approaches its axis value
            # linearly.  Floors: round-off of the nodal sums (~1e-14 absolute)
            # and, on curved maps, of the (s, z) map inversion, whose location
            # round-off ~1e-16 extent enters the gradient divided by r.
            inv = 1e-16 * extent / r
            assert abs(a - g0[0] * r) <= 1e-2 * scale * r * (r / 1e-3) + 1e-14 * scale, (
                f"{name}: A/r = {a / r:.9e} at r = {r:g}, axis dA/dr = {g0[0]:.9e}")
            assert abs(g[0] - g0[0]) <= 1e-2 * scale * (r / 1e-3) + (1e-12 + inv) * scale, (
                f"{name}: dA/dr = {g[0]:.9e} at r = {r:g}, axis value {g0[0]:.9e}")
            assert abs(g[1]) <= 1e-1 * scale * (r / 1e-3) + (1e-12 + inv) * scale, (
                f"{name}: dA/dz = {g[1]:.3e} at r = {r:g}, expected -> 0")


@pytest.mark.parametrize("factory,order,curvedquad", [
    (_curved_two_domain_mesh, 1, False),
    (_curved_two_domain_mesh, 2, False),
    (_curved_quad_sector_mesh, 2, True),
], ids=["p1", "p2", "q2"])
def test_curve3_mesh_is_refused(factory, order, curvedquad):
    """A Curve(3) map would be replaced by its quadratic interpolant: refused."""
    mesh = factory()
    mesh.Curve(3)
    gf = ng.GridFunction(H1Henrotte(mesh, order=order, curvedquad=curvedquad))
    with pytest.raises(Exception, match="degree > 2"):
        gf.Set(ng.x)


FD_CASES = [
    ("p1_curved", 1, False),
    ("p2_curved", 2, False),
    ("q2_curved", 2, True),
    ("p1_curved_axis", 1, False),
    ("p2_curved_axis", 2, False),
    ("p1_oblique_axis_triangle", 1, False),
    ("p2_oblique_axis_triangle", 2, False),
    ("p1_axis_triangles", 1, False),
    ("p2_axis_triangles", 2, False),
    ("q2_general_skew_axis", 2, True),
    ("p2_off_axis", 2, False),
]
FD_RTOL = 1.0e-6


@pytest.mark.parametrize("name,order,curvedquad", FD_CASES, ids=[c[0] for c in FD_CASES])
def test_gradient_is_the_derivative_of_the_value(name, order, curvedquad):
    """ng.grad equals central differences of the value, every DOF, interior points.

    Differences are taken along the reference directions through the element
    transformation (NGSolve's point location on curved elements is
    approximate): grad . dx/dxi = d value / d xi."""
    mesh = MATRIX_CASES[name][0]()
    fes = H1Henrotte(mesh, order=order, curvedquad=curvedquad)
    gf = ng.GridFunction(fes)
    grad = ng.grad(gf)
    worst = 0.0
    h = 1e-6
    for el in mesh.Elements(ng.VOL):
        ei = ng.ElementId(ng.VOL, el.nr)
        tr = mesh.GetTrafo(ei)
        refs = [(0.5, 0.5), (0.3, 0.2), (0.8, 0.6)] if el.type == ng.ET.QUAD else [(1 / 3, 1 / 3), (0.2, 0.55)]

        def pos(a, b):
            return np.array(tr(a, b).point[:2])

        for xi, eta in refs:
            jx = (pos(xi + h, eta) - pos(xi - h, eta)) / (2 * h)
            jy = (pos(xi, eta + h) - pos(xi, eta - h)) / (2 * h)
            for dof in fes.GetDofNrs(ei):
                gf.vec[:] = 0.0
                gf.vec[dof] = 1.0
                g = np.array(grad(tr(xi, eta)))
                fd = np.array([(gf(tr(xi + h, eta)) - gf(tr(xi - h, eta))) / (2 * h),
                               (gf(tr(xi, eta + h)) - gf(tr(xi, eta - h))) / (2 * h)])
                an = np.array([g @ jx, g @ jy])
                err = float(np.linalg.norm(an - fd)) / max(1.0, float(np.linalg.norm(fd)))
                worst = max(worst, err)
                assert err <= FD_RTOL, (
                    f"{name} element {el.nr} dof {dof} at ref ({xi:.3f}, {eta:.3f}): "
                    f"grad along the reference directions {an.tolist()} vs finite difference "
                    f"{fd.tolist()} (relative {err:.2e}, expected <= {FD_RTOL:.0e})")
    assert worst <= FD_RTOL
