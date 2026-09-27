"""Workpieces with holes: bores and cross holes through the IH thermal chain.

A bore turns the workpiece into a tube: the eddy-current field reaches the
bore through the wall, the meridian section of the surface is a closed loop
that never touches the axis, the bore wall is a heated or cooled surface of
nonzero revolved area, and depth rays can cross the hole.  A cross hole
breaks the rotational symmetry and must be refused wherever that symmetry is
assumed.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))

from radia import ih_axisym_coupled as C  # noqa: E402
from radia import ih_thermal, ih_thermal_post  # noqa: E402

MU0 = C.MU0
A_IN, A_OUT, L = 0.012, 0.02, 0.01          # tube radii, slice height
RC, TC, ROUT = 0.024, 5e-4, 0.03


def tube_axial_eddy_loss(H0, a, b, omega, sigma, mu_r):
    """Joule loss per unit length of an infinite tube a < r < b in an axial
    field of amplitude H0 outside it; the bore carries a uniform field.

    H_z = C1 I0(kr) + C2 K0(kr) in the wall (k^2 = j omega mu sigma), H_z = Hi
    in the bore with E_phi = -j omega mu0 Hi r / 2; H_z and E_phi are
    continuous at r = a and H_z(b) = H0.  Returns (P_joule, P_poynting).
    """
    from scipy import integrate, special

    k = np.sqrt(1j * omega * MU0 * mu_r * sigma)
    I0, I1 = special.iv(0, k * a), special.iv(1, k * a)
    K0, K1 = special.kv(0, k * a), special.kv(1, k * a)
    # unknowns C1, C2, Hi
    M = np.array([
        [I0, K0, -1.0],
        # E_phi(a) = -(1/sigma) dH/dr = -(k/sigma)(C1 I1 - C2 K1)
        [-(k / sigma) * I1, (k / sigma) * K1, 1j * omega * MU0 * a / 2],
        [special.iv(0, k * b), special.kv(0, k * b), 0.0]], complex)
    c1, c2, _hi = np.linalg.solve(M, np.array([0.0, 0.0, H0], complex))

    def J(r):                          # J_phi = -dH_z/dr
        return -k * (c1 * special.iv(1, k * r) - c2 * special.kv(1, k * r))

    P = integrate.quad(lambda r: abs(J(r)) ** 2 / (2 * sigma) * 2 * math.pi
                       * r, a, b, limit=400, epsabs=0, epsrel=1e-12)[0]
    E_b = J(b) / sigma
    P_s = -0.5 * (E_b * np.conj(H0)).real * 2 * math.pi * b
    return P, P_s


def _tube_slice(delta):
    """Axial slice of an infinite tube inside a current sheet; the bore air
    touches the axis, the outer air does not."""
    from netgen.occ import Glue, MoveTo, OCCGeometry, X
    from ngsolve import Mesh

    wall = A_OUT - A_IN
    sk = min(4 * delta, wall / 2)
    parts = []
    for r0, w, h in ((A_IN, sk, min(delta / 3, 1e-3)),
                     (A_IN + sk, wall - 2 * sk, 6e-4),
                     (A_OUT - sk, sk, min(delta / 3, 1e-3))):
        if w <= 0:
            continue
        f = MoveTo(r0, 0).Rectangle(w, L).Face()
        f.faces.name = "wp"
        f.maxh = h
        parts.append(f)
    coil = MoveTo(RC, 0).Rectangle(TC, L).Face()
    coil.faces.name = "coil"
    coil.maxh = TC / 2
    air = MoveTo(0, 0).Rectangle(ROUT, L).Face()
    air.faces.name = "air"
    air.edges.Min(X).name = "axis"
    for f in parts + [coil]:
        air = air - f
    return Mesh(OCCGeometry(Glue([air, *parts, coil]), dim=2)
                .GenerateMesh(maxh=2e-3))


@pytest.mark.parametrize("sigma, mu_r, rtol", [(1.0e6, 1.0, 1e-6),
                                               (5.8e6, 100.0, 1e-4)])
def test_tube_eddy_loss_matches_the_bessel_solution(sigma, mu_r, rtol):
    """The bore field (thick wall: 8 mm against a 6 mm skin depth) and the
    shielded case (0.26 mm skin depth, both wall faces resolved)."""
    f = 7000.0
    omega = 2 * math.pi * f
    delta = math.sqrt(2 / (omega * MU0 * mu_r * sigma))
    mesh = _tube_slice(delta)
    em = C.AxisymEddyCurrent(mesh, frequency=f, workpiece="wp",
                             coils={"coil": 1000.0}, dirichlet="axis")
    n = len(em.workpiece_elements)
    ratio, _ = em.skin_resolution(np.full(n, delta))
    assert ratio <= 1.0
    rec = em.solve(np.full(n, sigma), np.full(n, mu_r))
    P, P_s = tube_axial_eddy_loss(1000.0 / L, A_IN, A_OUT, omega, sigma,
                                  mu_r)
    assert P_s == pytest.approx(P, rel=1e-8)       # the reference itself
    # which tends to the solid cylinder as the bore closes
    from radia.analytical_formulas.induction_heating import (
        cylinder_axial_eddy_loss)
    P_solid, _ = tube_axial_eddy_loss(1.0, 1e-7, A_OUT, omega, sigma, mu_r)
    assert P_solid == pytest.approx(cylinder_axial_eddy_loss(
        1.0, A_OUT, omega, sigma, mu_r), rel=1e-10)
    assert rec["P_joule_W"] == pytest.approx(P * L, rel=rtol)
    assert abs(rec["power_balance_relative_error"]) < 1e-7


# --- thermal side ----------------------------------------------------------

H3 = 0.03


def _tube3d(maxh, bore_name="bore"):
    """3D tube A_IN < r < A_OUT, 0 < z < H3."""
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import Mesh

    outer = Cylinder(Axes(Pnt(0, 0, 0), Z), r=A_OUT, h=H3)
    outer.faces.name = "outer"
    outer.faces.Max(Z).name = "top"
    outer.faces.Min(Z).name = "bottom"
    inner = Cylinder(Axes(Pnt(0, 0, -0.001), Z), r=A_IN, h=H3 + 0.002)
    inner.faces.name = bore_name
    tube = outer - inner
    tube.solids.name = "wp"
    return Mesh(OCCGeometry(tube).GenerateMesh(maxh=maxh))


def _p1(mesh, expr):
    from ngsolve import GridFunction, H1
    gf = GridFunction(H1(mesh, order=1))
    gf.Set(expr)
    return gf


def test_tube_meridian_section_is_a_closed_loop():
    """The section of a tube surface never reaches the axis; the average
    of an axisymmetric source is reproduced on the bore and the outer wall,
    and the azimuthal part of 2 + x is removed, with the power kept."""
    from ngsolve import x, z
    mesh = _tube3d(0.003)
    names = ["outer", "bore", "top", "bottom"]
    assert set(names) <= set(mesh.GetBoundaries())
    pts = ih_thermal.mesh_vertices(mesh)
    vn = ih_thermal.boundary_vertex_numbers(mesh, names)
    for expr, expect in ((1.0 + 1e3 * z * z, 1.0 + 1e3 * pts[vn, 2] ** 2),
                         (2.0 + 10.0 * x, np.full(len(vn), 2.0))):
        prof = ih_thermal.AxisymmetricSurfaceProfile.from_gridfunction(
            mesh, _p1(mesh, expr), names)
        assert len(prof.chains) == 1 and prof.chains[0].closed
        assert prof.power == pytest.approx(prof.source_power, rel=1e-12)
        vals, dist = prof.evaluate_xyz(pts[vn], what="tube vertices")
        # the rims (kink between wall and annulus) carry the largest error
        np.testing.assert_allclose(vals, expect, atol=0.02 * np.max(expect))


def _write_bore_source(tmp_path, q0):
    from ngsolve import BND, CF, Integrate
    mesh = _tube3d(0.003)
    vol = str(tmp_path / "em_tube.vol")
    mesh.ngmesh.Save(vol)
    from ngsolve import Mesh
    mesh = Mesh(vol)
    gf = _p1(mesh, CF(q0))
    sol = str(tmp_path / "q.sol")
    gf.Save(sol)
    P = float(Integrate(gf, mesh, BND, definedon=mesh.Boundaries("bore")))
    ih_thermal.write_field_sidecar(
        sol, mesh_path=vol, mesh=mesh, fes_order=1,
        quantity=ih_thermal.QSURF_QUANTITY, unit=ih_thermal.QSURF_UNIT,
        boundaries=["bore"], extra={"P_wp_W": P})
    return sol, vol, P


def test_axisymmetric_heat_through_the_bore_wall(tmp_path):
    """Heating only the bore wall (r > 0, nonzero revolved area) of an
    annular meridian: the transferred power matches the EM power and all of
    it enters through the bore."""
    import calc_heat_axisym
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh

    q0 = 1.0e5
    sol, vol, P = _write_bore_source(tmp_path, q0)
    geo = SplineGeometry()
    geo.AddRectangle((A_IN, 0.0), (A_OUT, H3),
                     bcs=("bottom", "outer", "top", "bore"))
    m2 = Mesh(geo.GenerateMesh(maxh=0.0015))
    res = calc_heat_axisym.solve_heat_axisym(
        "<annulus>", material="custom", rho=7800.0, cp=500.0, k=40.0,
        h_conv=0.0, heat_flux_boundaries="bore", qsurf_sol=sol, em_vol=vol,
        em_heat_boundaries="bore", n_phi_samples=32, dt=0.5, t_end=2.0,
        fes_order=2, _wp_mesh=m2, _write_solution=False)
    assert "error" not in res, res
    first = res["qsurf_projection"]["initial"]
    assert abs(first["power_balance"]["relative_error"]) < 5e-3
    # the faceted EM bore has slightly less area than the true one
    exact = q0 * 2 * math.pi * A_IN * H3
    assert P == pytest.approx(exact, rel=5e-3)
    assert res["q_surf_int_W"] == pytest.approx(P, rel=5e-3)
    assert res["Q_input_J"] == pytest.approx(res["q_surf_int_W"] * 2.0,
                                             rel=1e-10)


def test_3d_heat_on_an_independently_meshed_tube(tmp_path):
    """Direct transfer onto a thermal tube meshed independently of the EM
    tube (different facets on the curved bore)."""
    import calc_heat

    q0 = 1.0e5
    sol, vol, P = _write_bore_source(tmp_path, q0)
    thermal = _tube3d(0.0022)
    res = calc_heat.solve_heat(
        "<tube>", material="custom", rho=7800.0, cp=500.0, k=40.0,
        h_conv=0.0, heat_flux_boundaries="bore", qsurf_sol=sol, em_vol=vol,
        em_heat_boundaries="bore", dt=0.5, t_end=1.0, fes_order=1,
        _wp_mesh=thermal, _write_solution=False)
    assert "error" not in res, res
    first = res["qsurf_projection"]["initial"]
    assert abs(first["power_balance"]["relative_error"]) < 5e-3
    assert first["max_transfer_distance_m"] < 2e-4      # chord sagitta


def test_phi_average_refuses_a_cross_hole():
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, X, Z
    from ngsolve import Mesh, z

    body = Cylinder(Axes(Pnt(0, 0, 0), Z), r=A_OUT, h=H3)
    body.faces.name = "surf"
    hole = Cylinder(Pnt(-0.03, 0, H3 / 2), X, r=0.004, h=0.06)
    hole.faces.name = "surf"
    mesh = Mesh(OCCGeometry(body - hole).GenerateMesh(maxh=0.003))
    with pytest.raises(ValueError, match="not a body of revolution|branches"):
        ih_thermal.AxisymmetricSurfaceProfile.from_gridfunction(
            mesh, _p1(mesh, 1.0 + z), ["surf"])


def test_depth_ray_stops_at_the_hole_wall():
    """A hot ray that reaches a hole ends there ('through') at the hole
    wall, not at the far side of the part; a cooling ray gives its depth."""
    from netgen.occ import Axes, Box, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import CF, Mesh, x

    block = Box(Pnt(0, 0, 0), Pnt(0.02, 0.02, 0.01))
    block.faces.name = "surf"
    hole = Cylinder(Axes(Pnt(0.01, 0.01, -0.001), Z), r=0.003, h=0.012)
    hole.faces.name = "hole"
    mesh = Mesh(OCCGeometry(block - hole).GenerateMesh(maxh=0.0015))
    hot = ih_thermal_post.case_depth(
        mesh, _p1_order2(mesh, CF(1000.0)), 850.0, span=0.02,
        origins=[[0.0, 0.01, 0.005]], normals=[[1.0, 0.0, 0.0]])
    assert hot["status"] == ["through"]
    # the hole wall is at x = 0.007 (faceted: within the chord sagitta)
    assert hot["depth_m"][0] == pytest.approx(0.007, abs=1.5e-4)
    cool = ih_thermal_post.case_depth(
        mesh, _p1_order2(mesh, 1000.0 - 1.0e5 * x), 850.0, span=0.02,
        origins=[[0.0, 0.01, 0.005]], normals=[[1.0, 0.0, 0.0]])
    assert cool["status"] == ["ok"]
    assert cool["depth_m"][0] == pytest.approx(1.5e-3, abs=1e-6)


def _p1_order2(mesh, expr):
    from ngsolve import GridFunction, H1
    gf = GridFunction(H1(mesh, order=2))
    gf.Set(expr)
    return gf


def test_coupled_transient_on_a_tube():
    """The staggered EM-heat run on a tube: the skin gate covers the bore
    interface, the energy balance closes, and the Joule power heats the
    wall from both faces' fields."""
    from radia import ih_heat_transient as iht
    from radia import ih_thermal_material as itm

    sigma = 6.0e6
    delta = math.sqrt(2 / (2 * math.pi * 7000.0 * MU0 * sigma))
    mesh = _tube_slice(delta)
    table = C.EMMaterialTable(T=[0.0, 2000.0], sigma=[sigma, sigma],
                              mu_r=[1.0, 1.0], source="constant")
    result, gfT, em = C.run_coupled(
        mesh, frequency=7000.0, workpiece="wp", coils={"coil": 2000.0},
        dirichlet="axis", em_material=table,
        thermal_material=itm.ThermalMaterial.constant(7800.0, 600.0, 35.0),
        boundaries=iht.HeatBoundaryTerms(), dt=0.1, t_end=0.3)
    nl = result["nonlinear_transient"]
    assert abs(nl["energy_balance_relative_error"]) < 1e-6
    assert all(h["skin_resolution_h_over_delta"] <= 1.0
               for h in result["em_history"])
    P, _ = tube_axial_eddy_loss(2000.0 / L, A_IN, A_OUT,
                                2 * math.pi * 7000.0, sigma, 1.0)
    assert result["em_history"][0]["P_joule_W"] == pytest.approx(P * L,
                                                                 rel=1e-4)
