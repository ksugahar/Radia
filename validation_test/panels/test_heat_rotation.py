"""Heat source on a rotating workpiece, and the circumferential average.

A coil fixed in the world illuminates the part; a point ``p`` of the part
sits at ``R(theta) p`` when the part has turned by ``theta``.  With the world
field ``q = 2 + x`` the body-frame source is exactly
``q(p, theta) = 2 + (R(theta) p)_x``, which the tests use as the reference.

* A body of revolution needs one EM solution, turned with the part.
* A part that is not (a cube, a cross hole) needs one EM solution per rotor
  angle over its symmetry period (``--rotor-states``); a single turned
  solution is refused.
* Each time step applies the exact average over the angles it sweeps, so a
  step of one revolution gives the revolution average rather than one
  stroboscopic sample.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))


def _write_em_pair(mesh, gf, td, name="q", boundary="default"):
    """Save an EM field with the sidecar the transfer requires."""
    from ngsolve import BND, Integrate
    from radia import ih_thermal

    sol = os.path.join(td, f"{name}.sol")
    vol = os.path.join(td, f"{name}.vol")
    gf.Save(sol)
    mesh.ngmesh.Save(vol)
    power = float(Integrate(gf, mesh, BND,
                            definedon=mesh.Boundaries(boundary)))
    ih_thermal.write_field_sidecar(
        sol, mesh_path=vol, mesh=mesh, fes_order=1,
        quantity=ih_thermal.QSURF_QUANTITY, unit=ih_thermal.QSURF_UNIT,
        boundaries=[boundary], extra={"P_wp_W": power})
    return sol, vol


def _world_field(mesh, expr):
    from ngsolve import GridFunction, H1
    gf = GridFunction(H1(mesh, order=1))
    gf.Set(expr)
    return gf


def _cylinder(maxh=0.25):
    from ngsolve import Mesh
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    return Mesh(OCCGeometry(Cylinder(Axes(Pnt(0, 0, -1), Z), r=1.0, h=2.0))
                .GenerateMesh(maxh=maxh))


def _cube(angle=0.0, maxh=0.4):
    """Cube [-1, 1]^3 turned by ``angle`` about z."""
    from ngsolve import Mesh
    from netgen.occ import Axis, Box, OCCGeometry, Pnt, Z
    box = Box(Pnt(-1, -1, -1), Pnt(1, 1, 1))
    if angle:
        box = box.Rotate(Axis((0, 0, 0), Z), math.degrees(angle))
    return Mesh(OCCGeometry(box).GenerateMesh(maxh=maxh))


def _exact(pts, theta):
    """2 + (R(theta) p)_x for points p of the part."""
    return 2.0 + pts[:, 0] * math.cos(theta) - pts[:, 1] * math.sin(theta)


def _exact_average(pts, t0, t1):
    c = math.sin(t1) - math.sin(t0)
    s = math.cos(t0) - math.cos(t1)
    return 2.0 + (pts[:, 0] * c - pts[:, 1] * s) / (t1 - t0)


def _single(tmp_path, mesh):
    """RotatingSurfaceSource from one EM solution (q = 2 + x) on ``mesh``."""
    import calc_heat
    from ngsolve import x
    from radia import ih_thermal

    sol, vol = _write_em_pair(mesh, _world_field(mesh, 2.0 + x),
                              str(tmp_path))
    args = calc_heat.qsurf_args(qsurf_sol=sol, em_vol=vol,
                                heat_flux_boundary_names=["default"])
    q_cf, _, audit = calc_heat._build_qsurf_source(mesh, args)
    return q_cf, ih_thermal.RotatingSurfaceSource(
        mesh, ["default"], q_cf, states=[(0.0, audit["_source"])],
        power_tolerance=audit["power_tolerance"])


def _surface(mesh, q_cf, rot):
    from radia import ih_thermal
    return ih_thermal.mesh_vertices(mesh)[rot.vnrs], rot.vnrs


def test_a_turned_body_of_revolution_follows_the_fixed_coil(tmp_path):
    mesh = _cylinder()
    q_cf, rot = _single(tmp_path, mesh)
    assert rot.kind == "turned-single-state" and len(rot.angles) >= 8
    pts, vn = _surface(mesh, q_cf, rot)
    for theta in (0.0, 0.3, math.pi / 2, 2.5, math.pi, 5.9):
        rot.at(theta)
        np.testing.assert_allclose(q_cf.vec.FV().NumPy()[vn],
                                   _exact(pts, theta), atol=0.03)


def test_a_step_applies_the_average_over_the_swept_angles(tmp_path):
    """The stroboscopic trap: a step of exactly one revolution sampled at its
    end sees the part where it started; the average sees the revolution."""
    mesh = _cylinder()
    q_cf, rot = _single(tmp_path, mesh)
    pts, vn = _surface(mesh, q_cf, rot)
    rot.average(0.0, 2 * math.pi)
    np.testing.assert_allclose(q_cf.vec.FV().NumPy()[vn], 2.0, atol=1e-3)
    rot.average(4 * math.pi, 8 * math.pi)            # two revolutions
    np.testing.assert_allclose(q_cf.vec.FV().NumPy()[vn], 2.0, atol=1e-3)
    for t0, t1 in ((0.2, 1.1), (5.5, 7.3), (1.0, 1.0 + 1e-3)):
        rot.average(t0, t1)
        np.testing.assert_allclose(q_cf.vec.FV().NumPy()[vn],
                                   _exact_average(pts, t0, t1), atol=0.03)
    with pytest.raises(ValueError, match="precede"):
        rot.average(1.0, 0.5)


def test_a_part_that_is_not_of_revolution_needs_rotor_states(tmp_path):
    with pytest.raises(ValueError, match="--rotor-states"):
        _single(tmp_path, _cube())


def _cube_states(tmp_path, angles, period):
    """World-frame EM solutions q = 2 + x on the cube turned by each angle."""
    from ngsolve import x
    states = []
    for i, a in enumerate(angles):
        mesh = _cube(a)
        sol, vol = _write_em_pair(mesh, _world_field(mesh, 2.0 + x),
                                  str(tmp_path), name=f"state{i}")
        states.append({"angle_rad": a, "qsurf_sol": os.path.basename(sol),
                       "em_vol": os.path.basename(vol)})
    man = tmp_path / "rotor.json"
    man.write_text(json.dumps({
        "schema": "radia.ih-rotor-states/1", "axis": "z",
        "period_rad": period, "frame": "world", "states": states}),
        encoding="utf-8")
    return str(man)


def _rotor(tmp_path, angles, period=math.pi / 2, tolerance=None):
    from radia import ih_thermal
    from ngsolve import GridFunction, H1
    man = _cube_states(tmp_path, angles, period)
    states, meta = ih_thermal.load_rotor_states(man)
    thermal = _cube(maxh=0.3)                 # meshed independently
    gf = GridFunction(H1(thermal, order=1))
    rot = ih_thermal.RotatingSurfaceSource(
        thermal, ["default"], gf, states=states, power_tolerance=0.02,
        axis=meta["axis"], period=meta["period_rad"], frame=meta["frame"],
        angle_step_tolerance=tolerance)
    return thermal, gf, rot


def test_rotor_states_of_a_four_fold_part(tmp_path):
    """Four states over the quarter-turn period of a cube give a knot every
    pi/8 over the revolution; values and step averages follow the exact
    body-frame source."""
    angles = [k * math.pi / 8 for k in range(4)]
    thermal, gf, rot = _rotor(tmp_path, angles)
    assert rot.kind == "rotor-states" and rot.n_fold == 4
    assert len(rot.angles) == 16
    pts, vn = _surface(thermal, gf, rot)
    # at knots (the states, and a state turned by one period) the transfer
    # is exact up to facets; between knots the source is the linear
    # interpolation of a cosine sampled every pi/8
    for theta, tol in [(a, 0.02) for a in angles] + [
            (math.pi / 2 + math.pi / 8, 0.02), (3.7, 0.08)]:
        rot.at(theta)
        np.testing.assert_allclose(gf.vec.FV().NumPy()[vn],
                                   _exact(pts, theta), atol=tol)
    rot.average(0.0, 2 * math.pi)
    np.testing.assert_allclose(gf.vec.FV().NumPy()[vn], 2.0, atol=1e-2)


def test_coarse_rotor_states_are_refused(tmp_path):
    with pytest.raises(ValueError, match="too coarse"):
        _rotor(tmp_path, [0.0, math.pi / 4], tolerance=0.1)


def test_a_declared_symmetry_the_part_does_not_have_is_refused(tmp_path):
    """Declaring the cube 8-fold makes the knot at pi/4 read the state at 0
    on the cube turned by pi/4, which is not on the cube."""
    with pytest.raises(ValueError, match="farther than"):
        _rotor(tmp_path, [0.0, math.pi / 8], period=math.pi / 4)


def test_rotor_states_through_the_heat_solver(tmp_path):
    """End to end: a rotating cube heated by a fixed source for exactly one
    revolution receives the revolution-average energy."""
    import calc_heat

    angles = [k * math.pi / 8 for k in range(4)]
    man = _cube_states(tmp_path, angles, math.pi / 2)
    common = dict(material="custom", rho=1.0, cp=1.0, k=1.0, h_conv=0.0,
                  heat_flux_boundaries="default", dt=0.25, t_end=1.0,
                  fes_order=1, _wp_mesh=_cube(maxh=0.3),
                  _write_solution=False)
    res = calc_heat.solve_heat("<cube>", rotor_states=man,
                               rotation_rpm=60.0, **common)
    assert "error" not in res, res
    rot = res["qsurf_projection"]["rotation"]
    assert rot["kind"] == "rotor-states" and rot["knots"] == 16
    assert rot["swept_per_step_rad"] == pytest.approx(math.pi / 2)
    # each quarter-turn step of 2 + x: the part's x-extent averages to a
    # field with the same power as 2 (the x term integrates to zero)
    area = 24.0
    assert res["Q_input_J"] == pytest.approx(2.0 * area * 1.0, rel=2e-2)
    with_rpm0 = calc_heat.solve_heat("<cube>", rotor_states=man,
                                     rotation_rpm=0.0, **common)
    assert "--rotation-rpm" in with_rpm0["error"]
    both = calc_heat.solve_heat("<cube>", rotor_states=man, qsurf_sol="x",
                                em_vol="y", rotation_rpm=60.0, **common)
    assert "does not apply" in both["error"]


# --- circumferential average ------------------------------------------------

@pytest.fixture(scope="module")
def cylinder_mesh():
    """A body of revolution about z: radius 1, height 2, centred on z=0."""
    return _cylinder()


def _phi_average_on(mesh, expression, td):
    import calc_heat

    sol, vol = _write_em_pair(mesh, _world_field(mesh, expression), td)
    args = calc_heat.qsurf_args(
        qsurf_sol=sol, em_vol=vol, heat_flux_boundary_names=["default"],
        q_phi_average=True)
    return calc_heat._build_qsurf_source(mesh, args)


def test_phi_average_removes_the_azimuthal_part(cylinder_mesh, tmp_path):
    """The circumferential average of q = 2 + x (x = r cos phi) is 2, and it
    is a static source."""
    from ngsolve import x
    from radia import ih_thermal
    q_cf, resample, audit = _phi_average_on(cylinder_mesh, 2.0 + x,
                                            str(tmp_path))
    assert resample is None and "_source" not in audit
    assert audit["mode"] == "phi-average"
    side = ih_thermal.boundary_vertex_numbers(cylinder_mesh, ["default"])
    vals = q_cf.vec.FV().NumPy()[side] - 2.0
    r = np.hypot(*ih_thermal.mesh_vertices(cylinder_mesh)[side, :2].T)
    # |x| reaches 1 on this surface.  Away from the axis the residual is
    # the facet asymmetry of each ring; within one element of the axis the
    # triangles straddle it and the meridian coordinate is not linear on
    # them, so the bound there is looser.
    assert np.max(np.abs(vals[r > 0.25])) < 0.02
    assert np.max(np.abs(vals)) < 0.07


def test_phi_average_keeps_an_axisymmetric_source(cylinder_mesh, tmp_path):
    """An already axisymmetric source is reproduced and its power kept."""
    from ngsolve import z
    from radia import ih_thermal
    q_cf, _, audit = _phi_average_on(cylinder_mesh, 1.0 + z * z,
                                     str(tmp_path))
    pts = ih_thermal.mesh_vertices(cylinder_mesh)
    side = ih_thermal.boundary_vertex_numbers(cylinder_mesh, ["default"])
    vals = q_cf.vec.FV().NumPy()[side]
    # The largest deviation (about 1.5 % of the value) is at the rims,
    # where one bin of width h/4 averages across the kink between the side
    # wall and the end cap.
    np.testing.assert_allclose(vals, 1.0 + pts[side, 2] ** 2, atol=0.04)
    balance = audit["initial"]["power_balance"]
    assert abs(balance["relative_error"]) < 5e-3
    assert audit["phi_average"]["profile_power_W"] == pytest.approx(
        audit["source_power_W"], rel=1e-12)


def test_phi_average_requires_spatial_not_uniform(cylinder_mesh):
    """--q-phi-average + --q-uniform is contradictory (a constant is
    already azimuthally uniform) and must fail loud."""
    import calc_heat

    args = calc_heat.qsurf_args(
        q_uniform=1.0e6, heat_flux_boundary_names=["default"],
        q_phi_average=True)
    with pytest.raises(ValueError, match="q-phi-average"):
        calc_heat._build_qsurf_source(cylinder_mesh, args)


def test_body_frame_states_match_the_world_frame(tmp_path):
    """Body frame: the part stays put and the coil turns by -a, so state k
    holds g(y) = f(R(a) y) on the unturned part.  It must give the same
    body-frame source as the world frame."""
    from ngsolve import GridFunction, H1, x, y
    from radia import ih_thermal

    angles = [k * math.pi / 8 for k in range(4)]
    states = []
    for i, a in enumerate(angles):
        mesh = _cube(0.0)
        field = 2.0 + x * math.cos(a) - y * math.sin(a)      # f(R(a) y)
        sol, vol = _write_em_pair(mesh, _world_field(mesh, field),
                                  str(tmp_path), name=f"body{i}")
        states.append({"angle_rad": a, "qsurf_sol": os.path.basename(sol),
                       "em_vol": os.path.basename(vol)})
    man = tmp_path / "rotor_body.json"
    man.write_text(json.dumps({
        "schema": "radia.ih-rotor-states/1", "axis": "z",
        "period_rad": math.pi / 2, "frame": "body", "states": states}),
        encoding="utf-8")
    loaded, meta = ih_thermal.load_rotor_states(str(man))
    thermal = _cube(maxh=0.3)
    gf = GridFunction(H1(thermal, order=1))
    rot = ih_thermal.RotatingSurfaceSource(
        thermal, ["default"], gf, states=loaded, power_tolerance=0.02,
        axis="z", period=meta["period_rad"], frame="body")
    pts, vn = _surface(thermal, gf, rot)
    for theta in (angles[2], math.pi + angles[1], 3 * math.pi / 2):
        rot.at(theta)
        np.testing.assert_allclose(gf.vec.FV().NumPy()[vn],
                                   _exact(pts, theta), atol=0.02)
