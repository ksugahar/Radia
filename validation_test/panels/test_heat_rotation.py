"""Unit tests for the 3D thermal solver's workpiece-rotation feature.

The feature (2026-05-20, v4.58.0) re-projects qsurf on the workpiece
body each timestep when ``--rotation-rpm > 0`` is set on calc_heat.py.
We verify the projection math (forward rotation around z) here without
running the full transient solve.  Spawning calc_heat end-to-end is
covered by validation_test/panels/test_heat_chain_golden.py (which runs slow).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

# Add src/radia/panels to import path for direct access to calc_heat.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PANEL_DIR = os.path.join(ROOT, "src", "radia", "panels")
sys.path.insert(0, PANEL_DIR)


def _surf_vnr_nearest(wp_mesh, target):
    """Find the surface (BND) vertex nr closest to ``target=(x,y,z)``."""
    tx, ty, tz = target
    surf_vnrs = set()
    for el in wp_mesh.Elements(__import__("ngsolve").BND):
        for v_nodeid in el.vertices:
            surf_vnrs.add(v_nodeid.nr)
    best_vnr, best_d = None, 1e9
    for vnr in surf_vnrs:
        p = wp_mesh.vertices[vnr].point
        d = (p[0] - tx) ** 2 + (p[1] - ty) ** 2 + (p[2] - tz) ** 2
        if d < best_d:
            best_d, best_vnr = d, vnr
    return best_vnr


@pytest.fixture(scope="module")
def synthetic_setup(tmp_path_factory):
    """Build a synthetic (em_mesh, gf_q_em, wp_mesh) trio:

    * EM mesh = unit cube
    * gf_q_em(x,y,z) = x      (linear in x so we can predict rotated values)
    * wp_mesh = same unit cube (so the surf vertices land inside em_mesh)

    A wp body point at (1,0,0) sees q_em(1,0,0)=1 at theta=0.
    At theta=pi/2 the same body point sits at world (0,1,0), and
    q_em(0,1,0)=0.  At theta=pi it sits at (-1,0,0), q_em=-1.
    """
    from ngsolve import (Mesh, H1, GridFunction, x as ng_x, TaskManager)
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import meshsize

    geo = OCCGeometry(Box(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
    mesh_ng = geo.GenerateMesh(maxh=0.4)
    em_mesh = Mesh(mesh_ng)
    wp_mesh = em_mesh  # share the mesh for simplicity

    fes_q = H1(em_mesh, order=1)
    gf_q_em = GridFunction(fes_q)
    with TaskManager():
        gf_q_em.Set(ng_x)  # q_em = x

        return em_mesh, gf_q_em, wp_mesh


def test_rotation_zero_matches_identity(synthetic_setup):
    """At theta=0 the rotated projection equals the original projection."""
    import calc_heat
    from types import SimpleNamespace

    em_mesh, gf_q_em, wp_mesh = synthetic_setup

    # Save gf_q_em to a temp .sol so _build_qsurf_cf can load it.
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        sol = os.path.join(td, "q.sol")
        vol = os.path.join(td, "wp.vol")
        gf_q_em.Save(sol)
        em_mesh.ngmesh.Save(vol)

        args = SimpleNamespace(
            q_uniform=None, qsurf_sol=sol, em_vol=vol, qsurf_order=1,
            heat_flux_boundary_names=["default"])
        # The helper walks BND vertices filtered by the resolved heat-flux role
        # to enumerate surf_vertex_nrs.
        q_cf, resample = calc_heat._build_qsurf_cf(wp_mesh, args)
        assert resample is not None

        # After initial projection at theta=0, the surface values
        # should approximate gf_q_em (which equals x) at the wp
        # surface.  Sample at a known point.
        from ngsolve import H1, GridFunction
        from ngsolve import TaskManager
        # gf_wp_q is q_cf itself.  Pick the surface vertex closest to (1,0,0).
        best_vnr = _surf_vnr_nearest(wp_mesh, (1.0, 0.0, 0.0))
        assert best_vnr is not None
        val_theta0 = q_cf.vec.FV()[best_vnr]
        # At theta=0, body point (~1,0,0) ↦ world (~1,0,0), q=x≈1.
        assert val_theta0 == pytest.approx(1.0, abs=0.05)


def test_rotation_pi_inverts_x_signed_value(synthetic_setup):
    """At theta=pi the q value at body (1,0,0) becomes q_em(-1,0,0)≈-1."""
    import calc_heat
    from types import SimpleNamespace
    import tempfile

    em_mesh, gf_q_em, wp_mesh = synthetic_setup

    with tempfile.TemporaryDirectory() as td:
        sol = os.path.join(td, "q.sol")
        vol = os.path.join(td, "wp.vol")
        gf_q_em.Save(sol)
        em_mesh.ngmesh.Save(vol)

        args = SimpleNamespace(
            q_uniform=None, qsurf_sol=sol, em_vol=vol, qsurf_order=1,
            heat_flux_boundary_names=["default"])
        q_cf, resample = calc_heat._build_qsurf_cf(wp_mesh, args)

        # Locate the surface vertex nearest to body (+1, 0, 0).
        best_vnr = _surf_vnr_nearest(wp_mesh, (1.0, 0.0, 0.0))

        # Resample with the body rotated by pi -- the body point
        # near (+1, 0, 0) now sits at world (-1, 0, 0) where q_em = -1.
        resample(math.pi)
        val_pi = q_cf.vec.FV()[best_vnr]
        assert val_pi == pytest.approx(-1.0, abs=0.05)


def test_rotation_pi_over_2_swaps_x_to_y(synthetic_setup):
    """At theta=pi/2 the q value at body (1,0,0) becomes q_em(0,1,0)≈0
    (since q_em=x and the world point at (0,1,0) has x=0)."""
    import calc_heat
    from types import SimpleNamespace
    import tempfile

    em_mesh, gf_q_em, wp_mesh = synthetic_setup

    with tempfile.TemporaryDirectory() as td:
        sol = os.path.join(td, "q.sol")
        vol = os.path.join(td, "wp.vol")
        gf_q_em.Save(sol)
        em_mesh.ngmesh.Save(vol)

        args = SimpleNamespace(
            q_uniform=None, qsurf_sol=sol, em_vol=vol, qsurf_order=1,
            heat_flux_boundary_names=["default"])
        q_cf, resample = calc_heat._build_qsurf_cf(wp_mesh, args)

        # Locate the surface vertex nearest to body (+1, 0, 0).
        best_vnr = _surf_vnr_nearest(wp_mesh, (1.0, 0.0, 0.0))

        resample(math.pi / 2.0)
        val_half = q_cf.vec.FV()[best_vnr]
        # At theta=pi/2, body (1,0,0) ↦ world (0,1,0), q_em = x = 0.
        assert val_half == pytest.approx(0.0, abs=0.05)


@pytest.fixture(scope="module")
def cylinder_mesh():
    """A body of revolution about z: radius 1, height 2, centred on z=0."""
    from ngsolve import Mesh
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z

    solid = Cylinder(Axes(Pnt(0, 0, -1), Z), r=1.0, h=2.0)
    return Mesh(OCCGeometry(solid).GenerateMesh(maxh=0.25))


def _phi_average_on(mesh, expression, td):
    import calc_heat
    from types import SimpleNamespace
    from ngsolve import H1, GridFunction

    gf = GridFunction(H1(mesh, order=1))
    gf.Set(expression)
    sol = os.path.join(td, "q.sol")
    vol = os.path.join(td, "em.vol")
    gf.Save(sol)
    mesh.ngmesh.Save(vol)
    args = SimpleNamespace(
        q_uniform=None, qsurf_sol=sol, em_vol=vol, qsurf_order=1,
        heat_flux_boundary_names=["default"], rotation_axis="z",
        q_phi_average=True)
    return calc_heat._build_qsurf_source(mesh, args)


def test_phi_average_removes_the_azimuthal_part(cylinder_mesh, tmp_path):
    """The circumferential average of q = x (= r cos phi) is zero, and it is
    a static source (no resampler)."""
    from ngsolve import x
    import ih_thermal

    q_cf, resample, audit = _phi_average_on(cylinder_mesh, x, str(tmp_path))
    assert resample is None
    assert audit["mode"] == "phi-average"
    side = ih_thermal.boundary_vertex_numbers(cylinder_mesh, ["default"])
    vals = q_cf.vec.FV().NumPy()[side]
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
    import ih_thermal

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


def test_phi_average_requires_spatial_not_uniform(synthetic_setup):
    """--q-phi-average + --q-uniform is contradictory (a constant is
    already azimuthally uniform) and must fail loud, per No-Fallback."""
    import calc_heat
    from types import SimpleNamespace

    _em, _gf, wp_mesh = synthetic_setup
    args = SimpleNamespace(
        q_uniform=1.0e6, qsurf_sol="", em_vol="", qsurf_order=1,
        heat_flux_boundary_names=["default"], rotation_axis="z",
        q_phi_average=True)
    with pytest.raises(ValueError, match="q-phi-average"):
        calc_heat._build_qsurf_cf(wp_mesh, args)
