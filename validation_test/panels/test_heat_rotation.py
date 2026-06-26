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

        args = SimpleNamespace(q_uniform=None, qsurf_sol=sol, em_vol=vol,
                               qsurf_order=1, surface_label="default")
        # The helper walks BND vertices (filtered by args.surface_label)
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

        args = SimpleNamespace(q_uniform=None, qsurf_sol=sol, em_vol=vol,
                               qsurf_order=1, surface_label="default")
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

        args = SimpleNamespace(q_uniform=None, qsurf_sol=sol, em_vol=vol,
                               qsurf_order=1, surface_label="default")
        q_cf, resample = calc_heat._build_qsurf_cf(wp_mesh, args)

        # Locate the surface vertex nearest to body (+1, 0, 0).
        best_vnr = _surf_vnr_nearest(wp_mesh, (1.0, 0.0, 0.0))

        resample(math.pi / 2.0)
        val_half = q_cf.vec.FV()[best_vnr]
        # At theta=pi/2, body (1,0,0) ↦ world (0,1,0), q_em = x = 0.
        assert val_half == pytest.approx(0.0, abs=0.05)


def test_phi_average_uniform_of_x_is_zero(synthetic_setup):
    """--q-phi-average (the 'uniform' / circumferential-average mode):
    the phi-average of q_em=x around z is 0 everywhere (x = r*cos(phi)
    integrates to 0 over a full turn).  Contrast with the no-rotation
    mode (theta=0 projection) which keeps the LOCAL value ~1 at (1,0,0)
    -- so this test, together with test_rotation_zero_matches_identity,
    pins both modes the panel offers (uniform vs no-rotation)."""
    import calc_heat
    from types import SimpleNamespace
    import tempfile

    em_mesh, gf_q_em, wp_mesh = synthetic_setup

    with tempfile.TemporaryDirectory() as td:
        sol = os.path.join(td, "q.sol")
        vol = os.path.join(td, "wp.vol")
        gf_q_em.Save(sol)
        em_mesh.ngmesh.Save(vol)

        args = SimpleNamespace(q_uniform=None, qsurf_sol=sol, em_vol=vol,
                               qsurf_order=1, surface_label="default",
                               rotation_axis="z",
                               q_phi_average=True, q_phi_average_n=48)
        q_cf, resample = calc_heat._build_qsurf_cf(wp_mesh, args)

        # phi-average produces a STATIC axisymmetric source -> NO resampler
        # (the time loop treats this like a uniform/constant source).
        assert resample is None

        # The face-centre vertex near (1,0,0): rotating it around z sweeps
        # the unit circle (cos t, sin t, 0) -- all inside the cube -- and
        # the mean of q_em=x=cos t over a full turn is 0.
        best_vnr = _surf_vnr_nearest(wp_mesh, (1.0, 0.0, 0.0))
        val_phi_avg = q_cf.vec.FV()[best_vnr]
        assert val_phi_avg == pytest.approx(0.0, abs=0.05)


def test_phi_average_requires_spatial_not_uniform(synthetic_setup):
    """--q-phi-average + --q-uniform is contradictory (a constant is
    already azimuthally uniform) and must fail loud, per No-Fallback."""
    import calc_heat
    from types import SimpleNamespace

    _em, _gf, wp_mesh = synthetic_setup
    args = SimpleNamespace(q_uniform=1.0e6, qsurf_sol="", em_vol="",
                           qsurf_order=1, surface_label="default",
                           rotation_axis="z",
                           q_phi_average=True, q_phi_average_n=48)
    with pytest.raises(ValueError, match="q-phi-average"):
        calc_heat._build_qsurf_cf(wp_mesh, args)
