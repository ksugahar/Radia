"""Contracts of radia.ih_thermal_post.thermal_exposure.

Overheating must be visible in the result: a molten skin that no vertex
falls inside, or a hot spot at one azimuth of a nominally axisymmetric
part, was reported as "0 % above melting" before.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "radia"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))

from radia import ih_thermal_post


def _cube(maxh=0.2):
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import Mesh
    return Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1)))
                .GenerateMesh(maxh=maxh))


def test_volume_above_threshold_matches_the_analytic_value():
    from ngsolve import GridFunction, H1, x

    mesh = _cube()
    gf = GridFunction(H1(mesh, order=2))
    gf.Set(1000.0 * x)
    ex = ih_thermal_post.thermal_exposure(mesh, gf, [600.0, 2000.0],
                                          ring_samples=0)
    above = {e["T_C"]: e for e in ex["thresholds"]}
    assert above[600.0]["volume_m3"] == pytest.approx(0.4, abs=2e-3)
    assert above[600.0]["bbox_min_m"][0] == pytest.approx(0.6, abs=0.05)
    assert above[2000.0]["volume_m3"] == 0.0
    assert ex["T_max_C"] == pytest.approx(1000.0, rel=1e-9)
    assert ex["T_max_location_m"][0] == pytest.approx(1.0)


def test_a_hot_layer_on_a_coarse_mesh_is_counted():
    """The volume comes from integration points, not vertices, so a layer
    narrower than the elements is still measured."""
    from ngsolve import GridFunction, H1, x

    mesh = _cube(maxh=0.5)
    gf = GridFunction(H1(mesh, order=2))
    gf.Set(1000.0 - 16000.0 * (x - 0.25) ** 2)
    ex = ih_thermal_post.thermal_exposure(mesh, gf, [900.0], ring_samples=0)
    exact = 2.0 * math.sqrt(100.0 / 16000.0)     # width of the T > 900 layer
    assert ex["thresholds"][0]["volume_m3"] == pytest.approx(exact, rel=0.1)
    vertex_max = max(gf(mesh(*v.point)) for v in mesh.vertices)
    assert ex["T_max_C"] >= vertex_max
    assert ex["T_max_C"] == pytest.approx(1000.0, abs=5.0)


def test_axisymmetric_volume_is_revolved():
    """T = 1000 r on a (r, z) rectangle: volume above 500 is pi (1 - 0.25)."""
    from netgen.geom2d import SplineGeometry
    from ngsolve import GridFunction, H1, Mesh, x

    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (1, 1), bcs=("b", "o", "t", "axis"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.1))
    gf = GridFunction(H1(mesh, order=2))
    gf.Set(1000.0 * x)
    ex = ih_thermal_post.thermal_exposure(mesh, gf, [500.0],
                                          axisymmetric=True)
    assert ex["domain_measure_m3"] == pytest.approx(math.pi, rel=1e-3)
    assert ex["thresholds"][0]["volume_m3"] == pytest.approx(
        math.pi * 0.75, rel=5e-3)
    assert "hottest_ring" not in ex


def test_hottest_ring_exposes_an_azimuthal_hot_spot():
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import GridFunction, H1, Mesh, x, z

    mesh = Mesh(OCCGeometry(Cylinder(Axes(Pnt(0, 0, 0), Z), r=1, h=1))
                .GenerateMesh(maxh=0.2))
    axisym = GridFunction(H1(mesh, order=1))
    axisym.Set(100.0 * z)
    ring = ih_thermal_post.thermal_exposure(mesh, axisym, [])["hottest_ring"]
    assert ring["spread_C"] < 1e-6 * 100
    hot = GridFunction(H1(mesh, order=1))
    hot.Set(100.0 * z + 50.0 * x)
    ring = ih_thermal_post.thermal_exposure(mesh, hot, [])["hottest_ring"]
    assert ring["n_inside"] == ring["n_requested"]
    assert ring["spread_C"] > 90.0


def test_limit_check_is_an_optimiser_constraint():
    from ngsolve import GridFunction, H1, x

    mesh = _cube()
    gf = GridFunction(H1(mesh, order=1))
    gf.Set(1000.0 * x)
    ex = ih_thermal_post.thermal_exposure(mesh, gf, [900.0], ring_samples=0)
    rec = ih_thermal_post.limit_check(ex, 900.0)
    assert rec["exceeded"] and rec["excess_C"] == pytest.approx(100.0)
    assert rec["volume_above_m3"] == pytest.approx(0.1, abs=3e-3)
    with pytest.raises(ValueError):
        ih_thermal_post.limit_check(ex, 800.0)


def test_heat_solver_reports_exposure_and_limit():
    import calc_heat
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import Mesh

    solid = Box(Pnt(0, 0, 0), Pnt(0.02, 0.02, 0.02))
    solid.faces.name = "side"
    solid.faces.Max((0, 0, 1)).name = "heated"
    mesh = Mesh(OCCGeometry(solid).GenerateMesh(maxh=0.005))
    result = calc_heat.solve_heat(
        "<in-memory-block>", material="steel", h_conv=0.0,
        heat_flux_boundaries="heated", q_uniform=5.0e6, dt=0.1, t_end=0.5,
        t_initial=20.0, fes_order=2, _wp_mesh=mesh, _write_solution=False,
        exposure_thresholds=[200.0], temperature_limit=400.0)
    assert "error" not in result, result
    ex = result["thermal_exposure"]
    assert [e["T_C"] for e in ex["thresholds"]] == [200.0, 400.0]
    assert result["temperature_limit"]["limit_C"] == 400.0
    assert len(result["T_max_history_C"]) == len(result["t_history_s"])
    assert np.all(np.diff(result["T_max_history_C"]) > 0)
    assert ex["T_max_C"] == pytest.approx(result["T_max_C"], rel=2e-3)
