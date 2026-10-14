"""Contracts of radia.ih_thermal_post.case_depth.

The depth probe must never report its own length as a depth: a ray that
leaves the material still hot is "through", one cut short by the probe is
"beyond_span", and both are distinguished from a real crossing.
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

from radia import ih_thermal_post


def _slab(thickness=0.02, maxh=0.004):
    from netgen.occ import Box, OCCGeometry, Pnt, X
    from ngsolve import Mesh

    solid = Box(Pnt(0, 0, 0), Pnt(thickness, 0.02, 0.02))
    solid.faces.name = "side"
    solid.faces.Min(X).name = "heated"
    return Mesh(OCCGeometry(solid).GenerateMesh(maxh=maxh))


def _field(mesh, expr, order=2):
    from ngsolve import GridFunction, H1
    gf = GridFunction(H1(mesh, order=order))
    gf.Set(expr)
    return gf


def test_linear_profile_gives_the_analytic_depth():
    from ngsolve import x
    mesh = _slab()
    gf = _field(mesh, 1000.0 - 1.0e5 * x)            # 850 C at x = 1.5 mm
    d = ih_thermal_post.case_depth(
        mesh, gf, 850.0, origins=[[0.0, 0.01, 0.01]],
        normals=[[1.0, 0.0, 0.0]], span=0.01)
    assert d["status"] == ["ok"]
    assert d["depth_m"][0] == pytest.approx(1.5e-3, abs=1e-6)


def test_status_distinguishes_every_non_crossing():
    from ngsolve import x
    thin = _slab(thickness=0.004, maxh=0.002)
    hot = _field(thin, 1000.0 + 0.0 * x)
    through = ih_thermal_post.case_depth(
        thin, hot, 850.0, origins=[[0.0, 0.01, 0.01]],
        normals=[[1.0, 0.0, 0.0]], span=0.01)
    assert through["status"] == ["through"]
    assert through["depth_m"][0] == pytest.approx(0.004, abs=5e-4)

    mesh = _slab()
    short = ih_thermal_post.case_depth(
        mesh, _field(mesh, 1000.0 + 0.0 * x), 850.0,
        origins=[[0.0, 0.01, 0.01]], normals=[[1.0, 0.0, 0.0]], span=0.005)
    assert short["status"] == ["beyond_span"]

    cold = ih_thermal_post.case_depth(
        mesh, _field(mesh, 500.0 + 0.0 * x), 850.0,
        origins=[[0.0, 0.01, 0.01]], normals=[[1.0, 0.0, 0.0]], span=0.005)
    assert cold["status"] == ["not_reached"] and cold["depth_m"] == [0.0]


def test_auto_stations_point_inward():
    from ngsolve import x
    mesh = _slab()
    gf = _field(mesh, 1000.0 - 1.0e5 * x)
    d = ih_thermal_post.case_depth(mesh, gf, 850.0,
                                   boundary_names=["heated"], span=0.01)
    normals = np.asarray(d["normals"])
    assert np.allclose(normals, [1.0, 0.0, 0.0], atol=1e-9)
    ok = np.asarray(d["status"]) == "ok"
    assert ok.all()
    np.testing.assert_allclose(np.asarray(d["depth_m"]), 1.5e-3, atol=1e-6)


def test_axisymmetric_mesh_depth_from_the_outer_surface():
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh, x

    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (0.02, 0.01), bcs=("b", "outer", "t", "axis"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.002))
    gf = _field(mesh, 1000.0 - 1.0e5 * (0.02 - x))
    d = ih_thermal_post.case_depth(mesh, gf, 850.0,
                                   boundary_names=["outer"], span=0.01)
    np.testing.assert_allclose(np.asarray(d["depth_m"]), 1.5e-3, atol=1e-6)


def test_meridian_stations_report_the_azimuthal_range():
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import Mesh, sqrt, x, y

    mesh = Mesh(OCCGeometry(Cylinder(Axes(Pnt(0, 0, 0), Z), r=0.02, h=0.01))
                .GenerateMesh(maxh=0.002))
    r = sqrt(x * x + y * y)
    # depth = 1.5 mm * (1 + 0.25 sin(phi)): 1.125 mm at 270 deg, 1.875 at 90
    gf = _field(mesh, 1000.0 - 1.0e5 * (0.02 - r) / (1.0 + 0.25 * y / r))
    d = ih_thermal_post.case_depth(
        mesh, gf, 850.0, origins=[[0.02, 0.005]], normals=[[-1.0, 0.0]],
        span=0.008, n_phi=8)
    assert d["status_counts"]["ok"] == 8
    assert d["depth_min_m"][0] == pytest.approx(1.125e-3, abs=2e-5)
    assert d["depth_max_m"][0] == pytest.approx(1.875e-3, abs=2e-5)


def test_span_is_required_and_positive():
    from ngsolve import x
    mesh = _slab()
    with pytest.raises(TypeError, match="span"):
        ih_thermal_post.case_depth(mesh, _field(mesh, x), 850.0,
                                   boundary_names=["heated"])
    with pytest.raises(ValueError, match="span"):
        ih_thermal_post.case_depth(mesh, _field(mesh, x), 850.0, span=0.0,
                                   boundary_names=["heated"])


def test_command_line_reads_the_sidecar(tmp_path):
    import json
    from radia import ih_thermal
    from ngsolve import x

    mesh = _slab()
    vol = tmp_path / "slab.vol"
    mesh.ngmesh.Save(str(vol))
    from ngsolve import Mesh
    mesh = Mesh(str(vol))
    gf = _field(mesh, 1000.0 - 1.0e5 * x)
    sol = tmp_path / "T.sol"
    gf.Save(str(sol))
    ih_thermal.write_field_sidecar(
        str(sol), mesh_path=str(vol), mesh=mesh, fes_order=2,
        quantity=ih_thermal.TEMPERATURE_QUANTITY,
        unit=ih_thermal.TEMPERATURE_UNIT, extra={"geometry": "3d"})
    out = tmp_path / "depth.json"
    csv_out = tmp_path / "depth.csv"
    assert ih_thermal_post.main([
        "depth", "--temperature", str(sol), "--threshold", "850",
        "--boundaries", "heated", "--span", "0.01",
        "--output", str(out), "--csv", str(csv_out)]) == 0
    result = json.loads(out.read_text(encoding="utf-8"))
    assert result["field"]["provenance"] == "sidecar-verified"
    np.testing.assert_allclose(result["depth_m"], 1.5e-3, atol=1e-6)
    assert csv_out.read_text(encoding="utf-8").startswith("s,")

def test_a_2d_mesh_takes_meridian_stations_only():
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh, x

    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (0.02, 0.01), bcs=("b", "outer", "t", "axis"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.002))
    gf = _field(mesh, 1000.0 - 1.0e5 * (0.02 - x))
    d = ih_thermal_post.case_depth(mesh, gf, 850.0, span=0.01,
                                   origins=[[0.02, 0.005]], normals=[[-1, 0]])
    assert d["depth_m"][0] == pytest.approx(1.5e-3, abs=1e-6)
    with pytest.raises(ValueError, match="do not fit a 2D mesh"):
        ih_thermal_post.case_depth(mesh, gf, 850.0, span=0.01,
                                   origins=[[0.02, 0.0, 0.005]],
                                   normals=[[-1, 0, 0]])


def test_a_station_off_the_surface_is_refused():
    from ngsolve import x
    mesh = _slab()
    with pytest.raises(ValueError, match="farther"):
        ih_thermal_post.case_depth(
            mesh, _field(mesh, 1000.0 - 1.0e5 * x), 850.0, span=0.01,
            origins=[[0.005, 0.005, 0.005]], normals=[[1, 0, 0]])


@pytest.mark.parametrize("header, ok", [
    ("x_m,y_m,z_m,nx,ny,nz", True),
    ("s_mm,r_mm,z_mm,n_r,n_z", True),
    ("x,y,z,nx,ny,nz", False),              # no unit: not guessed
    ("x_m,y_mm,z_m,nx,ny,nz", False),       # mixed units
    ("r_m,z_m,nr,nz", False),               # alias of n_r, n_z
])
def test_stations_csv_header_is_exact(tmp_path, header, ok):
    n = len(header.split(","))
    csv = tmp_path / "st.csv"
    csv.write_text(header + "\n" + ",".join(["1"] * n) + "\n", encoding="utf-8")
    if ok:
        o, nrm, s = ih_thermal_post._read_stations(str(csv))
        scale = 1e-3 if "_mm" in header else 1.0
        assert np.allclose(o, scale) and np.allclose(nrm, 1.0)
        assert (s is not None) == header.startswith("s_")
    else:
        with pytest.raises(ValueError):
            ih_thermal_post._read_stations(str(csv))
