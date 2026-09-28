"""Contracts of radia.ih_thermal: field artifacts and heat-source transfer.

These guard the failure modes that silently corrupted induction-heating
results before: a .sol read with the wrong mesh or order, a mesh modified
after loading, flux points that miss the source surface being set to zero,
and an "axisymmetric" average that still depends on the azimuth.
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

from radia import ih_thermal


@pytest.fixture(scope="module")
def cylinder(tmp_path_factory):
    """Cylinder r=1, h=2 about z with a P1 field q = 1 + z^2 saved on it."""
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import GridFunction, H1, Mesh, z

    d = tmp_path_factory.mktemp("cyl")
    solid = Cylinder(Axes(Pnt(0, 0, -1), Z), r=1.0, h=2.0)
    solid.faces.name = "side"
    solid.faces.Max(Z).name = "top"
    solid.faces.Min(Z).name = "bottom"
    vol = str(d / "cyl.vol")
    OCCGeometry(solid).GenerateMesh(maxh=0.3).Save(vol)
    mesh = Mesh(vol)
    gf = GridFunction(H1(mesh, order=1))
    gf.Set(1.0 + z * z)
    sol = str(d / "q.sol")
    gf.Save(sol)
    return {"dir": d, "vol": vol, "sol": sol, "mesh": mesh, "gf": gf}


# ---------------------------------------------------------------- artifacts

def test_sidecar_round_trip_and_load(cylinder):
    ih_thermal.write_field_sidecar(
        cylinder["sol"], mesh_path=cylinder["vol"], mesh=cylinder["mesh"],
        fes_order=1, quantity=ih_thermal.QSURF_QUANTITY,
        unit=ih_thermal.QSURF_UNIT, boundaries=["side"],
        extra={"P_wp_W": 1.0})
    mesh, gf, audit = ih_thermal.load_field(
        cylinder["sol"], quantity=ih_thermal.QSURF_QUANTITY)
    assert audit["provenance"] == "sidecar-verified"
    assert audit["mesh_curve_order"] == mesh.GetCurveOrder()
    np.testing.assert_array_equal(gf.vec.FV().NumPy(),
                                  cylinder["gf"].vec.FV().NumPy())


def test_wrong_order_is_rejected_by_size(cylinder):
    with pytest.raises(ValueError, match="do not belong together"):
        ih_thermal.verify_field_pair(cylinder["sol"], cylinder["vol"],
                                     cylinder["mesh"], 2,
                                     quantity=ih_thermal.QSURF_QUANTITY)


def test_modified_field_is_rejected(cylinder, tmp_path):
    sol = tmp_path / "q.sol"
    data = np.fromfile(cylinder["sol"], dtype="<f8")
    data.tofile(sol)
    ih_thermal.write_field_sidecar(
        str(sol), mesh_path=cylinder["vol"], mesh=cylinder["mesh"],
        fes_order=1, quantity="q", unit="W/m^2")
    (2.0 * data).tofile(sol)
    with pytest.raises(ValueError, match="changed after its sidecar"):
        ih_thermal.verify_field_pair(str(sol), cylinder["vol"],
                                     cylinder["mesh"], 1, quantity="q")


def test_a_different_mesh_file_is_rejected(cylinder, tmp_path):
    other = tmp_path / "other.vol"
    with open(cylinder["vol"], "rb") as src:
        other.write_bytes(src.read() + b"\n")
    sol = tmp_path / "q.sol"
    np.fromfile(cylinder["sol"], dtype="<f8").tofile(sol)
    ih_thermal.write_field_sidecar(
        str(sol), mesh_path=cylinder["vol"], mesh=cylinder["mesh"],
        fes_order=1, quantity="q", unit="W/m^2")
    with pytest.raises(ValueError, match="mesh file digest"):
        ih_thermal.verify_field_pair(str(sol), str(other), cylinder["mesh"], 1,
                                     quantity="q")


def test_scaled_artifact_carries_its_provenance(cylinder, tmp_path):
    src = tmp_path / "q.sol"
    np.fromfile(cylinder["sol"], dtype="<f8").tofile(src)
    ih_thermal.write_field_sidecar(
        str(src), mesh_path=cylinder["vol"], mesh=cylinder["mesh"],
        fes_order=1, quantity=ih_thermal.QSURF_QUANTITY,
        unit=ih_thermal.QSURF_UNIT, boundaries=["side"],
        extra={"P_wp_W": 10.0})
    out = tmp_path / "q_scaled.sol"
    ih_thermal.scale_field_artifact(str(src), str(out), 0.25)
    rec = ih_thermal.read_field_sidecar(str(out))
    assert rec["P_wp_W"] == pytest.approx(2.5)
    assert rec["scaled_from"]["factor"] == 0.25
    np.testing.assert_allclose(np.fromfile(out, dtype="<f8"),
                               0.25 * np.fromfile(src, dtype="<f8"))


# ---------------------------------------------------------------- transfer

def test_surface_field_rejects_points_off_the_source(cylinder):
    field = ih_thermal.SurfaceP1Field.from_gridfunction(
        cylinder["mesh"], cylinder["gf"], ["side"])
    on = np.array([[1.0, 0.0, 0.3], [0.0, -1.0, -0.5]])
    vals, dist = field.evaluate(on, max_distance=0.05)
    np.testing.assert_allclose(vals, 1.0 + on[:, 2] ** 2, atol=0.02)
    with pytest.raises(ValueError, match="no zero-flux fallback"):
        field.evaluate(np.array([[1.5, 0.0, 0.0]]), max_distance=0.05)


def test_profile_conserves_power_and_is_azimuth_free(cylinder):
    from ngsolve import GridFunction, H1, x

    gf = GridFunction(H1(cylinder["mesh"], order=1))
    gf.Set(2.0 + x)                     # azimuthal part averages out
    prof = ih_thermal.AxisymmetricSurfaceProfile.from_gridfunction(
        cylinder["mesh"], gf, ["side"])
    assert prof.power == pytest.approx(prof.source_power, rel=1e-12)
    phis = np.linspace(0, 2 * math.pi, 7, endpoint=False)
    ring = np.c_[np.cos(phis), np.sin(phis), np.full(7, 0.2)]
    vals, _ = prof.evaluate_xyz(ring)
    assert np.ptp(vals) == 0.0          # identical at every azimuth
    assert vals[0] == pytest.approx(2.0, abs=0.02)


def test_profile_rejects_a_body_that_is_not_of_revolution(tmp_path):
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import GridFunction, H1, Mesh

    mesh = Mesh(OCCGeometry(Box(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
                .GenerateMesh(maxh=0.4))
    gf = GridFunction(H1(mesh, order=1))
    gf.vec[:] = 1.0
    with pytest.raises(ValueError, match="not a body of revolution"):
        ih_thermal.AxisymmetricSurfaceProfile.from_gridfunction(
            mesh, gf, mesh.GetBoundaries())


def test_rotation_about_each_axis_is_right_handed():
    p = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    np.testing.assert_allclose(
        ih_thermal.rotate_about_axis(p[:1], math.pi / 2, "z"), [[0, 1, 0]],
        atol=1e-15)
    np.testing.assert_allclose(
        ih_thermal.rotate_about_axis(p[1:2], math.pi / 2, "x"), [[0, 0, 1]],
        atol=1e-15)
    np.testing.assert_allclose(
        ih_thermal.rotate_about_axis(p[2:], math.pi / 2, "y"), [[1, 0, 0]],
        atol=1e-15)


def test_power_balance_refuses_a_lossy_transfer():
    ok = ih_thermal.power_balance(1.001, 1.0, 0.02, what="t")
    assert ok["relative_error"] == pytest.approx(1e-3)
    with pytest.raises(ValueError, match="refusing to solve"):
        ih_thermal.power_balance(0.5, 1.0, 0.02, what="t")
    # a signed test field with no net power is judged against |q|
    assert ih_thermal.power_balance(1e-9, 0.0, 0.02, what="t",
                                    scale=1.0)["relative_error"] < 1e-8


def test_em_boundaries_come_from_the_request_or_the_sidecar(cylinder):
    """Nothing is borrowed from the thermal mesh's boundary names."""
    from ngsolve import GridFunction, H1, Mesh

    mesh = Mesh(cylinder["vol"])
    gf = GridFunction(H1(mesh, order=1))
    gf.Load(cylinder["sol"])
    with pytest.raises(ValueError, match="state the heated"):
        ih_thermal.resolve_em_heat_boundaries(mesh, gf, sidecar={})
    names, how = ih_thermal.resolve_em_heat_boundaries(
        mesh, gf, requested="side|top")
    assert names == ["side", "top"] and how == "explicit"
    names, how = ih_thermal.resolve_em_heat_boundaries(
        mesh, gf, sidecar={"boundaries": ["side"]})
    assert names == ["side"] and how == "sidecar"


def test_a_field_without_a_sidecar_is_refused(cylinder, tmp_path):
    """No byte-count-only acceptance: the pairing must be stated."""
    sol = tmp_path / "bare.sol"
    np.fromfile(cylinder["sol"], dtype="<f8").tofile(sol)
    with pytest.raises(ValueError, match="no .sol.json sidecar|no .json"):
        ih_thermal.verify_field_pair(str(sol), cylinder["vol"],
                                     cylinder["mesh"], 1,
                                     quantity=ih_thermal.QSURF_QUANTITY)
    with pytest.raises(ValueError, match="sidecar"):
        ih_thermal.load_field(str(sol), quantity=ih_thermal.QSURF_QUANTITY)
    # the explicit adoption path
    assert ih_thermal.main([
        "sidecar", "--sol", str(sol), "--mesh", cylinder["vol"],
        "--order", "1", "--quantity", ih_thermal.QSURF_QUANTITY,
        "--boundaries", "side", "--p-wp", "1.0"]) == 0
    assert ih_thermal.verify_field_pair(
        str(sol), cylinder["vol"], cylinder["mesh"], 1,
        quantity=ih_thermal.QSURF_QUANTITY)["provenance"] == \
        "sidecar-verified"


def test_a_field_of_another_quantity_is_refused(cylinder):
    with pytest.raises(ValueError, match="quantity"):
        ih_thermal.verify_field_pair(cylinder["sol"], cylinder["vol"],
                                     cylinder["mesh"], 1,
                                     quantity=ih_thermal.HT_QUANTITY)


def test_curved_em_mesh_passes_the_sidecar_power_check(tmp_path):
    """The EM power in the sidecar is the curved-geometry integral; the
    check must compare like with like (flat triangles differ by ~0.3 %)."""
    from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, Z
    from ngsolve import BND, CF, GridFunction, H1, Integrate, Mesh

    solid = Cylinder(Axes(Pnt(0, 0, 0), Z), r=0.01, h=0.02)
    solid.faces.name = "side"
    solid.faces.Max(Z).name = "top"
    solid.faces.Min(Z).name = "bottom"
    live = Mesh(OCCGeometry(solid).GenerateMesh(maxh=0.004))
    live.Curve(2)                       # curve while the geometry is live
    vol = str(tmp_path / "curved.vol")
    live.ngmesh.Save(vol)
    mesh = Mesh(vol)
    assert mesh.GetCurveOrder() == 2
    gf = GridFunction(H1(mesh, order=1))
    gf.Set(CF(1.0e6), definedon=mesh.Boundaries("side"))
    sol = str(tmp_path / "q.sol")
    gf.Save(sol)
    P = float(Integrate(gf, mesh, BND, definedon=mesh.Boundaries("side")))
    ih_thermal.write_field_sidecar(
        sol, mesh_path=vol, mesh=mesh, fes_order=1,
        quantity=ih_thermal.QSURF_QUANTITY, unit=ih_thermal.QSURF_UNIT,
        boundaries=["side"], extra={"P_wp_W": P})
    src = ih_thermal.EMHeatSource(sol, vol)
    assert src.sidecar_power["relative_error"] == pytest.approx(0.0, abs=1e-9)
    assert abs(src.power_flat / src.power - 1.0) > 1e-3

def test_axis_split_is_exact_for_p1_fields():
    """A triangle pierced by the axis is split at the piercing point; the
    P1 integral is unchanged and the new vertex carries the P1 value."""
    pts = np.array([[1.0, 0.0, 0.0], [-0.5, 0.8, 0.0], [-0.5, -0.9, 0.0],
                    [2.0, 2.0, 0.0]])
    tris = np.array([[0, 1, 2], [0, 3, 1]])
    vals = np.array([1.0, 2.0, 4.0, 3.0])
    p2, t2, v2 = ih_thermal._split_at_axis(pts, tris, vals, "z")
    assert len(t2) == 4 and len(p2) == 5          # only the pierced one split
    assert np.allclose(p2[4, :2], 0.0)
    assert ih_thermal.p1_surface_power(p2, t2, v2) == pytest.approx(
        ih_thermal.p1_surface_power(pts, tris, vals), rel=1e-14)
    # the P1 value at the origin of the first triangle
    lam = np.linalg.solve(np.r_[pts[tris[0], :2].T, np.ones((1, 3))],
                          [0.0, 0.0, 1.0])
    assert v2[4] == pytest.approx(lam @ vals[tris[0]], rel=1e-14)
