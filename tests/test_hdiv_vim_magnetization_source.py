"""Fast contracts for prescribed 3D magnetization in an independent HDiv space."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
import radia as rad
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve.meshes import MakeStructured3DMesh

from radia import vim


def _box_mesh(x0, x1, maxh=0.28):
    geometry = CSGeometry()
    geometry.Add(OrthoBrick(Pnt(x0, -0.2, -0.2), Pnt(x1, 0.2, 0.2)).mat("body"))
    return ng.Mesh(geometry.GenerateMesh(maxh=maxh))


def test_prescribed_magnetization_is_l2_projected_and_native_field_matches_direct():
    pm_mesh = _box_mesh(-0.5, -0.1)
    iron_mesh = _box_mesh(0.1, 0.5)
    prescribed = ng.CoefficientFunction((
        1.0e5*(1.0+ng.x), 2.0e5*ng.y, 3.0e5*(1.0+ng.z)))

    with ng.TaskManager():
        source = vim.MagnetizationSource(pm_mesh, prescribed, order=1)
        error2 = ng.Integrate(
            ng.InnerProduct(source.magnetization-prescribed,
                            source.magnetization-prescribed), pm_mesh)

    assert source.stats["hmatrix_built"] is False
    assert source.stats["permanent_magnet_model"] == "fixed-given"
    assert source.stats["permanent_magnet_level"] == 1
    assert source.stats["projection_relative_residual"] < 1.0e-12
    assert np.sqrt(error2) < 1.0e-8

    points = np.array([[0.2, 0.0, 0.0], [0.3, 0.05, 0.02]])
    direct = source.Field(points, algorithm="direct")
    with ng.TaskManager():
        native_cf = np.asarray(
            source.field_cf(iron_mesh(points[:, 0], points[:, 1], points[:, 2])), float)
    assert np.allclose(native_cf, direct, rtol=5.0e-14, atol=1.0e-8)


def test_prescribed_source_coupling_equals_explicit_native_field_and_stays_fixed():
    pm_mesh = _box_mesh(-0.5, -0.1)
    iron_mesh = _box_mesh(0.1, 0.5)
    with ng.TaskManager():
        source = vim.MagnetizationSource(pm_mesh, (0.0, 0.0, 8.0e5), order=1)
    original = source._coefficients.copy()

    with ng.TaskManager():
        coupled = vim.Solve(
            iron_mesh, mu_r=20.0, magnetization_sources=source,
            order=1, tol=1.0e-9)
        explicit = vim.Solve(
            iron_mesh, mu_r=20.0, H_ext=source.field_cf,
            order=1, tol=1.0e-9)

    assert coupled["magnetization_source_count"] == 1
    assert coupled["_magnetization_sources"] == (source,)
    assert np.linalg.norm(coupled["M_avg"]) > 1.0
    assert np.allclose(coupled["M_avg"], explicit["M_avg"], rtol=2.0e-13, atol=1.0e-8)
    assert np.array_equal(original, source._coefficients)
    assert np.array_equal(original, source.magnetization.vec.FV().NumPy())


def test_prescribed_sources_superpose_and_same_mesh_is_rejected():
    pm_mesh = _box_mesh(-0.5, -0.1)
    iron_mesh = _box_mesh(0.1, 0.5)
    with ng.TaskManager():
        source_z = vim.MagnetizationSource(
            pm_mesh, ng.CoefficientFunction((0.0, 0.0, 5.0e5)), order=1)
        source_y = vim.MagnetizationSource(
            pm_mesh, ng.CoefficientFunction((0.0, 2.0e5, 0.0)), order=1)

    points = np.array([[0.2, 0.0, 0.0], [0.3, 0.05, 0.02]])
    expected = source_z.Field(points, "direct") + source_y.Field(points, "direct")
    with ng.TaskManager():
        actual = np.asarray(
            (source_z.field_cf+source_y.field_cf)(
                iron_mesh(points[:, 0], points[:, 1], points[:, 2])), float)
    assert np.allclose(actual, expected, rtol=5.0e-14, atol=1.0e-8)

    with ng.TaskManager(), pytest.raises(ValueError, match="separate mesh"):
        vim.Solve(pm_mesh, mu_r=20.0, magnetization_sources=source_z, order=1)


def test_prescribed_source_geometry_only_path_covers_production_3d_elements():
    def mapping(x, y, z):
        return (0.6*x-0.3, 0.4*y-0.2, 0.2*z-0.1)

    cases = [
        ("tet-rt2", _box_mesh(-0.3, 0.3, maxh=0.35), 2, None, "analytic-tet-rt2"),
        ("tet-curved", _box_mesh(-0.3, 0.3, maxh=0.35), 1, 2, "curved-element-exact-rt1"),
        ("hex-rt1", MakeStructured3DMesh(
            hexes=True, nx=1, ny=1, nz=1, mapping=mapping), 1, None, "element-cloud-rt1"),
        ("wedge-rt1", MakeStructured3DMesh(
            prism=True, nx=1, ny=1, nz=1, mapping=mapping), 1, None, "element-cloud-rt1"),
    ]
    prescribed = ng.CoefficientFunction((1.0e5, 2.0e5, 3.0e5))
    for name, mesh, order, curve_order, expected_kind in cases:
        with ng.TaskManager():
            source = vim.MagnetizationSource(
                mesh, prescribed, order=order, curve_order=curve_order)
        assert source.stats["hmatrix_built"] is False, name
        assert source.stats["field_evaluator"]["source_kind"] == expected_kind, name
        assert np.isfinite(source.Field([[1.0, 0.0, 0.0]], "direct")).all(), name


def test_tet_source_near_touching_interface_matches_uniform_cuboid_field():
    pm_mesh = _box_mesh(-0.5, -0.1)
    magnetization = (0.0, 0.0, 8.0e5)
    with ng.TaskManager():
        source = vim.MagnetizationSource(pm_mesh, magnetization, order=1)

    points = np.array([
        [-0.099, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.2, 0.0, 0.0],
        [0.2, 0.1, 0.1],
    ])
    analytic_magnet = vim.magnet_box(
        center=(-0.3, 0.0, 0.0), size=(0.4, 0.4, 0.4), M=magnetization)
    try:
        analytic = np.asarray(rad.Fld(analytic_magnet, "h", points), float)
    finally:
        rad.UtiDel(analytic_magnet)
    projected = source.Field(points, "direct")
    assert np.allclose(projected, analytic, rtol=2.0e-13, atol=1.0e-7)


def test_prescribed_source_image_reconstructs_full_cuboid_to_roundoff():
    half_mesh = _box_mesh(0.0, 0.2, maxh=0.2)
    magnetization = (0.0, 0.0, 8.0e5)
    with ng.TaskManager():
        source = vim.MagnetizationSource(
            half_mesh, magnetization, order=1, image="+x")
    assert source.stats["field_evaluator"]["image_count"] == 1

    points = np.array([[0.3, 0.0, 0.0], [0.4, 0.1, 0.05], [-0.3, 0.0, 0.0]])
    analytic_magnet = vim.magnet_box(
        center=(0.0, 0.0, 0.0), size=(0.4, 0.4, 0.4), M=magnetization)
    try:
        analytic = np.asarray(rad.Fld(analytic_magnet, "h", points), float)
    finally:
        rad.UtiDel(analytic_magnet)
    reduced = source.Field(points, "direct")
    assert np.allclose(reduced, analytic, rtol=2.0e-13, atol=1.0e-7)
