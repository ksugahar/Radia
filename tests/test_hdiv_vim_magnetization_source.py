"""Fast contracts for prescribed 3D magnetization in an independent HDiv space."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
import radia as rad
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve.meshes import MakeStructured3DMesh

from radia import vim
from radia.vim._magnetization_source import _field_coefficient_algorithm
from tests._ngsolve_2606 import curve_mesh


def _box_mesh(x0, x1, maxh=0.28):
    geometry = CSGeometry()
    geometry.Add(OrthoBrick(Pnt(x0, -0.2, -0.2), Pnt(x1, 0.2, 0.2)).mat("body"))
    return ng.Mesh(geometry.GenerateMesh(maxh=maxh))


def test_large_fixed_magnetization_source_uses_native_tree_coefficient_path():
    assert _field_coefficient_algorithm({"source_count": 511}, None) == "direct"
    assert _field_coefficient_algorithm({"source_count": 512}, None) == "tree"
    assert _field_coefficient_algorithm({"source_count": 10000}, "direct") == "direct"
    assert _field_coefficient_algorithm({"source_count": 1}, "tree") == "tree"
    with pytest.raises(ValueError, match="field_cf_algorithm"):
        _field_coefficient_algorithm({"source_count": 1}, "auto")


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
        ("tet-bdm2", _box_mesh(-0.3, 0.3, maxh=0.35), 2, None, "analytic-tet-bdm2"),
        ("tet-curved", _box_mesh(-0.3, 0.3, maxh=0.35), 1, 2, "curved-element-exact-bdm1"),
        ("tet-curved-bdm2", _box_mesh(-0.3, 0.3, maxh=0.35), 2, 2,
         "curved-element-exact-bdm2"),
        ("hex-bdm1", MakeStructured3DMesh(
            hexes=True, nx=1, ny=1, nz=1, mapping=mapping), 1, None, "analytic-tet-bdm1"),
        ("hex-curved-bdm2", MakeStructured3DMesh(
            hexes=True, nx=1, ny=1, nz=1, mapping=mapping), 2, 2, "analytic-tet-bdm2"),
        ("wedge-bdm1", MakeStructured3DMesh(
            prism=True, nx=1, ny=1, nz=1, mapping=mapping), 1, None, "element-cloud-bdm1"),
        ("wedge-curved-bdm2", MakeStructured3DMesh(
            prism=True, nx=1, ny=1, nz=1, mapping=mapping), 2, 2, "element-cloud-bdm2"),
    ]
    prescribed = ng.CoefficientFunction((1.0e5, 2.0e5, 3.0e5))
    for name, mesh, order, curve_order, expected_kind in cases:
        with ng.TaskManager():
            source = vim.MagnetizationSource(
                mesh, prescribed, order=order, curve_order=curve_order)
        assert source.stats["hmatrix_built"] is False, name
        assert source.stats["field_evaluator"]["source_kind"] == expected_kind, name
        assert source.stats["curve_order"] == curve_order, name
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


@pytest.mark.parametrize("kind, mesh_options", [
    ("hex", {"hexes": True}),
    ("wedge", {"prism": True}),
])
def test_curved_hex_wedge_rt2_field_matches_independent_ngsolve_boundary_integral(
        kind, mesh_options):
    mesh = MakeStructured3DMesh(
        nx=1, ny=1, nz=1, **mesh_options,
        mapping=lambda x, y, z: (x-0.5, y-0.5, z-0.5))
    expected_vertices = 8 if kind == "hex" else 6
    if {len(element.vertices) for element in mesh.Elements(ng.VOL)} != {expected_vertices}:
        pytest.skip(f"structured {kind} generator returned a different topology")
    curve_mesh(mesh, 2)
    with ng.TaskManager():
        source = vim.MagnetizationSource(
            mesh, (1.0e5, 2.0e5, 3.0e5), order=2, curve_order=2)

    points = np.asarray([[2.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    normal = ng.specialcf.normal(3)
    position = ng.CoefficientFunction((ng.x, ng.y, ng.z))
    sigma = ng.InnerProduct(source.magnetization, normal)
    reference = []
    with ng.TaskManager():
        for point in points:
            delta = ng.CoefficientFunction(tuple(point)) - position
            radius = ng.sqrt(ng.InnerProduct(delta, delta))
            reference.append(ng.Integrate(
                sigma*delta/(4.0*np.pi*radius**3), mesh, ng.BND, order=12))

    reference = np.asarray(reference, dtype=float)
    actual = source.Field(points, "direct")
    relative = np.linalg.norm(actual-reference, axis=1) / np.linalg.norm(reference, axis=1)
    assert relative.max() < 1.0e-5
