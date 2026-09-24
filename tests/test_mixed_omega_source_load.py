"""Surface-flux source loads for mixed total/reduced Omega.

For a divergence-free coil field ``int_Omega H_s . grad(v) dx`` equals
``oint (H_s . n_out) v ds``.  These tests pin the boundary orientation the
surface form relies on, check that the two discrete loads and the resulting
fields agree for a non-polynomial divergence-free source, and check that the
surface form refuses the cases it cannot represent (a linked source, a
permeability that is not the declared constant, a missing gate).
"""
from __future__ import annotations

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")


def _mesh(maxh):
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X

    reduced = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    reduced.mat("reduced")
    reduced.faces.name = "reduced_outer"
    reduced.faces.Max(X).name = "source_total_interface"
    total = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    total.mat("total")
    total.faces.name = "total_outer"
    total.faces.Min(X).name = "source_total_interface"
    return ng.Mesh(OCCGeometry(Glue([reduced, total])).GenerateMesh(maxh=maxh))


def _divergence_free_source():
    # A gradient of the exterior point-charge potential (curl- and div-free
    # in the box) plus curl((0, 0, f)), which is div-free everywhere and has
    # current only in the reduced half x < 0, as a coil does (C1 at x = 0).
    r_x = ng.x - 2.5
    r2 = r_x * r_x + ng.y * ng.y + ng.z * ng.z
    gradient = ng.CF((r_x / r2**1.5, ng.y / r2**1.5, ng.z / r2**1.5))
    a = ng.IfPos(-ng.x, ng.x**2 * (ng.x + 1.0)**2, 0.0)
    b = (1.0 - ng.y**2)**2
    c = (1.0 - ng.z**2)**2
    da_dx = ng.IfPos(-ng.x, 2.0 * ng.x * (ng.x + 1.0) * (2.0 * ng.x + 1.0), 0.0)
    db_dy = -4.0 * ng.y * (1.0 - ng.y**2)
    return gradient + ng.CF((a * db_dy * c, -da_dx * b * c, 0.0))


def test_boundary_normal_points_from_domin_to_domout():
    from radia.kelvin_solver import material_set_boundary_signs

    mesh = _mesh(0.5)
    position = ng.CF((ng.x, ng.y, ng.z))
    for material, volume in (("reduced", 4.0), ("total", 4.0)):
        signs = material_set_boundary_signs(mesh, (material,))
        flux = ng.Integrate(
            mesh.BoundaryCF(signs, default=0.0) * ng.InnerProduct(position, ng.specialcf.normal(3))
            * ng.ds(definedon=mesh.Boundaries("|".join(signs))), mesh)
        assert flux / 3.0 == pytest.approx(volume, rel=1e-12)
    assert material_set_boundary_signs(mesh, ("reduced",))["source_total_interface"] == 1.0
    assert material_set_boundary_signs(mesh, ("total",))["source_total_interface"] == -1.0


def test_shared_label_cannot_carry_one_orientation():
    from radia.kelvin_solver import material_set_boundary_signs
    from netgen.occ import Box, Glue, OCCGeometry, Pnt

    left = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1)); left.mat("a"); left.faces.name = "outer"
    right = Box(Pnt(0, -1, -1), Pnt(1, 1, 1)); right.mat("b"); right.faces.name = "outer"
    mesh = ng.Mesh(OCCGeometry(Glue([left, right])).GenerateMesh(maxh=0.7))
    with pytest.raises(ValueError, match="exclusive labels|both orientations"):
        material_set_boundary_signs(mesh, ("a",))


@pytest.mark.parametrize("order", [1, 2])
def test_surface_flux_load_matches_the_volume_load(order):
    from radia.kelvin_solver import (
        project_source_total_hodge, solve_magnetostatic_mixed_total_reduced_omega_kelvin)

    mesh = _mesh(0.4)
    source = _divergence_free_source()
    common = dict(mu_r_by_material={"reduced": 1.0, "total": 50.0},
                  reduced_materials=("reduced",), total_materials=("total",),
                  interface_boundary="source_total_interface", order=order,
                  dirichlet_bbbnd="outer", return_system=True,
                  total_source_materials=("total",))
    with ng.TaskManager():
        volume_hodge = project_source_total_hodge(mesh, source, ("total",), order=order)
        surface_hodge = project_source_total_hodge(
            mesh, source, ("total",), order=order, source_load="surface_flux",
            tangential_tolerance=0.2)
        volume = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, volume_hodge["potential"], 1.0, (3.0, 0.0, 0.0),
            total_source_h=volume_hodge["harmonic_field"], **common)
        surface = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, surface_hodge["potential"], 1.0, (3.0, 0.0, 0.0),
            total_source_h=source, total_source_potential=surface_hodge["potential"],
            reduced_source_load="surface_flux", total_source_load="surface_flux", **common)

    gradient_difference = ng.sqrt(ng.Integrate(
        ng.InnerProduct(ng.grad(volume_hodge["potential"]) - ng.grad(surface_hodge["potential"]),
                        ng.grad(volume_hodge["potential"]) - ng.grad(surface_hodge["potential"]))
        * ng.dx(definedon=mesh.Materials("total")), mesh))
    gradient_norm = ng.sqrt(ng.Integrate(
        ng.InnerProduct(ng.grad(volume_hodge["potential"]), ng.grad(volume_hodge["potential"]))
        * ng.dx(definedon=mesh.Materials("total")), mesh))
    assert gradient_difference / gradient_norm < 1e-7
    assert surface_hodge["relative_harmonic_norm"] is None
    assert 0.0 <= surface_hodge["relative_tangential_residual"] < 0.2

    # The test mesh has no point gauge, so compare the reduced-block load and
    # physical fields rather than whole vectors that carry additive constants.
    n_reduced = volume["fes_reduced"].ndof
    loads = [result["system"]["linear_form"].vec.FV().NumPy()[:n_reduced]
             for result in (volume, surface)]
    assert np.linalg.norm(loads[1] - loads[0]) / np.linalg.norm(loads[0]) < 1e-6
    balance = surface["source_loads"]["surface_flux_balance"]
    assert set(balance) == {"reduced", "total:total"}
    assert all(item["relative_flux_imbalance"] < 1e-6 for item in balance.values())
    for point in ((-0.5, 0.15, -0.1), (-0.15, -0.2, 0.25), (0.5, 0.1, 0.2), (0.2, -0.3, -0.1)):
        expected = np.asarray(volume["H_cf"](mesh(*point)))
        actual = np.asarray(surface["H_cf"](mesh(*point)))
        assert np.linalg.norm(actual - expected) <= 1e-7 * max(np.linalg.norm(expected), 1.0)


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("on_boundary", [False, True])
def test_nodal_interpolant_reproduces_its_polynomial_space(order, on_boundary):
    from radia.kelvin_solver import interpolate_source_field

    mesh = _mesh(0.5)  # straight-sided, so physical polynomials are in the space
    region = (mesh.Boundaries("source_total_interface|total_outer") if on_boundary
              else mesh.Materials("total"))
    polynomial = (ng.CF((1 + ng.x + 2 * ng.y - ng.z, 3 * ng.x - ng.y, ng.z - 0.5))
                  if order == 1 else
                  ng.CF((ng.x * ng.x - ng.y * ng.z, ng.x * ng.y + ng.z, 1 - ng.z * ng.z)))
    result = interpolate_source_field(mesh, polynomial, region, order=order)
    measure = (ng.ds(definedon=region) if on_boundary else ng.dx(definedon=region))
    difference = result["field"] - polynomial
    error = ng.Integrate(ng.InnerProduct(difference, difference) * measure, mesh)
    norm = ng.Integrate(ng.InnerProduct(polynomial, polynomial) * measure, mesh)
    assert (error / norm) ** 0.5 < 1e-12
    vertices = {v.nr for element in mesh.Elements(ng.BND if on_boundary else ng.VOL)
                if (element.mat in ("source_total_interface", "total_outer")
                    if on_boundary else element.mat == "total")
                for v in element.vertices}
    assert result["nodes"] >= len(vertices)
    if order == 1:
        assert result["nodes"] == len(vertices)


@pytest.fixture(scope="module")
def nonlinear_case():
    import math
    from netgen.meshing import Element0D
    from radia.kelvin_solver import project_source_interface_potential

    mesh = _mesh(0.5)
    material = mesh.GetMaterials().index("total") + 1
    vertex = next(v for e in mesh.ngmesh.Elements3D() if e.index == material
                  for v in e.vertices if mesh.ngmesh.Points()[v].p[0] > 0.)
    gauge = len(mesh.GetBBBoundaries()) + 1
    mesh.ngmesh.Add(Element0D(vertex, index=gauge))
    mesh.ngmesh.SetCD3Name(gauge, "GND")
    mesh = ng.Mesh(mesh.ngmesh)
    source = _divergence_free_source()
    mu0 = 4.0e-7 * math.pi
    # Saturating, so the constitutive update is genuinely nonlinear here.
    table = ((0.0, 0.0), (0.5, 0.5 * mu0 * 2000.0), (2.0, 1.4e-3), (8.0, 2.0e-3), (40.0, 2.4e-3))
    with ng.TaskManager():
        potential = project_source_interface_potential(
            mesh, source, "source_total_interface", order=2, relative_tolerance=0.05)["potential"]
    return mesh, source, potential, table


_LANES = {
    "picard_centroid_p1": dict(order=1, relaxation=0.3),
    "picard_projected_p2": dict(order=2, material_update_order=1, relaxation=0.3),
    "picard_pointwise_p1": dict(order=1, material_sampling="integration_point",
                                relaxation=1.0, anderson_depth=0),
    "newton_p2": dict(order=2),
}


def _nonlinear_solve(case, lane, load):
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin
    from radia.mixed_omega_newton import solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin

    mesh, source, potential, table = case
    common = dict(bh_table=table, nonlinear_materials=("total",), reduced_materials=("reduced",),
                  total_materials=("total",), interface_boundary="source_total_interface",
                  dirichlet_bbbnd="GND", kelvin_mats=(), reduced_source_load=load)
    with ng.TaskManager():
        if lane.startswith("newton"):
            return solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin(
                mesh, source, potential, 1.0, (3.0, 0.0, 0.0), tolerance=1e-8,
                residual_tolerance=1e-10, **_LANES[lane], **common)
        return solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh, source, potential, 1.0, (3.0, 0.0, 0.0), tolerance=1e-9,
            max_iterations=400, **_LANES[lane], **common)


@pytest.mark.parametrize("lane", sorted(_LANES))
def test_nonlinear_reduced_surface_flux_matches_the_volume_load(nonlinear_case, lane):
    """Air keeps mu0 under a B(H) iron, so the reduced surface identity holds."""
    mesh = nonlinear_case[0]
    volume = _nonlinear_solve(nonlinear_case, lane, "volume")
    surface = _nonlinear_solve(nonlinear_case, lane, "surface_flux")
    assert volume["nonlinear_stats"]["converged"] and surface["nonlinear_stats"]["converged"]
    points = ((-0.5, 0.15, -0.1), (-0.15, -0.2, 0.25), (-0.8, 0.5, 0.5),
              (0.5, 0.1, 0.2), (0.2, -0.3, -0.1), (0.05, 0.0, 0.0))
    expected = np.asarray([volume["B_cf"](mesh(*point)) for point in points])
    actual = np.asarray([surface["B_cf"](mesh(*point)) for point in points])
    # Relative to the field scale: the iron points see a field far below the air's.
    scale = np.linalg.norm(expected, axis=1).max()
    assert np.linalg.norm(actual - expected, axis=1).max() <= 1e-6 * scale


def test_surface_flux_hodge_refuses_a_linked_source():
    from netgen.occ import Box, OCCGeometry, Pnt
    from radia.kelvin_solver import project_source_total_hodge

    # One ring-shaped iron body around the x axis, linked by the line current.
    ring = Box(Pnt(-0.5, -1.0, -1.0), Pnt(0.5, 1.0, 1.0)) - Box(
        Pnt(-0.6, -0.4, -0.4), Pnt(0.6, 0.4, 0.4))
    ring.mat("total")
    ring.faces.name = "ring"
    mesh = ng.Mesh(OCCGeometry(ring).GenerateMesh(maxh=0.35))
    radius2 = ng.y * ng.y + ng.z * ng.z
    linked = ng.CF((0.0, -ng.z / radius2, ng.y / radius2))
    with ng.TaskManager():
        with pytest.raises(RuntimeError, match="tangential residual"):
            project_source_total_hodge(mesh, linked, ("total",), order=2,
                                       source_load="surface_flux", tangential_tolerance=0.05)


def test_surface_flux_needs_its_gate_and_its_permeability():
    from radia.kelvin_solver import (
        MU_0, project_source_total_hodge, solve_magnetostatic_mixed_total_reduced_omega_kelvin)

    mesh = _mesh(0.6)
    source = _divergence_free_source()
    with pytest.raises(ValueError, match="tangential_tolerance"):
        project_source_total_hodge(mesh, source, ("total",), source_load="surface_flux")
    with pytest.raises(ValueError, match="source_load"):
        project_source_total_hodge(mesh, source, ("total",), source_load="nodal")
    common = dict(reduced_materials=("reduced",), total_materials=("total",),
                  interface_boundary="source_total_interface", order=1, dirichlet_bbbnd="outer")
    with ng.TaskManager(), pytest.raises(ValueError, match="assembled permeability"):
        solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, 0.0, 1.0, (3.0, 0.0, 0.0),
            mu_cf=ng.CF(5.0 * MU_0), mu_r_by_material={"reduced": 1.0, "total": 5.0},
            reduced_source_load="surface_flux", **common)
    with pytest.raises(ValueError, match="total_source_potential"):
        solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, 0.0, 1.0, (3.0, 0.0, 0.0), mu_r_by_material={"total": 5.0},
            total_source_load="surface_flux", **common)
