"""TOSCA-style total/reduced Omega coupling tests."""

from __future__ import annotations

import numpy as np


def _two_region_mesh(maxh):
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X
    import ngsolve as ng

    reduced = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    reduced.mat("reduced")
    reduced.faces.name = "outer"
    reduced.faces.Max(X).name = "source_total_interface"
    total = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    total.mat("total")
    total.faces.name = "outer"
    total.faces.Min(X).name = "source_total_interface"
    return ng.Mesh(OCCGeometry(Glue([reduced, total])).GenerateMesh(maxh=maxh))


def test_mixed_total_reduced_omega_keeps_source_out_of_high_mu_total_region():
    """The interface jump retains a reduced source without iron cancellation."""
    import ngsolve as ng
    from radia.kelvin_solver import (
        project_source_interface_potential,
        solve_magnetostatic_mixed_total_reduced_omega_kelvin,
    )

    mesh = _two_region_mesh(maxh=0.32)
    r_x = ng.x - 2.5
    r2 = r_x * r_x + ng.y * ng.y + ng.z * ng.z
    h_gradient = ng.CoefficientFunction((
        r_x / r2**1.5, ng.y / r2**1.5, ng.z / r2**1.5))

    # curl((0, 0, f)): a compact source contribution entirely in the
    # reduced enclosure.  It vanishes at the source/total interface, while
    # h_gradient supplies a nonzero, analytically known potential jump.
    a = ng.x**2 * (ng.x + 1.0)**2
    b = (1.0 - ng.y**2)**2
    c = (1.0 - ng.z**2)**2
    da_dx = 2.0 * ng.x * (ng.x + 1.0) * (2.0 * ng.x + 1.0)
    db_dy = -4.0 * ng.y * (1.0 - ng.y**2)
    h_curl = ng.CoefficientFunction((a * db_dy * c, -da_dx * b * c, 0.0))
    h_source = h_gradient + h_curl

    with ng.TaskManager():
        source_trace = project_source_interface_potential(
            mesh, h_source, "source_total_interface", order=2,
            relative_tolerance=0.03)
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, h_source, source_trace["potential"], 1.0, (3.0, 0.0, 0.0),
            mu_r_by_material={"reduced": 1.0, "total": 1000.0},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=2,
            dirichlet_bbnd="outer")
    assert source_trace["relative_tangential_residual"] < 0.03

    h_field = result["H_cf"]
    for point in ((-0.5, 0.15, -0.10), (-0.15, -0.20, 0.25)):
        mip = mesh(*point)
        expected = np.asarray(h_curl(mip), dtype=float)
        actual = np.asarray(h_field(mip), dtype=float)
        assert np.linalg.norm(actual - expected) / np.linalg.norm(expected) < 0.04
    for point in ((0.5, 0.10, 0.20), (0.15, -0.20, 0.25)):
        assert np.linalg.norm(np.asarray(h_field(mesh(*point)), dtype=float)) < 2.0e-3

    phi_reduced = result["phi_reduced"]
    phi_total = result["phi_total"]
    d_interface = ng.ds(definedon=mesh.Boundaries("source_total_interface"))
    jump_error = ng.sqrt(ng.Integrate(
        (phi_total.Trace() - phi_reduced.Trace() - source_trace["potential"])**2 * d_interface,
        mesh))
    assert jump_error < 5.0e-3


def test_mixed_total_reduced_omega_requires_an_exhaustive_material_partition():
    """An omitted Kelvin/iron material must fail before an invalid solve starts."""
    import ngsolve as ng
    import pytest
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin

    mesh = _two_region_mesh(maxh=0.7)
    with pytest.raises(ValueError, match="total_materials"):
        solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, ng.CoefficientFunction((0.0, 0.0, 0.0)), 0.0,
            1.0, (3.0, 0.0, 0.0), mu_r_by_material={"reduced": 1.0},
            reduced_materials=("reduced",), total_materials=(),
            interface_boundary="source_total_interface")


def test_source_trace_rejects_a_non_exact_tangential_field_without_a_cut():
    """A linked/non-exact trace must not be silently turned into an Omega source."""
    import ngsolve as ng
    import pytest
    from radia.kelvin_solver import project_source_interface_potential

    mesh = _two_region_mesh(maxh=0.35)
    # On x=0 this has nonzero surface curl, so it is not -grad_Gamma(Phi).
    non_exact_trace = ng.CoefficientFunction((0.0, -ng.z, ng.y))
    with ng.TaskManager(), pytest.raises(RuntimeError, match="cut/cohomology"):
        project_source_interface_potential(
            mesh, non_exact_trace, "source_total_interface", order=2,
            relative_tolerance=1.0e-3)


def test_mixed_total_reduced_omega_picard_uses_the_same_interface_contract():
    """The nonlinear driver keeps source topology separate from B(H) updates."""
    import math
    import ngsolve as ng
    from radia.kelvin_solver import (
        project_source_interface_potential,
        solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin,
    )

    mesh = _two_region_mesh(maxh=0.32)
    r_x = ng.x - 2.5
    r2 = r_x * r_x + ng.y * ng.y + ng.z * ng.z
    h_gradient = ng.CoefficientFunction((
        r_x / r2**1.5, ng.y / r2**1.5, ng.z / r2**1.5))
    # This curl contribution has a nonzero normal flux at the interface, so
    # the total-potential material sees a genuine response while its
    # tangential trace remains the projected scalar source trace.
    a = (1.0 - ng.x**2)**2
    b = (1.0 - ng.y**2)**2
    c = (1.0 - ng.z**2)**2
    da_dx = -4.0 * ng.x * (1.0 - ng.x**2)
    db_dy = -4.0 * ng.y * (1.0 - ng.y**2)
    h_source = h_gradient + ng.CoefficientFunction((
        a * db_dy * c, -da_dx * b * c, 0.0))
    mu0 = 4.0e-7 * math.pi

    with ng.TaskManager():
        trace = project_source_interface_potential(
            mesh, h_source, "source_total_interface", order=2,
            relative_tolerance=0.04)
        result = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh, h_source, trace["potential"], 1.0, (3.0, 0.0, 0.0),
            bh_table=((0.0, 0.0), (10.0, 10.0 * mu0 * 1000.0)),
            nonlinear_materials=("total",), reduced_materials=("reduced",),
            total_materials=("total",), interface_boundary="source_total_interface",
            order=2, dirichlet_bbnd="outer", tolerance=1.0e-9,
            max_iterations=5, relaxation=0.5)

    assert result["nonlinear_stats"]["converged"]
    assert result["nonlinear_stats"]["iterations"] == 2
    assert np.linalg.norm(np.asarray(result["H_cf"](mesh(0.5, 0.1, 0.2)))) > 1.0e-5
