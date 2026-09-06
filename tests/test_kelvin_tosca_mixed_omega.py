"""TOSCA-style total/reduced Omega coupling tests."""

from __future__ import annotations

import numpy as np
import pytest


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
            dirichlet_bbbnd="outer")
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


def test_fixed_magnetization_source_uses_one_global_physical_potential():
    """A PM source preserves potential offsets across both material sides."""
    import ngsolve as ng
    from radia.kelvin_solver import project_source_physical_potential

    mesh = _two_region_mesh(maxh=0.35)
    # H_s = -grad(x) is globally exact across the reduced/total material
    # interface.  A volume projection has one additive gauge, unlike two
    # independently gauged surface traces.
    source_h = ng.CoefficientFunction((-1.0, 0.0, 0.0))
    with ng.TaskManager():
        source = project_source_physical_potential(
            mesh,
            source_h,
            ("reduced", "total"),
            order=2,
            relative_tolerance=1.0e-10,
        )
    assert source["relative_volume_residual"] < 1.0e-10
    potential = source["potential"]
    left = float(potential(mesh(-0.75, 0.0, 0.0)))
    right = float(potential(mesh(0.75, 0.0, 0.0)))
    assert abs((right - left) - 1.5) < 1.0e-8


def test_global_physical_source_potential_rejects_current_linked_field():
    """A non-exact source remains on the explicit cut/cohomology route."""
    import ngsolve as ng
    import pytest
    from radia.kelvin_solver import project_source_physical_potential

    mesh = _two_region_mesh(maxh=0.35)
    with ng.TaskManager(), pytest.raises(RuntimeError, match="globally exact"):
        project_source_physical_potential(
            mesh,
            ng.CoefficientFunction((0.0, -ng.z, ng.y)),
            ("reduced", "total"),
            order=2,
            relative_tolerance=1.0e-3,
        )


def test_total_hodge_projection_retains_a_linked_harmonic_source():
    """A linked curl-free field is split, not rejected or scalarized away."""
    from netgen.occ import Box, Glue, OCCGeometry, Pnt
    import ngsolve as ng
    from radia.kelvin_solver import project_source_total_hodge

    bars = [
        Box(Pnt(-0.5, -1.0, 0.4), Pnt(0.5, 1.0, 1.0)),
        Box(Pnt(-0.5, -1.0, -1.0), Pnt(0.5, 1.0, -0.4)),
        Box(Pnt(-0.5, -1.0, -0.4), Pnt(0.5, -0.4, 0.4)),
        Box(Pnt(-0.5, 0.4, -0.4), Pnt(0.5, 1.0, 0.4)),
    ]
    for bar in bars:
        bar.mat("total")
    mesh = ng.Mesh(OCCGeometry(Glue(bars)).GenerateMesh(maxh=0.35))
    radius2 = ng.y * ng.y + ng.z * ng.z
    linked_h = ng.CoefficientFunction((0.0, -ng.z / radius2, ng.y / radius2))

    with ng.TaskManager():
        source = project_source_total_hodge(
            mesh, linked_h, ("total",), order=2)

    assert source["relative_harmonic_norm"] > 0.5
    reconstructed = -ng.grad(source["potential"]) + source["harmonic_field"]
    error = ng.sqrt(ng.Integrate(
        ng.InnerProduct(reconstructed - linked_h, reconstructed - linked_h),
        mesh,
    ))
    assert float(error) < 1.0e-12


def test_mixed_omega_accepts_the_total_hodge_source_components():
    """The Hodge scalar and harmonic terms both enter the total region."""
    import ngsolve as ng
    from radia.kelvin_solver import (
        project_source_total_hodge,
        solve_magnetostatic_mixed_total_reduced_omega_kelvin,
    )

    mesh = _two_region_mesh(maxh=0.55)
    source_h = ng.CoefficientFunction((-1.0, 0.0, 0.0))
    with ng.TaskManager():
        source = project_source_total_hodge(
            mesh, source_h, ("total",), order=1)
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,
            source_h,
            source["potential"],
            1.0,
            (3.0, 0.0, 0.0),
            mu_r_by_material={"reduced": 1.0, "total": 2.0},
            reduced_materials=("reduced",),
            total_materials=("total",),
            interface_boundary="source_total_interface",
            order=1,
            dirichlet_bbbnd="outer",
            total_source_h=source["harmonic_field"],
            total_source_materials=("total",),
        )

    value = np.asarray(result["H_cf"](mesh(0.5, 0.1, 0.1)), dtype=float)
    assert np.isfinite(value).all()
    assert result["total_source_materials"] == ("total",)


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
            order=2, dirichlet_bbbnd="outer", tolerance=1.0e-9,
            max_iterations=5, relaxation=0.5)

    assert result["nonlinear_stats"]["converged"]
    assert result["nonlinear_stats"]["iterations"] == 2
    assert np.linalg.norm(np.asarray(result["H_cf"](mesh(0.5, 0.1, 0.2)))) > 1.0e-5


def _picard_case(maxh=0.32):
    """Shared source/interface setup of the Picard contract test above."""
    import math
    import ngsolve as ng
    from radia.kelvin_solver import project_source_interface_potential

    mesh = _two_region_mesh(maxh=maxh)
    r_x = ng.x - 2.5
    r2 = r_x * r_x + ng.y * ng.y + ng.z * ng.z
    h_gradient = ng.CoefficientFunction((
        r_x / r2**1.5, ng.y / r2**1.5, ng.z / r2**1.5))
    a = (1.0 - ng.x**2)**2
    b = (1.0 - ng.y**2)**2
    c = (1.0 - ng.z**2)**2
    da_dx = -4.0 * ng.x * (1.0 - ng.x**2)
    db_dy = -4.0 * ng.y * (1.0 - ng.y**2)
    h_source = h_gradient + ng.CoefficientFunction((
        a * db_dy * c, -da_dx * b * c, 0.0))
    mu0 = 4.0e-7 * math.pi
    # A saturating table so the Picard map is genuinely nonlinear at this source.
    bh_table = ((0.0, 0.0), (0.5, 0.5 * mu0 * 2000.0), (2.0, 1.4e-3),
                (8.0, 2.0e-3), (40.0, 2.4e-3))
    with ng.TaskManager():
        trace = project_source_interface_potential(
            mesh, h_source, "source_total_interface", order=2,
            relative_tolerance=0.04)
    return mesh, h_source, trace["potential"], bh_table


def _picard_solve(mesh, h_source, potential, bh_table, **overrides):
    import ngsolve as ng
    from radia.kelvin_solver import (
        solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin,
    )

    settings = dict(
        bh_table=bh_table,
        nonlinear_materials=("total",), reduced_materials=("reduced",),
        total_materials=("total",), interface_boundary="source_total_interface",
        order=2, dirichlet_bbbnd="outer", tolerance=1.0e-6,
        max_iterations=60, relaxation=0.3)
    settings.update(overrides)
    with ng.TaskManager():
        return solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh, h_source, potential, 1.0, (3.0, 0.0, 0.0), **settings)


def test_mixed_omega_picard_records_history_and_per_element_state():
    import ngsolve as ng

    mesh, h_source, potential, bh_table = _picard_case()
    observation = np.array([[0.5, 0.1, 0.2], [0.7, -0.2, 0.1]])
    result = _picard_solve(mesh, h_source, potential, bh_table,
                           observation_points=observation)
    stats = result["nonlinear_stats"]
    assert stats["converged"]
    assert stats["anderson_depth"] == 0
    assert stats["warm_start"] is False
    assert len(stats["history"]) == stats["iterations"]
    assert stats["history"][0]["iteration"] == 1
    assert "relative_B_change" not in stats["history"][0]
    assert all("relative_B_change" in row for row in stats["history"][1:])
    assert stats["history"][-1]["relative_B_change"] == stats["relative_B_change"]
    assert all("observation_relative_change" in row for row in stats["history"][1:])
    assert np.asarray(stats["observation_field_T"]).shape == (2, 3)
    iron = [el.nr for el in mesh.Elements(ng.VOL) if str(el.mat) == "total"]
    assert stats["element_numbers"] == iron
    mu_r = np.asarray(stats["mu_r_elements"])
    assert mu_r.shape == (len(iron),)
    assert np.all(mu_r >= 1.0)
    assert stats["contraction_rate_estimate"] is None or 0.0 < stats["contraction_rate_estimate"] < 1.0


def test_mixed_omega_picard_warm_start_resumes_from_the_converged_state():
    mesh, h_source, potential, bh_table = _picard_case()
    cold = _picard_solve(mesh, h_source, potential, bh_table)
    warm = _picard_solve(mesh, h_source, potential, bh_table,
                         mu_r_initial=np.asarray(cold["nonlinear_stats"]["mu_r_elements"]))
    assert warm["nonlinear_stats"]["warm_start"] is True
    assert warm["nonlinear_stats"]["converged"]
    # A converged warm start needs one solve to reproduce B and one to confirm it.
    assert warm["nonlinear_stats"]["iterations"] <= 3
    assert warm["nonlinear_stats"]["iterations"] < cold["nonlinear_stats"]["iterations"]
    probe = mesh(0.5, 0.1, 0.2)
    np.testing.assert_allclose(
        np.asarray(warm["H_cf"](probe)), np.asarray(cold["H_cf"](probe)), rtol=1.0e-4)
    with pytest.raises(ValueError, match="one value per nonlinear element"):
        _picard_solve(mesh, h_source, potential, bh_table, mu_r_initial=np.ones(3))
    with pytest.raises(ValueError, match=">= 1"):
        bad = np.asarray(cold["nonlinear_stats"]["mu_r_elements"]).copy()
        bad[0] = 0.5
        _picard_solve(mesh, h_source, potential, bh_table, mu_r_initial=bad)


def test_mixed_omega_picard_constrained_anderson_reaches_the_same_solution():
    mesh, h_source, potential, bh_table = _picard_case()
    plain = _picard_solve(mesh, h_source, potential, bh_table)
    mixed = _picard_solve(mesh, h_source, potential, bh_table, anderson_depth=2)
    stats = mixed["nonlinear_stats"]
    assert stats["converged"]
    assert stats["anderson_depth"] == 2
    assert stats["anderson"]["accelerated_steps"] >= 1
    assert stats["iterations"] <= plain["nonlinear_stats"]["iterations"]
    mu_r = np.asarray(stats["mu_r_elements"])
    assert np.all(np.isfinite(mu_r)) and np.all(mu_r >= 1.0)
    probe = mesh(0.5, 0.1, 0.2)
    np.testing.assert_allclose(
        np.asarray(mixed["H_cf"](probe)), np.asarray(plain["H_cf"](probe)), rtol=2.0e-3)


def test_mixed_omega_picard_non_convergence_raises_with_the_state():
    from radia.kelvin_solver import MixedOmegaPicardNotConverged

    mesh, h_source, potential, bh_table = _picard_case()
    with pytest.raises(MixedOmegaPicardNotConverged, match="did not converge") as excinfo:
        _picard_solve(mesh, h_source, potential, bh_table, max_iterations=2, tolerance=1.0e-12)
    state = excinfo.value.state
    assert state["nonlinear_stats"]["converged"] is False
    assert state["nonlinear_stats"]["iterations"] == 2
    assert len(state["nonlinear_stats"]["history"]) == 2
    assert len(state["mu_r_elements"]) == len(state["element_numbers"]) > 0
    # The partial state is a valid warm start.
    resumed = _picard_solve(mesh, h_source, potential, bh_table,
                            mu_r_initial=np.asarray(state["mu_r_elements"]))
    assert resumed["nonlinear_stats"]["converged"]
