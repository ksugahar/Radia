"""mixed total/reduced Omega coupling tests."""

from __future__ import annotations

import numpy as np
import pytest


def test_reduced_normal_boundary_is_not_total_normal_boundary():
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import Element0D
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin
    iron=Box(Pnt(-.3,-.3,-.3),Pnt(.3,.3,.3))
    iron.mat('total'); iron.faces.name='interface'
    outer=Box(Pnt(-1,-1,-1),Pnt(1,1,1));outer.faces.name='outer'
    air=outer-iron;air.mat('reduced')
    mesh=ng.Mesh(OCCGeometry(Glue([iron,air])).GenerateMesh(maxh=.8))
    material=mesh.GetMaterials().index('total')+1
    e=next(e for e in mesh.ngmesh.Elements3D() if e.index==material)
    mesh.ngmesh.Add(Element0D(e.vertices[0],index=1));mesh.ngmesh.SetCD3Name(1,'GND')
    mesh=ng.Mesh(mesh.ngmesh)
    source=ng.CoefficientFunction((0.,0.,1.))
    kwargs=dict(mu_r_by_material={'total':1.,'reduced':1.},
                reduced_materials=('reduced',),total_materials=('total',),
                interface_boundary='interface',kelvin_mats=(),order=1)
    with ng.TaskManager():
        correction=solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,source,-ng.z,1.,(0,0,0),reduced_zero_normal_boundary='outer',**kwargs)
        shifted=solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,source,-ng.z,1.,(0,0,0),reduced_zero_normal_boundary='outer',
            total_dirichlet_cf=ng.CoefficientFunction(2.),**kwargs)
        total=solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,source,-ng.z,1.,(0,0,0),**kwargs)
        with pytest.raises(ValueError,match='exterior'):
            solve_magnetostatic_mixed_total_reduced_omega_kelvin(
                mesh,source,-ng.z,1.,(0,0,0),reduced_zero_normal_boundary='interface',**kwargs)
    for point in ((.1,.1,.1),(.6,.1,.1)):
        assert np.linalg.norm(np.asarray(correction['H_cf'](mesh(*point)))-[0,0,1])<1e-10
        assert np.linalg.norm(np.asarray(shifted['H_cf'](mesh(*point)))-[0,0,1])<1e-10
        assert np.linalg.norm(np.asarray(total['H_cf'](mesh(*point))))<1e-10
    point=mesh(.1,.1,.1)
    assert float(shifted['phi_total'](point)-correction['phi_total'](point))==pytest.approx(2.)


def test_lift_only_drive_is_judged_by_its_effective_load():
    """With zero source the free-row load f is zero; only the lift A u_D drives.

    The residual gate must divide by ||P(f - A u_D)||, not by ||P f|| (which
    turned roundoff into a spurious failure).
    """
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from netgen.meshing import Element0D
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin
    iron=Box(Pnt(-.3,-.3,-.3),Pnt(.3,.3,.3))
    iron.mat('total'); iron.faces.name='interface'
    outer=Box(Pnt(-1,-1,-1),Pnt(1,1,1));outer.faces.name='outer'
    air=outer-iron;air.mat('reduced')
    mesh=ng.Mesh(OCCGeometry(Glue([iron,air])).GenerateMesh(maxh=.8))
    material=mesh.GetMaterials().index('total')+1
    e=next(e for e in mesh.ngmesh.Elements3D() if e.index==material)
    mesh.ngmesh.Add(Element0D(e.vertices[0],index=1));mesh.ngmesh.SetCD3Name(1,'GND')
    mesh=ng.Mesh(mesh.ngmesh)
    with ng.TaskManager():
        result=solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,ng.CoefficientFunction((0.,0.,0.)),ng.CoefficientFunction(0.),1.,(0,0,0),
            mu_r_by_material={'total':1.,'reduced':1.},
            reduced_materials=('reduced',),total_materials=('total',),
            interface_boundary='interface',kelvin_mats=(),order=1,
            reduced_zero_normal_boundary='outer',
            total_dirichlet_cf=ng.CoefficientFunction(3.))
    point=mesh(.1,.1,.1)
    assert float(result['phi_total'](point))==pytest.approx(3.)
    assert np.linalg.norm(np.asarray(result['H_cf'](point)))<1e-9


def test_realized_bh_response_binds_pchip_tangent_energy_and_vacuum_tail():
    from radia.scalar_potential_solver import MU_0, sample_bh_constitutive_response

    table = [(0.0, 0.0), (100.0, 0.5), (1000.0, 1.4), (10000.0, 1.7)]
    grid = [0.0, 50.0, 100.0, 500.0, 1000.0, 10000.0, 20000.0]
    response = sample_bh_constitutive_response(table, grid)

    assert response["identity"]["constitutive_interpolation"] == "monotone_pchip"
    assert response["identity"]["constitutive_extrapolation"] == "vacuum_slope"
    assert response["B_T"][-1] == pytest.approx(1.7 + MU_0 * 10000.0)
    assert response["differential_permeability_H_per_m"][-1] == pytest.approx(MU_0)
    for h, b, energy, coenergy in zip(
        response["H_A_per_m"],
        response["B_T"],
        response["energy_density_J_per_m3"],
        response["coenergy_density_J_per_m3"],
    ):
        assert energy + coenergy == pytest.approx(h * b, abs=1.0e-10)


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


def test_mixed_omega_symmetry_potential_boundaries_constrain_both_spaces():
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X, Y
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin

    reduced = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    reduced.mat("reduced")
    reduced.faces.name = "outer"
    reduced.faces.Max(X).name = "source_total_interface"
    reduced.faces.Max(Y).name = "norm_boundary"
    total = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    total.mat("total")
    total.faces.name = "outer"
    total.faces.Min(X).name = "source_total_interface"
    total.faces.Max(Y).name = "norm_boundary"
    mesh = ng.Mesh(OCCGeometry(Glue([reduced, total])).GenerateMesh(maxh=0.7))
    source = ng.CoefficientFunction((0.0, 0.0, 0.0))

    with ng.TaskManager():
        natural = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, 0.0, 1.0, (3.0, 0.0, 0.0),
            mu_r_by_material={"reduced": 1.0, "total": 1.0},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=1,
            dirichlet_bbbnd="outer")
        fixed = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, 0.0, 1.0, (3.0, 0.0, 0.0),
            mu_r_by_material={"reduced": 1.0, "total": 1.0},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=1,
            dirichlet_bbbnd=None,
            reduced_dirichlet_boundary="norm_boundary",
            total_dirichlet_boundary="norm_boundary")

    for part in ("fes_reduced", "fes_total"):
        assert sum(fixed[part].FreeDofs()) < sum(natural[part].FreeDofs())
    assert np.linalg.norm(np.asarray(fixed["H_cf"](mesh(-0.5, 0.0, 0.0)))) < 1e-10
    assert np.linalg.norm(np.asarray(fixed["H_cf"](mesh(0.5, 0.0, 0.0)))) < 1e-10

    with ng.TaskManager():
        loaded = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, ng.CoefficientFunction((0.0, 0.0, 1.0)), 0.0,
            1.0, (3.0, 0.0, 0.0),
            mu_r_by_material={"reduced": 1.0, "total": 1.0},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=1,
            dirichlet_bbbnd=None,
            reduced_dirichlet_boundary="norm_boundary",
            total_dirichlet_boundary="norm_boundary")
    assert loaded["linear_residual"]["free_dofs"]["relative"] < 1e-8
    with ng.TaskManager():
        edge_fixed = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, ng.CoefficientFunction((0.0, 0.0, 1.0)), 0.0,
            1.0, (3.0, 0.0, 0.0),
            mu_r_by_material={"reduced": 1.0, "total": 1.0},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=1,
            dirichlet_bbbnd=None,
            reduced_dirichlet_boundary="norm_boundary",
            total_dirichlet_boundary="norm_boundary",
            interface_multiplier_dirichlet_boundary="norm_boundary")
    assert edge_fixed["linear_residual"]["free_dofs"]["relative"] < 1e-8
    for point in ((-0.4, 0.0, 0.0), (0.4, 0.0, 0.0)):
        assert np.linalg.norm(
            np.asarray(edge_fixed["H_cf"](mesh(*point)))
            - np.asarray(loaded["H_cf"](mesh(*point)))) < 1e-9


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
    phases = []

    with ng.TaskManager():
        source_trace = project_source_interface_potential(
            mesh, h_source, "source_total_interface", order=2,
            relative_tolerance=0.03)
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, h_source, source_trace["potential"], 1.0, (3.0, 0.0, 0.0),
            mu_r_by_material={"reduced": 1.0, "total": 1000.0},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=2,
            dirichlet_bbbnd="outer", phase_callback=phases.append)
    assert source_trace["relative_tangential_residual"] < 0.03
    names = ["matrix_assembly", "source_rhs_assembly", "factorization", "backsolve"]
    assert [event["phase"] for event in phases if event["event"] == "complete"] == names
    assert all(np.isfinite(result["phase_timings_seconds"][name])
               and result["phase_timings_seconds"][name] >= 0 for name in names)

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


@pytest.mark.parametrize("bonus", [0, 6])
def test_total_hodge_projection_retains_a_linked_harmonic_source(bonus):
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
            mesh, linked_h, ("total",), order=2, bonus_intorder=bonus)
        test = source["fes"].TestFunction()
        measure = ng.dx(definedon=mesh.Materials("total"), bonus_intorder=bonus)
        orthogonality = ng.LinearForm(source["fes"])
        orthogonality += (ng.InnerProduct(source["harmonic_field"], ng.grad(test))
                          + 1e-12 * source["potential"] * test) * measure
        orthogonality.Assemble()
        load = ng.LinearForm(source["fes"])
        load += ng.InnerProduct(linked_h, ng.grad(test)) * measure
        load.Assemble()
        assert ng.Norm(orthogonality.vec) / ng.Norm(load.vec) < 1e-10

    assert source["bonus_intorder"] == bonus
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


def test_mixed_omega_accepts_preassembled_reduced_source_load():
    """Replacing only the reduced volume load preserves the full system."""
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin

    mesh = _two_region_mesh(maxh=0.55)
    source_h = ng.CoefficientFunction((ng.x * ng.x, ng.y, ng.z))
    trace = 0.0
    args = dict(mu_r_by_material={"reduced": 1.0, "total": 2.0},
                reduced_materials=("reduced",), total_materials=("total",),
                interface_boundary="source_total_interface", order=1,
                dirichlet_bbbnd="outer", return_system=True)
    with ng.TaskManager():
        original = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source_h, trace, 1.0, (3.0, 0.0, 0.0), **args)
        reduced = original["fes_reduced"]
        test = reduced.TestFunction()
        load = ng.LinearForm(reduced)
        load += (4e-7 * np.pi) * ng.InnerProduct(source_h, ng.grad(test)) * ng.dx(
            definedon=mesh.Materials("reduced"), bonus_intorder=4)
        load.Assemble()
        supplied = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source_h, trace, 1.0, (3.0, 0.0, 0.0),
            source_rhs_reduced=load.vec.FV().NumPy().copy(), **args)

    assert np.allclose(supplied["system"]["linear_form"].vec.FV().NumPy(),
                       original["system"]["linear_form"].vec.FV().NumPy(),
                       rtol=1e-11, atol=1e-11)
    # The test mesh has no BBBND point gauge; compare physical fields rather
    # than additive constants in its scalar potentials.
    for point in ((-0.5, 0.15, 0.1), (-0.2, -0.1, -0.2),
                  (0.2, 0.1, 0.2), (0.6, -0.2, 0.1)):
        assert np.allclose(supplied["H_cf"](mesh(*point)),
                           original["H_cf"](mesh(*point)),
                           rtol=1e-9, atol=1e-9)
    with pytest.raises(ValueError, match="source_rhs_reduced"):
        solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source_h, trace, 1.0, (3.0, 0.0, 0.0),
            source_rhs_reduced=np.zeros(2), **args)


@pytest.mark.parametrize("order", [1, 2])
def test_matching_trace_condensation_matches_multiplier_and_cg(order):
    """The restricted SPD path is the saddle system with its trace eliminated."""
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X
    import ngsolve as ng
    from radia.kelvin_solver import (
        solve_magnetostatic_matching_trace_total_reduced_omega,
        solve_magnetostatic_mixed_total_reduced_omega_kelvin,
    )

    reduced = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    reduced.mat("reduced")
    reduced.faces.name = "natural"
    reduced.faces.Min(X).name = "norm_boundary"
    reduced.faces.Max(X).name = "source_total_interface"
    total = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    total.mat("total")
    total.faces.name = "natural"
    total.faces.Max(X).name = "norm_boundary"
    total.faces.Min(X).name = "source_total_interface"
    mesh = ng.Mesh(OCCGeometry(Glue([reduced, total])).GenerateMesh(maxh=0.65))
    source_h = ng.CoefficientFunction((
        ng.x * ng.x + 0.1 * ng.y,
        0.2 * ng.y,
        -0.15 * ng.z,
    ))
    trace_space = ng.H1(
        mesh, order=order,
        definedon=mesh.Boundaries("source_total_interface"))
    trace = ng.GridFunction(trace_space)
    trace.Set(-ng.z, definedon=mesh.Boundaries("source_total_interface"))
    common = dict(
        mu_r_by_material={"reduced": 1.0, "total": 2.0},
        reduced_materials=("reduced",), total_materials=("total",),
        interface_boundary="source_total_interface", order=order,
    )
    with ng.TaskManager():
        saddle = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source_h, trace, 1.0, (3.0, 0.0, 0.0),
            dirichlet_bbbnd=None,
            reduced_dirichlet_boundary="norm_boundary",
            total_dirichlet_boundary="norm_boundary",
            kelvin_mats=(), **common)
        direct = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, source_h, trace, dirichlet_boundary="norm_boundary",
            solver="direct", **common)
        iterative = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, source_h, trace, dirichlet_boundary="norm_boundary",
            solver="cg", cg_tolerance=1.0e-12, **common)

    assert direct["linear_residual_relative"] < 1.0e-10
    assert iterative["linear_residual_relative"] < 1.0e-9
    for point in ((-0.65, 0.1, 0.15), (-0.2, -0.2, 0.25),
                  (0.2, 0.15, -0.1), (0.65, -0.1, 0.2)):
        saddle_h = np.asarray(saddle["H_cf"](mesh(*point)), dtype=float)
        direct_h = np.asarray(direct["H_cf"](mesh(*point)), dtype=float)
        iterative_h = np.asarray(iterative["H_cf"](mesh(*point)), dtype=float)
        assert np.allclose(direct_h, saddle_h, rtol=1.0e-10, atol=1.0e-11)
        assert np.allclose(iterative_h, direct_h, rtol=1.0e-9, atol=1.0e-10)


def test_matching_trace_condensation_maps_preassembled_reduced_load():
    import ngsolve as ng
    from radia.kelvin_solver import (
        solve_magnetostatic_matching_trace_total_reduced_omega,
    )

    mesh = _two_region_mesh(maxh=0.65)
    trace_space = ng.H1(
        mesh, order=2,
        definedon=mesh.Boundaries("source_total_interface"))
    trace = ng.GridFunction(trace_space)
    trace.vec[:] = 0.0
    args = dict(
        mu_r_by_material={"reduced": 1.0, "total": 1.0},
        reduced_materials=("reduced",), total_materials=("total",),
        interface_boundary="source_total_interface", order=2,
        dirichlet_boundary="outer", solver="direct", return_system=True,
    )
    source_h = ng.CoefficientFunction((ng.x * ng.x, ng.y, ng.z))
    with ng.TaskManager():
        native = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, source_h, trace, **args)
        source_space = native["fes_reduced_source"]
        test = source_space.TestFunction()
        source_load = ng.LinearForm(source_space)
        source_load += ((4.0e-7 * np.pi) * source_h * ng.grad(test)
                        * ng.dx(definedon=mesh.Materials("reduced"),
                                bonus_intorder=4))
        source_load.Assemble()
        supplied = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, source_h, trace,
            source_rhs_reduced=source_load.vec.FV().NumPy().copy(), **args)

        global_space = ng.H1(mesh, order=2)
        global_test = global_space.TestFunction()
        global_load = ng.LinearForm(global_space)
        global_load += ((4.0e-7 * np.pi) * source_h * ng.grad(global_test)
                        * ng.dx(definedon=mesh.Materials("reduced"),
                                bonus_intorder=4))
        global_load.Assemble()
        supplied_global = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, source_h, trace,
            source_rhs_reduced=global_load.vec.FV().NumPy().copy(), **args)

    assert np.allclose(
        native["system"]["linear_form"].vec.FV().NumPy(),
        supplied["system"]["linear_form"].vec.FV().NumPy(),
        rtol=1.0e-12, atol=1.0e-12)
    for point in ((-0.5, 0.1, 0.1), (0.5, -0.1, 0.2)):
        assert np.allclose(native["H_cf"](mesh(*point)),
                           supplied["H_cf"](mesh(*point)),
                           rtol=1.0e-10, atol=1.0e-11)
        assert np.allclose(native["H_cf"](mesh(*point)),
                           supplied_global["H_cf"](mesh(*point)),
                           rtol=1.0e-10, atol=1.0e-11)


def test_matching_trace_condensation_rejects_unsupported_trace_and_order():
    import ngsolve as ng
    from radia.kelvin_solver import (
        solve_magnetostatic_matching_trace_total_reduced_omega,
    )

    mesh = _two_region_mesh(maxh=0.8)
    args = dict(
        mu_r_by_material={"reduced": 1.0, "total": 1.0},
        reduced_materials=("reduced",), total_materials=("total",),
        interface_boundary="source_total_interface",
    )
    with pytest.raises(TypeError, match="projected NGSolve GridFunction"):
        solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, ng.CoefficientFunction((0.0, 0.0, 0.0)), ng.x,
            order=1, **args)
    trace_space = ng.H1(
        mesh, order=2,
        definedon=mesh.Boundaries("source_total_interface"))
    trace = ng.GridFunction(trace_space)
    with pytest.raises(ValueError, match="only order 1 or 2"):
        solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, ng.CoefficientFunction((0.0, 0.0, 0.0)), trace,
            order=3, **args)


def test_mixed_omega_prescribed_reduced_normal_flux_recovers_uniform_source():
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin

    mesh = _two_region_mesh(maxh=0.55)
    source_h = ng.CoefficientFunction((0., 0., 1.))
    flux = (4e-7 * np.pi) * ng.InnerProduct(source_h, ng.specialcf.normal(3))
    with ng.TaskManager():
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source_h, -ng.z, 1.0, (3., 0., 0.),
            mu_r_by_material={"reduced": 1., "total": 1.},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=1,
            dirichlet_bbbnd="outer", kelvin_mats=(),
            reduced_normal_flux=flux, reduced_flux_boundary="outer",
            total_normal_flux=flux, total_flux_boundary="outer")
    for point in ((-.5, .1, .1), (.5, -.1, -.1)):
        assert np.allclose(result["H_cf"](mesh(*point)), (0., 0., 1.),
                           rtol=1e-8, atol=1e-8)


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


def test_ngsolve_bh_coefficient_function_matches_scalar_pchip_and_vacuum_tail():
    import ngsolve as ng
    from radia.scalar_potential_solver import (
        _build_bh_coefficient_function,
        _build_bh_coenergy_coefficient_function,
        _build_bh_coenergy_interpolator,
        _build_bh_interpolator,
    )

    mesh, _, _, bh_table = _picard_case()
    parameter = ng.Parameter(0.0)
    coefficient = _build_bh_coefficient_function(parameter, bh_table)
    coenergy_coefficient = _build_bh_coenergy_coefficient_function(parameter, bh_table)
    scalar = _build_bh_interpolator(bh_table)
    coenergy_scalar = _build_bh_coenergy_interpolator(bh_table)
    point = mesh(0.5, 0.1, 0.2)
    for H_value in (0.0, 0.25, 0.5, 1.1, 5.0, 20.0, 40.0, 100.0):
        parameter.Set(H_value)
        assert float(coefficient(point)) == pytest.approx(
            scalar(H_value), rel=2.0e-12, abs=1.0e-15
        )
        assert float(coenergy_coefficient(point)) == pytest.approx(
            coenergy_scalar(H_value), rel=2.0e-12, abs=1.0e-15
        )
    for H_value in (0.25, 1.1, 20.0, 100.0):
        step = 1.0e-6 * max(H_value, 1.0)
        derivative = (
            coenergy_scalar(H_value + step) - coenergy_scalar(H_value - step)
        ) / (2.0 * step)
        assert derivative == pytest.approx(scalar(H_value), rel=2.0e-8, abs=1.0e-10)


@pytest.mark.parametrize("name,kwargs", [
    ("project_source_interface_potential", {"interface_boundary": "interface"}),
    ("project_source_physical_potential", {"physical_materials": ("total",)}),
    ("project_source_total_hodge", {"total_source_materials": ("total",)}),
])
def test_source_projections_reject_retired_direct_backend_before_mesh_access(name, kwargs):
    import radia.kelvin_solver as kelvin
    with pytest.raises(ValueError, match="sparsecholesky"):
        getattr(kelvin, name)(None, None, inverse="pardiso", **kwargs)


def test_source_projection_diagnostics_survive_many_threads():
    """Exercise the fixed heap share without changing the suite's thread count."""
    import subprocess
    import sys
    from pathlib import Path

    # NGSolve has no public getter for the caller's configured thread count.
    # Keep the stress setting in a child instead of guessing how to restore it.
    script = r"""
import runpy, sys
import ngsolve as ng
import numpy as np
from radia.kelvin_solver import project_source_total_hodge
ng.SetNumThreads(1)
# Bound assembly heap reservation in this child (76 * 10 MB). The regression
# still exercises 76 threads and the same high-order integral; Integrate's
# fixed per-thread diagnostic heap is unaffected by SetHeapSize.
ng.SetHeapSize(10_000_000)
case = runpy.run_path(sys.argv[1])["_picard_case"]
mesh, source, _, _ = case(maxh=0.35)
values = []
for threads in (1, 76):
    ng.SetNumThreads(threads)
    with ng.TaskManager():
        result = project_source_total_hodge(mesh, source, ("total",), order=2, bonus_intorder=12)
        assert result["inverse"] == "sparsecholesky"
        assert result["linear_true_relative_residual"] < 1e-8
    values.append(result["relative_harmonic_norm"])
np.testing.assert_allclose(values[0], values[1], rtol=1e-10, atol=0.)
"""
    result = subprocess.run(
        [sys.executable, "-c", script, str(Path(__file__).resolve())],
        capture_output=True, text=True, timeout=120, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


def test_trace_without_tolerance_is_not_accepted():
    import ngsolve as ng
    from radia.kelvin_solver import project_source_interface_potential
    mesh = _two_region_mesh(.8)
    with ng.TaskManager():
        result = project_source_interface_potential(
            mesh, ng.CF((0, 0, 1)), 'source_total_interface')
    assert result['gate_enabled'] is False
    assert result['acceptance'] == 'not_evaluated'
    for tolerance in (float('nan'), float('inf'), -1):
        with pytest.raises(ValueError, match='positive and finite'):
            project_source_interface_potential(
                mesh, ng.CF((0, 0, 1)), 'source_total_interface', relative_tolerance=tolerance)


@pytest.mark.parametrize("order", [1, 2, 3])
def test_sparse_direct_saddle_solve_satisfies_original_system(order):
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin
    mesh, source, potential, _ = _picard_case(maxh=0.7)
    with ng.TaskManager():
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, potential, 1., (3., 0., 0.),
            mu_r_by_material={"reduced": 1., "total": 1000.},
            reduced_materials=("reduced",), total_materials=("total",),
            interface_boundary="source_total_interface", order=order,
            dirichlet_bbbnd="outer", inverse="sparsecholesky", return_system=True)
    assert np.isfinite(result["solution"].vec.FV().NumPy()).all()
    assert result["linear_residual"]["free_dofs"]["relative"] < 1e-8
    # A manufactured full-system RHS exercises nonzero multiplier loads and
    # multiplier recovery independently of the physical source above.
    from radia.kelvin_solver import _matching_trace_direct_inverse
    matrix = result["system"]["bilinear_form"].mat
    fes = result["fes"]
    free = np.asarray(list(fes.FreeDofs()), dtype=bool)
    exact = matrix.CreateColVector()
    exact.FV().NumPy()[:] = np.random.default_rng(813).normal(size=fes.ndof)
    exact.FV().NumPy()[~free] = 0.
    rhs = matrix.CreateColVector()
    rhs.data = matrix * exact
    with ng.TaskManager():
        computed = _matching_trace_direct_inverse(matrix, fes, order=order) * rhs
    residual = matrix.CreateColVector()
    residual.data = matrix * computed - rhs
    assert np.linalg.norm(residual.FV().NumPy()[free]) / np.linalg.norm(
        rhs.FV().NumPy()[free]) < 1e-10
    # Reuse the same space while both constitutive and interface coefficients
    # change. Stale numeric factors would fail the original-system residual.
    from scipy.sparse import coo_matrix, diags
    rows, cols, values = matrix.COO()
    assembled = coo_matrix((np.asarray(values), (np.asarray(rows), np.asarray(cols))),
                           shape=(fes.ndof, fes.ndof)).tocsr()
    if isinstance(matrix, ng.la.SparseMatrixSymmetricdouble):
        assembled = assembled + assembled.T - diags(assembled.diagonal())
    entries = assembled.tocoo()
    primal_end = fes.Range(2).start
    is_primal = (entries.row < primal_end) & (entries.col < primal_end)
    cache = {}
    for primal_scale, trace_scale in ((1., 1.), (2., 1.), (.5, 3.)):
        scaled = entries.data * np.where(is_primal, primal_scale, trace_scale)
        changed = ng.la.SparseMatrixdouble.CreateFromCOO(
            entries.row, entries.col, scaled, *assembled.shape)
        rhs.data = changed * exact
        with ng.TaskManager():
            computed = _matching_trace_direct_inverse(
                changed, fes, order=order, cache=cache) * rhs
        residual.data = changed * computed - rhs
        assert np.linalg.norm(residual.FV().NumPy()[free]) / np.linalg.norm(
            rhs.FV().NumPy()[free]) < 1e-10
