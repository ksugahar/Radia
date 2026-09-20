"""TOSCA-style total/reduced Omega coupling tests."""

from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.parametrize('harmonic', [False, True, 'varying'])
def test_picard_fixed_rhs_cache_matches_full_reassembly(harmonic):
    import ngsolve as ng
    mesh, source, potential, table = _picard_case(maxh=0.7)
    options = dict(tolerance=1e-8, max_iterations=100)
    if harmonic:
        options.update(total_source_h=ng.CoefficientFunction((0.01, 0.02, 0.03)),
                       total_source_materials=('total',))
    if harmonic == 'varying':
        options['total_source_h'] = ng.CoefficientFunction((
            .03*ng.exp(3*ng.x)*ng.cos(3*ng.y),
            -.03*ng.exp(3*ng.x)*ng.sin(3*ng.y), 0.))
    fresh = _picard_solve(mesh, source, potential, table, cache_fixed_rhs=False, **options)
    cached = _picard_solve(mesh, source, potential, table, cache_fixed_rhs=True, **options)
    for point in ((-.5,.13,.17),(.5,.13,.17)):
        a=np.asarray(fresh['H_cf'](mesh(*point)))
        b=np.asarray(cached['H_cf'](mesh(*point)))
        assert np.linalg.norm(a-b) < 1e-10 * max(np.linalg.norm(a),1.)
    assert fresh['nonlinear_stats']['iterations'] == cached['nonlinear_stats']['iterations']
    assert fresh['nonlinear_stats']['rhs_reuse'] == 'none'
    assert cached['nonlinear_stats']['rhs_reuse'] == (
        'fixed_vector_and_material_operator' if harmonic else 'fixed_vector')


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
        order=1, dirichlet_bbbnd="outer", tolerance=1.0e-6,
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


def test_mixed_omega_picard_rejects_unmatched_high_order_nonlinear_material_update():
    mesh, h_source, potential, bh_table = _picard_case()
    with pytest.raises(ValueError, match="order-0 centroid value per element"):
        _picard_solve(mesh, h_source, potential, bh_table, order=2)


def test_mixed_omega_picard_projects_positive_order_matched_material_state():
    mesh, h_source, potential, bh_table = _picard_case()
    observation = np.array([[0.5, 0.1, 0.2]])
    result = _picard_solve(
        mesh,
        h_source,
        potential,
        bh_table,
        order=2,
        material_update_order=1,
        anderson_depth=0,
        observation_points=observation,
    )
    stats = result["nonlinear_stats"]
    assert stats["converged"]
    assert stats["response_order"] == 2
    assert stats["material_update_order"] == 1
    assert stats["material_sampling"] == "L2_log_secant_projection"
    assert stats["physical_permeability_bounds"][0] == 1.0
    assert stats["final_material_state_resolved"]
    assert len(stats["material_log_state_dofs"]) == len(
        stats["material_state_dof_numbers"]
    )
    assert len(stats["material_log_state_dofs"]) > 0
    np.testing.assert_allclose(
        stats["observation_field_T"],
        [np.asarray(result["B_cf"](mesh(*observation[0])))],
        rtol=1.0e-12,
        atol=1.0e-12,
    )

    resumed = _picard_solve(
        mesh,
        h_source,
        potential,
        bh_table,
        order=2,
        material_update_order=1,
        anderson_depth=0,
        material_log_state_initial=stats["material_restart_state"],
    )
    resumed_stats = resumed["nonlinear_stats"]
    assert resumed_stats["converged"]
    assert resumed_stats["warm_start"]
    assert resumed_stats["final_material_state_resolved"]


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


def test_mixed_omega_projected_material_state_validates_resume_shape():
    mesh, h_source, potential, bh_table = _picard_case()
    solved = _picard_solve(
        mesh,
        h_source,
        potential,
        bh_table,
        order=2,
        material_update_order=1,
        anderson_depth=0,
    )
    restart = dict(solved["nonlinear_stats"]["material_restart_state"])
    restart["values"] = [0.0]
    with pytest.raises(ValueError, match="one value per active"):
        _picard_solve(
            mesh,
            h_source,
            potential,
            bh_table,
            order=2,
            material_update_order=1,
            anderson_depth=0,
            material_log_state_initial=restart,
        )


def test_mixed_omega_projected_material_state_rejects_identity_mismatch():
    mesh, h_source, potential, bh_table = _picard_case()
    solved = _picard_solve(
        mesh,
        h_source,
        potential,
        bh_table,
        order=2,
        material_update_order=1,
        anderson_depth=0,
    )
    restart = solved["nonlinear_stats"]["material_restart_state"]
    changed_table = list(bh_table)
    changed_table[-1] = (changed_table[-1][0], changed_table[-1][1] * 1.01)
    with pytest.raises(ValueError, match="does not match the current mesh"):
        _picard_solve(
            mesh,
            h_source,
            potential,
            changed_table,
            order=2,
            material_update_order=1,
            anderson_depth=0,
            material_log_state_initial=restart,
        )


def test_mixed_omega_projected_energy_uses_current_solution_and_material_state():
    mesh, h_source, potential, bh_table = _picard_case()
    solved = _picard_solve(
        mesh,
        h_source,
        potential,
        bh_table,
        order=2,
        material_update_order=1,
        anderson_depth=0,
        refinement_parent_identity={
            "geometry": "two-box-interface-v1",
            "material": "saturating-table-v1",
            "excitation": "fixed-source-field-v1",
        },
    )
    energy = solved["energy_observables"]
    stats = solved["nonlinear_stats"]
    assert energy["energy_J"] > 0.0
    assert energy["coenergy_J"] > 0.0
    assert energy["identity"]["material_state_identity_sha256"] == stats[
        "material_state_identity_sha256"
    ]
    assert energy["identity"]["bh_table_sha256"] == stats[
        "material_state_identity"
    ]["bh_table_sha256"]
    assert len(energy["identity"]["solution_sha256"]) == 64
    assert len(energy["identity"]["refinement_parent_identity_sha256"]) == 64
    assert energy["identity"]["constitutive_interpolation"] == "monotone_pchip"
    assert energy["identity"]["constitutive_extrapolation"] == "vacuum_slope"
    assert energy["identity"]["magnetic_anisotropy"] == "isotropic"
    assert energy["h_dot_b_integral_J"] > 0.0
    assert energy["legendre_residual_relative"] < 1.0e-12
    np.testing.assert_allclose(
        energy["energy_J"] + energy["coenergy_J"],
        energy["h_dot_b_integral_J"],
        rtol=1.0e-12,
        atol=1.0e-14,
    )


def test_mixed_omega_projected_rejects_invalid_refinement_parent_identity():
    mesh, h_source, potential, bh_table = _picard_case()
    with pytest.raises(ValueError, match="must be a mapping"):
        _picard_solve(
            mesh,
            h_source,
            potential,
            bh_table,
            order=2,
            material_update_order=1,
            anderson_depth=0,
            refinement_parent_identity="stale-parent",
        )


def test_mixed_omega_envelope_includes_interpolated_material_targets(monkeypatch):
    import math
    from radia.kelvin_solver import MixedOmegaPicardNotConverged
    from radia.picard_acceleration import ConstrainedAndersonAccelerator

    mesh, h_source, potential, _ = _picard_case()
    mu0 = 4e-7 * math.pi
    table = [(0, 0), (1e-4, mu0 * 0.2), (2e-4, mu0 * 0.4),
             (3e-4, mu0 * 0.402), (1.0, mu0 * 200)]
    maxima = []
    original = ConstrainedAndersonAccelerator.step

    def checked_step(self, current, target):
        maxima.append(float(np.max(target)))
        assert self.upper >= np.max(target)
        result = original(self, current, target)
        np.testing.assert_array_equal(result, 0.3 * target + 0.7 * current)
        return result

    monkeypatch.setattr(ConstrainedAndersonAccelerator, "step", checked_step)
    try:
        _picard_solve(mesh, h_source, potential, table, max_iterations=4,
                     tolerance=1e-12, anderson_depth=0)
    except MixedOmegaPicardNotConverged:
        pass
    assert maxima and max(maxima) > 2000.0


@pytest.mark.parametrize('table', [((1., 0.), (2., .01)), ((0., .001), (1., .01))])
def test_picard_requires_soft_magnetic_origin(table):
    mesh, source, potential, _ = _picard_case(maxh=.8)
    with pytest.raises(ValueError, match='start at'):
        _picard_solve(mesh, source, potential, table)


def test_picard_returns_assembled_material_state(monkeypatch):
    import radia.kelvin_solver as solver
    mesh, source, potential, table = _picard_case(maxh=.8)
    original = solver.solve_magnetostatic_mixed_total_reduced_omega_kelvin
    samples = []
    point = mesh(.5, .13, .17)

    def record(*args, **kwargs):
        result = original(*args, **kwargs)
        samples.append(float(result['mu_cf'](point)))
        return result

    monkeypatch.setattr(solver, 'solve_magnetostatic_mixed_total_reduced_omega_kelvin', record)
    result = _picard_solve(mesh, source, potential, table, tolerance=.01, max_iterations=100)
    assert float(result['mu_cf'](point)) == samples[-1]
    assert result['nonlinear_stats']['relative_constitutive_change'] <= .01


def test_picard_zero_field_material_is_independent_of_initial_guess():
    import ngsolve as ng
    from scipy.interpolate import PchipInterpolator
    from radia.kelvin_material import MU_0
    mesh, _, _, table = _picard_case(maxh=.8)
    expected = float(PchipInterpolator(*np.asarray(table).T).derivative()(0))
    for initial in (10., 1000.):
        result = _picard_solve(mesh, ng.CF((0, 0, 0)), ng.CF(0), table,
                               mu_r_initial=initial)
        assert result['mu_cf'](mesh(.5, .13, .17)) == pytest.approx(expected)
        assert result['nonlinear_stats']['mu_r_elements'][0] == pytest.approx(expected / MU_0)


@pytest.mark.parametrize('order', [1, 2])
def test_picard_hybrid_linear_material_is_preserved(order):
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    from radia.kelvin_material import MU_0
    air = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    iron = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    pole = Box(Pnt(1, -1, -1), Pnt(2, 1, 1))
    air.mat('reduced'); iron.mat('total'); pole.mat('pole')
    air.faces.name = iron.faces.name = pole.faces.name = 'source_total_interface'
    mesh = ng.Mesh(OCCGeometry(Glue([air, iron, pole])).GenerateMesh(maxh=1))
    table = [(0, 0), (1, 100*MU_0), (2, 150*MU_0)]
    options = dict(total_materials=('total', 'pole'), order=order,
                   material_update_order=order-1, tolerance=1e-6, max_iterations=100)
    with pytest.raises(ValueError, match='linear total materials'):
        _picard_solve(mesh, ng.CF((0, 0, 0)), ng.CF(0), table, **options)
    result = _picard_solve(mesh, ng.CF((0, 0, 0)), ng.CF(0), table,
                           mu_r_by_material={'pole': 5000}, **options)
    assert result['mu_cf'](mesh(1.5, .13, .17)) / MU_0 == pytest.approx(5000)
    with pytest.raises(ValueError, match='override nonlinear'):
        _picard_solve(mesh, ng.CF((0, 0, 0)), ng.CF(0), table,
                       mu_r_by_material={'pole': 5000, 'total': 1}, **options)


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
