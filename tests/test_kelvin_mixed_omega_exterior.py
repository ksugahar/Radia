"""The mixed total/reduced Omega Kelvin exterior must carry the field.

A space restricted to the total materials owns no degree of freedom on the
physical sphere, because the material touching that sphere is the reduced
source enclosure.  The periodic identification then pairs nothing, the whole
exterior block is eliminated, and the open boundary silently degenerates into a
Dirichlet truncation at the sphere.  These tests lock the coupling that the
:func:`radia.kelvin_solver.solve_magnetostatic_mixed_total_reduced_omega_kelvin`
formulation actually needs.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

RADIUS = 1.0
OFFSET = (3.0, 0.0, 0.0)
IRON_RADIUS = 0.35
COIL_RADIUS = 0.60
COIL_CURRENT = 1000.0


def _kelvin_mesh(maxh=0.26):
    """Iron ball, reduced air shell, and the offset Kelvin exterior sphere."""
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
    from netgen.meshing import IdentificationType
    import ngsolve as ng

    iron = Sphere(Pnt(0, 0, 0), IRON_RADIUS)
    iron.mat("iron")
    iron.maxh = 0.10
    for face in iron.faces:
        face.name = "iron_air_interface"
    shell = Sphere(Pnt(0, 0, 0), RADIUS)
    shell.mat("air")
    for face in shell.faces:
        face.name = "kelvin_int"
    air = shell - iron
    air.mat("air")

    exterior = Sphere(Pnt(*OFFSET), RADIUS)
    exterior.mat("kelvin")
    for face in exterior.faces:
        face.name = "kelvin_ext"

    # The gauge must sit on a real mesh vertex.  A free-floating point inside
    # a solid is registered as a BBBND element but stays unconnected, so its H1
    # degree of freedom is unused and constraining it changes nothing.  The
    # sphere pole is a genuine topological vertex.
    for vertex in exterior.vertices:
        if vertex.p[2] > 0.0:
            vertex.name = "GND"

    geometry = Glue([air, iron, exterior])
    inner_face = outer_face = None
    for solid in geometry.solids:
        for face in solid.faces:
            if face.name == "kelvin_int":
                inner_face = face
            elif face.name == "kelvin_ext":
                outer_face = face
    inner_face.Identify(outer_face, "periodic", IdentificationType.PERIODIC)
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(geometry).GenerateMesh(maxh=maxh, grading=0.6))


def _coil():
    """A closed circular loop inside the reduced air shell.

    A current source is what the reduced/total split is built for: it is
    divergence free everywhere and its circulation lives entirely inside the
    reduced enclosure, so both interface traces admit a single-valued scalar
    potential and the exterior stays current free.
    """
    import radia as rad
    from radia.coil_builder import CoilBuilder

    builder = CoilBuilder(COIL_CURRENT).set_start((COIL_RADIUS, 0.0, 0.0))
    builder = builder.set_cross_section(0.05, 0.05)
    for _ in range(4):
        builder = builder.add_arc(COIL_RADIUS, 90.0)
    if not builder.is_closed or builder.gap > 1.0e-12:
        raise RuntimeError("test coil path is open by %.3e m" % builder.gap)
    return rad.ObjCnt(builder.to_radia(arc_max_segment_length=0.05))


def _solve(mesh, coil, mu_r_iron=1.0, kelvin_scale=1.0,
           kelvin_interface="kelvin_int", exact_exterior_source=False,
           return_system=False):
    import ngsolve as ng
    import radia as rad
    from radia.kelvin_material import make_kelvin_mu_cf, MU_0
    from radia.kelvin_solver import (
        project_source_interface_potential,
        solve_magnetostatic_mixed_total_reduced_omega_kelvin,
    )

    h_s = rad.RadiaField(coil, "h")
    with ng.TaskManager():
        source = project_source_interface_potential(
            mesh, h_s, "iron_air_interface", order=2, relative_tolerance=0.08)
        kelvin_source = project_source_interface_potential(
            mesh, h_s, "kelvin_int", order=2, relative_tolerance=0.08)
        kelvin_mu = make_kelvin_mu_cf(mesh, RADIUS, OFFSET,
                                      kelvin_mats=("kelvin",),
                                      mu_r_by_material={})
        mu_cf = mesh.MaterialCF({
            "iron": mu_r_iron * MU_0,
            "air": MU_0,
            "kelvin": kelvin_scale * kelvin_mu,
        }, default=MU_0)
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, h_s, source["potential"], RADIUS, OFFSET,
            mu_cf=mu_cf, reduced_materials=("air",),
            total_materials=("iron", "kelvin"),
            interface_boundary="iron_air_interface", order=2,
            kelvin_interface_boundary=kelvin_interface,
            kelvin_source_potential=(
                None if exact_exterior_source else kelvin_source["potential"]),
            kelvin_source_h=(
                rad.KelvinRadiaFieldStrength(coil, OFFSET, RADIUS, (0.0, 0.0, 0.0))
                if exact_exterior_source else None),
            total_source_h=None, total_source_materials=(),
            return_system=return_system)
    return result


def test_kelvin_exterior_carries_the_field():
    """A dead exterior integrates to exactly zero; a coupled one does not."""
    import ngsolve as ng

    mesh = _kelvin_mesh()
    result = _solve(mesh, _coil())
    with ng.TaskManager():
        exterior = float(ng.Integrate(
            result["mu_cf"] * ng.InnerProduct(result["H_cf"], result["H_cf"]),
            mesh, definedon=mesh.Materials("kelvin"), order=6))
        physical = float(ng.Integrate(
            result["mu_cf"] * ng.InnerProduct(result["H_cf"], result["H_cf"]),
            mesh, definedon=mesh.Materials("iron|air"), order=6))
    assert physical > 0.0
    assert exterior / physical > 1.0e-4, (
        "the Kelvin exterior holds no energy, so the open boundary has "
        f"degenerated into a truncation: exterior={exterior:.6e}, "
        f"physical={physical:.6e}")


def test_kelvin_exterior_influences_the_physical_solution():
    """Perturbing the exterior medium must move the physical field."""
    import ngsolve as ng

    mesh = _kelvin_mesh()
    points = ((0.5, 0.1, 0.2), (-0.45, 0.2, -0.15), (0.2, -0.55, 0.1))
    coil = _coil()
    base = _solve(mesh, coil)
    perturbed = _solve(mesh, coil, kelvin_scale=4.0)
    reference = np.asarray([
        np.asarray(base["B_cf"](mesh(*point)), dtype=float) for point in points])
    changed = np.asarray([
        np.asarray(perturbed["B_cf"](mesh(*point)), dtype=float) for point in points])
    relative = float(np.linalg.norm(reference - changed) / np.linalg.norm(reference))
    assert relative > 1.0e-6, (
        "scaling the Kelvin permeability left the physical field unchanged, so "
        f"the exterior transmits nothing: relative change {relative:.3e}")


def test_kelvin_exterior_reproduces_the_pulled_back_source_field():
    """With no iron contrast the exterior must hold the source field itself.

    This is the sign-sensitive check.  The magnetic scalar potential is a
    twisted 0-form, so an exterior coupled with the wrong orientation returns
    the negated field rather than a slightly inaccurate one.
    """
    import radia as rad
    from radia.kelvin_source import evaluate_kelvin_exterior

    mesh = _kelvin_mesh(maxh=0.20)
    coil = _coil()
    result = _solve(mesh, coil, mu_r_iron=1.0)
    directions = np.asarray([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                             (0.6, -0.6, 0.52)], dtype=float)
    directions = directions / np.linalg.norm(directions, axis=1)[:, None]
    physical_points = (RADIUS / 0.62) * directions
    recovered = evaluate_kelvin_exterior(
        result["H_cf"], mesh, physical_points, kelvin_center=OFFSET,
        physical_center=(0.0, 0.0, 0.0), radius=RADIUS, form="field_strength")
    expected = np.asarray(
        rad.Fld(coil, "h", physical_points.tolist()), dtype=float).reshape(-1, 3)
    errors = list(np.linalg.norm(recovered - expected, axis=1)
                  / np.linalg.norm(expected, axis=1))
    # A dead exterior returns exactly the source magnitude as the error (1.0);
    # a sign-flipped one returns about 2.0.  The band below is the coarse-mesh
    # discretisation error of a correctly coupled exterior.
    assert max(errors) < 0.25, (
        "the exterior field does not reproduce the pulled-back source: "
        f"relative errors {['%.3f' % value for value in errors]}")


def test_kelvin_source_lift_is_shared_across_the_identification():
    """The lift must reach the Kelvin sphere through the periodic pairing."""
    import ngsolve as ng

    mesh = _kelvin_mesh()
    result = _solve(mesh, _coil())
    lift = result["kelvin_source_lift"]
    assert lift is not None
    with ng.TaskManager():
        inner = float(ng.Integrate(lift * lift, mesh,
                                   definedon=mesh.Boundaries("kelvin_int"),
                                   order=4))
        outer = float(ng.Integrate(lift * lift, mesh,
                                   definedon=mesh.Boundaries("kelvin_ext"),
                                   order=4))
    assert inner > 0.0
    assert abs(outer / inner - 1.0) < 1.0e-10


def test_direct_solve_reports_its_own_residual():
    """A factorisation has no iteration history, but it still has r = b - A x.

    Reporting nothing because there is no iteration count is the same gap as
    treating a Krylov solve that ran out of iterations as converged.  The free
    degrees of freedom carry the system that was actually solved, and the
    multiplier block carries the interface jump condition in weak form, so the
    two are reported apart: a constraint defect must not be able to hide
    inside a healthy PDE residual.
    """
    mesh = _kelvin_mesh()
    result = _solve(mesh, _coil())
    residual = result["linear_residual"]

    assert residual["free_dofs"]["relative"] < 1.0e-8
    assert residual["free_dofs"]["rhs_l2"] > 0.0
    assert set(residual["blocks"]) == {
        "phi_reduced", "phi_total", "interface_constraint"}
    constraint = residual["blocks"]["interface_constraint"]
    assert constraint is not None
    assert constraint["relative"] < 1.0e-8


def test_assembled_system_is_not_retained_by_default():
    """The assembled matrix must be freed after a production solve.

    The system is exposed for diagnostics that embed another order's solution,
    but the wrappers return the solver's result dict to the caller unchanged;
    carrying the matrix in it by default would pin it in memory for every
    solve.  The scalar diagnostics stay: they cost nothing to keep.
    """
    mesh = _kelvin_mesh()
    result = _solve(mesh, _coil())
    assert result["system"] is None
    assert set(result["assembled_energy"]) >= {"energy", "half_xAx", "b_dot_x"}
    assert result["linear_residual"]["free_dofs"]["relative"] < 1.0e-8


def test_assembled_system_is_available_on_explicit_request():
    result = _solve(_kelvin_mesh(), _coil(), return_system=True)
    system = result["system"]
    assert system["bilinear_form"].mat.height == result["fes"].ndof
    assert len(system["linear_form"].vec) == result["fes"].ndof


def test_kelvin_material_without_its_interface_is_rejected():
    """Silently uncoupling the exterior is a fail-loud configuration error."""
    mesh = _kelvin_mesh()
    with pytest.raises(ValueError, match="kelvin_interface_boundary"):
        _solve(mesh, _coil(), kelvin_interface=None)


def test_exact_pulled_back_exterior_source_matches_the_lift():
    """The two exterior source representations must agree, exactness aside.

    The lift only reproduces the projected interface trace; the exact field is
    the analytic pullback throughout the ball and needs no projection at all.
    They therefore differ by discretisation, not by formulation, and both must
    reproduce the source field in the exterior.
    """
    import radia as rad
    from radia.kelvin_source import kelvin_pullback_vector

    mesh = _kelvin_mesh(maxh=0.20)
    coil = _coil()
    lifted = _solve(mesh, coil, mu_r_iron=1.0)
    exact = _solve(mesh, coil, mu_r_iron=1.0, exact_exterior_source=True)
    assert exact["kelvin_source_lift"] is None
    assert lifted["kelvin_source_lift"] is not None

    centre = np.asarray(OFFSET, dtype=float)
    for direction in ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.6, -0.6, 0.52)):
        unit = np.asarray(direction, dtype=float)
        unit = unit / np.linalg.norm(unit)
        computational = centre + 0.62 * RADIUS * unit
        physical = (RADIUS * RADIUS / (0.62 * RADIUS)) * unit
        expected = np.asarray(
            rad.Fld(coil, "h", [physical.tolist()]), dtype=float).reshape(3)
        errors = []
        for result in (lifted, exact):
            h_comp = np.asarray(
                result["H_cf"](mesh(*map(float, computational))), dtype=float)
            h_physical = -kelvin_pullback_vector(
                h_comp, physical, np.zeros(3), RADIUS)
            errors.append(float(np.linalg.norm(h_physical - expected)
                                / np.linalg.norm(expected)))
        assert max(errors) < 0.25, (
            "exterior field wrong for one representation: lift %.3f, exact %.3f"
            % (errors[0], errors[1]))


def test_solution_pullback_inverts_the_native_kelvin_source_exactly():
    """The pullback API must undo the native forward transform to round-off.

    The native ``KelvinRadia*`` coefficients ARE the forward pullback of a
    known Radia field, so composing them with the solution pullback has to
    return the physical field itself.  That pins the metric factor, the
    Householder reflection, the two-sphere offset, and the orientation sign of
    each form at once, without needing a solve.
    """
    import radia as rad
    from radia.kelvin_source import (
        evaluate_kelvin_exterior,
        kelvin_computational_to_physical,
        kelvin_physical_to_computational,
    )

    mesh = _kelvin_mesh(maxh=0.34)
    coil = _coil()
    physical_center = (0.0, 0.0, 0.0)
    directions = np.asarray([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                             (0.5, -0.7, 0.51)], dtype=float)
    directions = directions / np.linalg.norm(directions, axis=1)[:, None]
    points = (RADIUS / 0.55) * directions

    computational = kelvin_physical_to_computational(
        points, OFFSET, physical_center, RADIUS)
    back = kelvin_computational_to_physical(
        computational, OFFSET, physical_center, RADIUS)
    assert np.allclose(back, points, rtol=1e-12, atol=1e-14)

    for form, factory, component in (
            ("field_strength", rad.KelvinRadiaFieldStrength, "h"),
            ("flux_density", rad.KelvinRadiaFluxDensity, "b")):
        forward = factory(coil, OFFSET, RADIUS, physical_center)
        recovered = evaluate_kelvin_exterior(
            forward, mesh, points, kelvin_center=OFFSET,
            physical_center=physical_center, radius=RADIUS, form=form)
        expected = np.asarray(
            rad.Fld(coil, component, points.tolist()), dtype=float).reshape(-1, 3)
        relative = (np.linalg.norm(recovered - expected, axis=1)
                    / np.linalg.norm(expected, axis=1))
        assert float(np.max(relative)) < 1.0e-9, (
            "%s pullback does not invert the native transform: %s"
            % (form, relative))


def test_solution_pullback_rejects_interior_points_and_unknown_forms():
    """Fail loud rather than transform something the formula does not cover."""
    from radia.kelvin_source import evaluate_kelvin_exterior, kelvin_solution_to_physical

    mesh = _kelvin_mesh(maxh=0.34)
    with pytest.raises(ValueError, match="EXTERIOR"):
        kelvin_solution_to_physical(
            np.zeros((1, 3)), [[0.1, 0.0, 0.0]], kelvin_center=OFFSET,
            physical_center=(0.0, 0.0, 0.0), radius=RADIUS,
            form="field_strength")
    with pytest.raises(ValueError, match="form must be one of"):
        kelvin_solution_to_physical(
            np.zeros((1, 3)), [[2.0, 0.0, 0.0]], kelvin_center=OFFSET,
            physical_center=(0.0, 0.0, 0.0), radius=RADIUS, form="B")


def test_every_form_degree_and_orientation_round_trips():
    """Kelvin -> real -> Kelvin must be the identity for all eight cases.

    Degree and orientation are independent axes, so the table has to be closed
    under the involution for 0, 1, 2 and 3 forms, straight and twisted alike.
    """
    from radia.kelvin_source import (
        kelvin_computational_to_physical,
        kelvin_physical_to_computational,
        kelvin_solution_to_computational,
        kelvin_solution_to_physical,
    )

    rng = np.random.default_rng(20260908)
    physical_center = np.array([0.0, 0.0, 0.0])
    kelvin_center = np.asarray(OFFSET, dtype=float)
    directions = rng.normal(size=(6, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    physical = (RADIUS / rng.uniform(0.2, 0.9, size=(6, 1))) * directions
    computational = kelvin_physical_to_computational(
        physical, kelvin_center, physical_center, RADIUS)
    assert np.allclose(
        kelvin_computational_to_physical(
            computational, kelvin_center, physical_center, RADIUS),
        physical, rtol=1e-12, atol=1e-14)

    for degree in (0, 1, 2, 3):
        for twisted in (False, True):
            form = (degree, twisted)
            shape = (6,) if degree in (0, 3) else (6, 3)
            values = rng.normal(size=shape)
            pushed = kelvin_solution_to_computational(
                values, computational, kelvin_center=kelvin_center,
                physical_center=physical_center, radius=RADIUS, form=form)
            pulled = kelvin_solution_to_physical(
                pushed, physical, kelvin_center=kelvin_center,
                physical_center=physical_center, radius=RADIUS, form=form)
            assert np.allclose(pulled.reshape(shape), values,
                               rtol=1e-11, atol=1e-13), (
                "round trip failed for degree %d twisted=%s" % (degree, twisted))


def test_named_forms_agree_with_their_degree_and_orientation():
    """The physics names must resolve to the premetric classification."""
    from radia.kelvin_source import KELVIN_FORMS, kelvin_form_degree_and_twist

    assert KELVIN_FORMS["electric_potential"] == (0, False)
    assert KELVIN_FORMS["scalar_potential"] == (0, True)
    assert KELVIN_FORMS["vector_potential"] == (1, False)
    assert KELVIN_FORMS["field_strength"] == (1, True)
    assert KELVIN_FORMS["flux_density"] == (2, False)
    assert KELVIN_FORMS["current_density"] == (2, True)
    assert KELVIN_FORMS["charge_density"] == (3, True)
    assert kelvin_form_degree_and_twist("flux_density") == (2, False)
    assert kelvin_form_degree_and_twist((1, True)) == (1, True)
    with pytest.raises(ValueError, match="degree must be"):
        kelvin_form_degree_and_twist((4, False))

