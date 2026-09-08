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
           kelvin_interface="kelvin_int"):
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
            kelvin_source_potential=kelvin_source["potential"],
            total_source_h=None, total_source_materials=())
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
    from radia.kelvin_source import kelvin_pullback_vector

    mesh = _kelvin_mesh(maxh=0.20)
    coil = _coil()
    result = _solve(mesh, coil, mu_r_iron=1.0)
    centre = np.asarray(OFFSET, dtype=float)
    errors = []
    for direction in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                      (0.6, -0.6, 0.52)):
        unit = np.asarray(direction, dtype=float)
        unit = unit / np.linalg.norm(unit)
        computational = centre + 0.62 * RADIUS * unit
        physical = (RADIUS * RADIUS / (0.62 * RADIUS)) * unit
        h_comp = np.asarray(
            result["H_cf"](mesh(*map(float, computational))), dtype=float)
        # Twisted 1-form: h = -k* h', so undoing the pullback of the stored
        # computational field carries the same minus.
        h_physical = -kelvin_pullback_vector(
            h_comp, physical, np.zeros(3), RADIUS)
        expected = np.asarray(
            rad.Fld(coil, 'h', [physical.tolist()]), dtype=float).reshape(3)
        errors.append(float(np.linalg.norm(h_physical - expected)
                            / np.linalg.norm(expected)))
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


def test_kelvin_material_without_its_interface_is_rejected():
    """Silently uncoupling the exterior is a fail-loud configuration error."""
    mesh = _kelvin_mesh()
    with pytest.raises(ValueError, match="kelvin_interface_boundary"):
        _solve(mesh, _coil(), kelvin_interface=None)
