"""Fast contracts for connected cyclic pure-HEX HDiv sectors."""
import math

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from netgen.csg import Pnt  # noqa: E402
from netgen.meshing import (  # noqa: E402
    Element2D, Element3D, FaceDescriptor, IdentificationType,
    Mesh as NetgenMesh, MeshPoint,
)

import radia.vim as vim  # noqa: E402
from radia.vim._solve import _hdiv_space_with_image_constraints  # noqa: E402
from radia.vim._vim import (  # noqa: E402
    _broken_hex_face_charge_basis,
    _charge_basis_hex,
    build_charge_gram,
)
from tests._ngsolve_2606 import curve_mesh  # noqa: E402


def _connected_sector_mesh(*, identification_count=4):
    fold = 6
    inner, outer = 0.03, 0.05
    bottom, top = -0.004, 0.004
    mesh = NetgenMesh(dim=3)
    mesh.SetMaterial(1, "yoke")
    descriptors = []
    for index, name in enumerate(
            ("skin", "periodic_min", "periodic_max"), start=1):
        descriptors.append(mesh.Add(FaceDescriptor(
            surfnr=index, domin=1, domout=0, bc=index)))
        mesh.SetBCName(index - 1, name)
    raw_points = []
    for angle in (0.0, 2.0 * math.pi / fold):
        for radius in (inner, outer):
            for z_value in (bottom, top):
                raw_points.append(mesh.Add(MeshPoint(Pnt(
                    radius * math.cos(angle), radius * math.sin(angle),
                    z_value))))
    points = [raw_points[index] for index in (0, 2, 6, 4, 1, 3, 7, 5)]
    mesh.Add(Element3D(1, points))
    for descriptor, face in (
            (descriptors[0], (0, 3, 2, 1)),
            (descriptors[0], (4, 5, 6, 7)),
            (descriptors[0], (1, 2, 6, 5)),
            (descriptors[0], (3, 0, 4, 7)),
            (descriptors[1], (0, 1, 5, 4)),
            (descriptors[2], (2, 3, 7, 6))):
        mesh.Add(Element2D(
            descriptor, [points[index] for index in face]))
    for master, slave in list(zip(
            raw_points[:4], raw_points[4:]))[:identification_count]:
        mesh.AddPointIdentification(
            master, slave, identnr=1, type=IdentificationType.PERIODIC)
    return ng.Mesh(mesh)


def _curve_cyclic_sector_q2(mesh, *, equivariant=True):
    """Install a genuine Q2 deformation while preserving the sector rotation."""
    curve_mesh(mesh, 2)
    space = ng.VectorH1(mesh, order=2)
    deformation = ng.GridFunction(space)
    radius_squared = ng.x * ng.x + ng.y * ng.y
    coefficient = (
        ng.CF((2.5 * ng.x * radius_squared,
               2.5 * ng.y * radius_squared,
               1.5 * ng.z * radius_squared))
        if equivariant else ng.CF((2.5 * ng.x * ng.x, 0.0, 0.0)))
    deformation.Interpolate(coefficient)
    mesh.SetDeformation(deformation)
    return deformation


@pytest.mark.parametrize(
    "order, expected_ndof, expected_slaves, expected_charges",
    ((1, 32, 4, 24), (2, 99, 9, 63)),
)
def test_connected_cyclic_hex_compresses_trace_and_removes_seam_charge(
        order, expected_ndof, expected_slaves, expected_charges):
    mesh = _connected_sector_mesh()
    with ng.TaskManager():
        fes, constrained, slave_dofs = _hdiv_space_with_image_constraints(
            mesh, order, (), ("periodic_min", "periodic_max"))
        charge = _charge_basis_hex(
            fes, cob_quad=max(3, order + 1), materialize_mass=False,
            excluded_boundaries=("periodic_min", "periodic_max"))

    assert constrained == ()
    assert fes.ndof == expected_ndof
    assert len(slave_dofs) == expected_slaves
    assert charge["B"].shape == (expected_charges, expected_ndof)
    assert charge["n_bf"] == 4


def test_conforming_fem_curved_q2_keeps_periodic_trace_and_curved_skin():
    mesh = _connected_sector_mesh()
    _deformation = _curve_cyclic_sector_q2(mesh)
    try:
        with ng.TaskManager():
            fes, constrained, slave_dofs = _hdiv_space_with_image_constraints(
                mesh, 2, (), ("periodic_min", "periodic_max"))
            charge = _charge_basis_hex(
                fes, cob_quad=3, materialize_mass=False,
                excluded_boundaries=("periodic_min", "periodic_max"))
    finally:
        mesh.UnsetDeformation()

    nodes = np.asarray(charge["face_nodes"], dtype=float).reshape(-1, 9, 3)
    bilinear_centres = 0.25 * (
        nodes[:, 0] + nodes[:, 2] + nodes[:, 6] + nodes[:, 8])
    assert constrained == ()
    assert len(slave_dofs) == 9
    assert fes.ndof == 99
    assert charge["B"].shape == (63, 99)
    assert charge["n_bf"] == 4
    assert np.max(np.linalg.norm(
        nodes[:, 4] - bilinear_centres, axis=1)) > 1.0e-6


def test_labeled_cyclic_seam_cannot_silently_use_images_only():
    mesh = _connected_sector_mesh()
    with ng.TaskManager(), pytest.raises(
            ValueError, match="cyclic_periodic_boundaries"):
        vim.Solve(
            mesh, mu_r=100.0, H_ext=ng.CF((0.0, 0.0, 1.0e5)),
            order=1, image_cyclic=6)


def test_cyclic_labels_require_identifications_between_the_named_faces():
    mesh = _connected_sector_mesh(identification_count=1)
    with ng.TaskManager(), pytest.raises(
            ValueError, match="pair every vertex"):
        _hdiv_space_with_image_constraints(
            mesh, 1, (), ("periodic_min", "periodic_max"))


def test_connected_cyclic_hdiv_solver_reuses_the_periodic_operator():
    mesh = _connected_sector_mesh()
    solver = vim.HDivSolver(
        mesh, order=1, image_cyclic=6,
        cyclic_periodic_boundaries=("periodic_min", "periodic_max"),
        gram_eps=1.0e-10)
    with ng.TaskManager():
        first = solver.Solve(
            mu_r=100.0, H_ext=ng.CF((0.0, 0.0, 1.0e5)), tol=1.0e-10)
        second = solver.Solve(
            mu_r=80.0, H_ext=ng.CF((0.0, 0.0, 8.0e4)), tol=1.0e-10)

    assert solver.operator_build_count == 1
    assert first["periodic_slave_dofs"] == 4
    assert second["operator_reused"]
    assert second["cyclic_periodic_boundaries"] == (
        "periodic_min", "periodic_max")


@pytest.mark.parametrize(
    "order, raw_rows, paired_rows, seam_tolerance",
    ((1, 24, 20, 1.0e-15), (2, 54, 45, 1.0e-14)),
)
def test_broken_vim_pairs_periodic_charge_jump_on_trapezoid_hex(
        order, raw_rows, paired_rows, seam_tolerance):
    mesh = _connected_sector_mesh()
    fes = ng.HDiv(mesh, order=order, discontinuous=True)
    radius = ng.sqrt(ng.x * ng.x + ng.y * ng.y)
    tangent = ng.CF((-ng.y / radius, ng.x / radius, 0.0))
    field = ng.GridFunction(fes)
    with ng.TaskManager():
        field.Set(tangent)
        raw = _broken_hex_face_charge_basis(fes, order)
        paired = _broken_hex_face_charge_basis(
            fes, order,
            cyclic_periodic_boundaries=("periodic_min", "periodic_max"),
            image_rot_angle=(2.0 * math.pi / 6.0,))

    assert raw["B"].shape == (raw_rows, fes.ndof)
    assert paired["B"].shape == (paired_rows, fes.ndof)
    assert paired["periodic_face_pair_count"] == 1
    block_size = (order + 1) ** 2
    block = paired["facet_numbers"].index(
        paired["periodic_master_facets"][0])
    charge = np.asarray(
        paired["B"] @ field.vec.FV().NumPy(), dtype=float).reshape(-1)
    seam = charge[block * block_size:(block + 1) * block_size]
    assert np.max(np.abs(seam)) < seam_tolerance


@pytest.mark.parametrize("order, expected_rows", ((1, 28), (2, 72)))
def test_broken_vim_public_chargegram_builds_cyclic_face_pair(
        order, expected_rows):
    mesh = _connected_sector_mesh()
    fes = ng.HDiv(mesh, order=order, discontinuous=True)
    fold = 6
    with ng.TaskManager():
        charge, gram, _ = build_charge_gram(
            fes, eps=1.0e-10, leafsize=32,
            internal_interfaces=True,
            image_masks=(0,) * (fold - 1),
            image_signs=(1.0,) * (fold - 1),
            image_rot_angle=tuple(
                2.0 * math.pi * index / fold
                for index in range(1, fold)),
            cyclic_periodic_boundaries=("periodic_min", "periodic_max"),
            _build_hmatrix=False)

    assert charge.shape == (expected_rows, fes.ndof)
    assert gram.ndof() == expected_rows
    assert build_charge_gram.last_timings[
        "charge_basis_periodic_face_pair_count"] == 1


def test_broken_vim_curved_q2_periodic_faces_use_ngsolve_geometry():
    mesh = _connected_sector_mesh()
    _deformation = _curve_cyclic_sector_q2(mesh)
    fes = ng.HDiv(mesh, order=2, discontinuous=True)
    radius = ng.sqrt(ng.x * ng.x + ng.y * ng.y)
    tangent = ng.CF((-ng.y / radius, ng.x / radius, 0.0))
    field = ng.GridFunction(fes)
    fold = 6
    try:
        with ng.TaskManager():
            field.Set(tangent)
            paired = _broken_hex_face_charge_basis(
                fes, 2,
                cyclic_periodic_boundaries=(
                    "periodic_min", "periodic_max"),
                image_rot_angle=(2.0 * math.pi / fold,))
            charge, gram, _ = build_charge_gram(
                fes, eps=1.0e-10, leafsize=32,
                internal_interfaces=True,
                image_masks=(0,) * (fold - 1),
                image_signs=(1.0,) * (fold - 1),
                image_rot_angle=tuple(
                    2.0 * math.pi * index / fold
                    for index in range(1, fold)),
                cyclic_periodic_boundaries=(
                    "periodic_min", "periodic_max"),
                _build_hmatrix=False)
    finally:
        mesh.UnsetDeformation()

    nodes = np.asarray(paired["face_nodes"], dtype=float)
    bilinear_centres = 0.25 * (
        nodes[:, 0] + nodes[:, 2] + nodes[:, 6] + nodes[:, 8])
    assert np.max(np.linalg.norm(
        nodes[:, 4] - bilinear_centres, axis=1)) > 1.0e-6
    assert paired["periodic_geometry_relative_residual"] < 1.0e-12
    block = paired["facet_numbers"].index(
        paired["periodic_master_facets"][0])
    seam = np.asarray(
        paired["B"] @ field.vec.FV().NumPy(), dtype=float).reshape(-1)[
            block * 9:(block + 1) * 9]
    assert np.max(np.abs(seam)) < 1.0e-14
    assert charge.shape == (72, fes.ndof)
    assert gram.ndof() == 72
    assert build_charge_gram.last_timings[
        "charge_basis_periodic_geometry_relative_residual"] < 1.0e-12


def test_broken_vim_curved_q2_rejects_nonperiodic_geometry():
    mesh = _connected_sector_mesh()
    _deformation = _curve_cyclic_sector_q2(mesh, equivariant=False)
    fes = ng.HDiv(mesh, order=2, discontinuous=True)
    try:
        with ng.TaskManager(), pytest.raises(
                ValueError, match="cyclic curved HEX facets"):
            _broken_hex_face_charge_basis(
                fes, 2,
                cyclic_periodic_boundaries=(
                    "periodic_min", "periodic_max"),
                image_rot_angle=(2.0 * math.pi / 6.0,))
    finally:
        mesh.UnsetDeformation()
