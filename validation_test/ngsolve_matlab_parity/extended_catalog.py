"""Deterministic catalogs and mesh fixtures for extended MATLAB parity."""

from __future__ import annotations

from typing import Any

import ngsolve as ng
from netgen.csg import Pnt
from netgen.geom2d import unit_square
from netgen.meshing import (
    Element1D,
    Element2D,
    Element3D,
    FaceDescriptor,
    Mesh as NetgenMesh,
    MeshPoint,
)
from netgen.occ import unit_cube
from ngsolve.meshes import MakeStructured3DMesh


DOF_LIMIT = 1_000_000
SCALE_DOF_FLOOR = 10_000
SMALL_MESH_IDS = (
    "tri_single", "tri_meshed", "quad_single", "quad_warped",
    "mixed_tri_quad", "tet_single", "tet_meshed", "hex_single",
    "hex_warped", "wedge_single",
)
SCALE_MESH_IDS = ("tri_scale", "tet_scale", "hex_scale", "wedge_scale")


def _single_cell(dimension: int, vertices, volume, boundary_faces):
    mesh = NetgenMesh(dim=dimension)
    mesh.SetMaterial(1, "domain")
    if dimension == 2:
        mesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    else:
        boundary = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    mesh.SetBCName(0, "outer")
    points = [mesh.Add(MeshPoint(Pnt(*point))) for point in vertices]
    if dimension == 2:
        mesh.Add(Element2D(1, [points[index] for index in volume]))
        for first, second in boundary_faces:
            mesh.Add(Element1D([points[first], points[second]], index=1))
    else:
        mesh.Add(Element3D(1, [points[index] for index in volume]))
        for face in boundary_faces:
            mesh.Add(Element2D(boundary, [points[index] for index in face]))
    mesh.Update()
    return mesh


def _mixed_tri_quad():
    mesh = NetgenMesh(dim=2)
    mesh.SetMaterial(1, "domain")
    mesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    mesh.SetBCName(0, "outer")
    coordinates = ((0, 0, 0), (1, 0, 0), (2, 0, 0),
                   (0, 1, 0), (1, 1, 0), (2, 1, 0))
    points = [mesh.Add(MeshPoint(Pnt(*point))) for point in coordinates]
    mesh.Add(Element2D(1, [points[index] for index in (0, 1, 4, 3)]))
    mesh.Add(Element2D(1, [points[index] for index in (1, 2, 5)]))
    mesh.Add(Element2D(1, [points[index] for index in (1, 5, 4)]))
    for first, second in ((0, 1), (1, 2), (2, 5),
                          (5, 4), (4, 3), (3, 0)):
        mesh.Add(Element1D([points[first], points[second]], index=1))
    mesh.Update()
    return mesh


def build_mesh(mesh_id: str):
    triangle = (
        ((0, 0, 0), (1, 0, 0), (0, 1, 0)),
        (0, 1, 2),
        ((0, 1), (1, 2), (2, 0)),
    )
    quad_faces = ((0, 1), (1, 2), (2, 3), (3, 0))
    tetra = (
        ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
        (0, 1, 2, 3),
        ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)),
    )
    hex_faces = (
        (0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
        (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
    )
    wedge_faces = (
        (0, 2, 1), (3, 4, 5), (0, 1, 4, 3),
        (1, 2, 5, 4), (2, 0, 3, 5),
    )
    if mesh_id == "tri_single":
        return _single_cell(2, *triangle)
    if mesh_id == "tri_meshed":
        return unit_square.GenerateMesh(maxh=0.55)
    if mesh_id == "quad_single":
        return _single_cell(
            2, ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)),
            (0, 1, 2, 3), quad_faces)
    if mesh_id == "quad_warped":
        return _single_cell(
            2, ((0, 0, 0), (1.1, 0, 0), (0.9, 1.0, 0), (-0.1, 0.85, 0)),
            (0, 1, 2, 3), quad_faces)
    if mesh_id == "mixed_tri_quad":
        return _mixed_tri_quad()
    if mesh_id == "tet_single":
        return _single_cell(3, *tetra)
    if mesh_id == "tet_meshed":
        return unit_cube.GenerateMesh(maxh=0.8)
    if mesh_id == "hex_single":
        return _single_cell(
            3,
            ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
             (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),
            tuple(range(8)), hex_faces)
    if mesh_id == "hex_warped":
        return _single_cell(
            3,
            ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
             (0.04, -0.02, 1), (1.03, 0.03, 1.02),
             (0.97, 1.04, 0.96), (-0.03, 0.98, 1.01)),
            tuple(range(8)), hex_faces)
    if mesh_id == "wedge_single":
        return _single_cell(
            3,
            ((0, 0, 0), (1, 0, 0), (0, 1, 0),
             (0, 0, 1), (1, 0, 1), (0, 1, 1)),
            tuple(range(6)), wedge_faces)
    if mesh_id == "tri_scale":
        return unit_square.GenerateMesh(maxh=0.025)
    if mesh_id == "tet_scale":
        return unit_cube.GenerateMesh(maxh=0.16)
    if mesh_id == "hex_scale":
        return MakeStructured3DMesh(hexes=True, nx=10, ny=10, nz=10)
    if mesh_id == "wedge_scale":
        return MakeStructured3DMesh(prism=True, nx=10, ny=10, nz=10)
    raise ValueError(f"unknown extended mesh id {mesh_id!r}")


def build_small_meshes() -> dict[str, Any]:
    return {mesh_id: build_mesh(mesh_id) for mesh_id in SMALL_MESH_IDS}


def build_scale_meshes() -> dict[str, ng.Mesh]:
    return {mesh_id: build_mesh(mesh_id) for mesh_id in SCALE_MESH_IDS}


def _operator(space: str, order: int, form: str, *, weight=1.0,
              dirichlet="", solve=True, coefficient_kind="constant",
              boundary=False, variant="base") -> dict[str, Any]:
    return {
        "space": space,
        "order": order,
        "form": form,
        "weight_real": float(complex(weight).real),
        "weight_imag": float(complex(weight).imag),
        "dirichlet": dirichlet,
        "solve": bool(solve),
        "coefficient_kind": coefficient_kind,
        "boundary": bool(boundary),
        "variant": variant,
    }


def _base_operators() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    specs.extend(_operator("h1", p, "mass") for p in range(1, 5))
    specs.extend(_operator("h1", p, "stiffness", dirichlet=".*")
                 for p in range(1, 5))
    specs.extend(_operator("hcurl", p, "mass") for p in range(1, 4))
    specs.extend(_operator("hcurl", p, "curlcurl", solve=False)
                 for p in range(1, 4))
    specs.extend(_operator("hdiv", p, "mass") for p in range(1, 4))
    specs.extend(_operator("hdiv", p, "divdiv", solve=False)
                 for p in range(1, 4))
    specs.extend((
        _operator("h1", 3, "mass", weight=0.75, variant="weighted075"),
        _operator("h1", 3, "stiffness", weight=2.0, dirichlet=".*",
                  variant="weighted200"),
        _operator("hcurl", 2, "mass", weight=1.25, variant="weighted125"),
        _operator("hdiv", 2, "mass", weight=0.5, variant="weighted050"),
        _operator("hdiv", 2, "divdiv", weight=1.5, solve=False,
                  variant="weighted150"),
    ))
    assert len(specs) == 25
    return specs


def breadth_operator_specs() -> list[dict[str, Any]]:
    base = _base_operators()
    spatial_indices = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 14, 15, 17, 20)
    spatial = [
        {**base[index], "solve": False, "coefficient_kind": "affine_x",
         "variant": f"affine_x_{index + 1:02d}"}
        for index in spatial_indices
    ]
    complex_specs = [
        _operator("h1", 2, "mass", weight=1.0 + 0.15j, solve=False,
                  variant="complex"),
        _operator("h1", 2, "stiffness", weight=1.0 + 0.15j,
                  dirichlet=".*", solve=False, variant="complex"),
        _operator("hcurl", 2, "mass", weight=1.0 + 0.15j, solve=False,
                  variant="complex"),
        _operator("hcurl", 2, "curlcurl", weight=1.0 + 0.15j,
                  solve=False, variant="complex"),
        _operator("hdiv", 2, "mass", weight=1.0 + 0.15j, solve=False,
                  variant="complex"),
        _operator("hdiv", 2, "divdiv", weight=1.0 + 0.15j,
                  solve=False, variant="complex"),
    ]
    boundary = [
        _operator("h1", p, "boundary_mass", solve=False,
                  coefficient_kind="affine_x", boundary=True,
                  variant="boundary_affine_x")
        for p in range(1, 5)
    ]
    result = base + spatial + complex_specs + boundary
    assert len(result) == 50
    return result


def build_breadth_catalog(mesh_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for mesh in mesh_rows:
        for operator in breadth_operator_specs():
            number = len(cases) + 1
            cases.append({
                "case_id": (
                    f"breadth_{number:03d}_{mesh['id']}_{operator['space']}_"
                    f"p{operator['order']}_{operator['form']}_{operator['variant']}"
                ),
                "oracle_key": f"breadth_{number:03d}",
                "mesh_id": mesh["id"],
                "mesh_path": mesh["path"],
                "dimension": mesh["dimension"],
                **operator,
            })
    assert len(cases) == 500
    return cases


def build_scale_catalog(mesh_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    operators = (
        _operator("h1", 4, "mass", solve=False),
        _operator("h1", 4, "stiffness", dirichlet=".*", solve=False),
        _operator("hcurl", 3, "mass", solve=False),
        _operator("hcurl", 3, "curlcurl", solve=False),
        _operator("hdiv", 2, "divdiv", solve=False),
    )
    cases: list[dict[str, Any]] = []
    for mesh in mesh_rows:
        for operator in operators:
            number = len(cases) + 1
            cases.append({
                "case_id": (
                    f"scale_{number:02d}_{mesh['id']}_{operator['space']}_"
                    f"p{operator['order']}_{operator['form']}"
                ),
                "oracle_key": f"scale_{number:02d}",
                "mesh_id": mesh["id"],
                "mesh_path": mesh["path"],
                "dimension": mesh["dimension"],
                **operator,
            })
    assert len(cases) == 20
    return cases
