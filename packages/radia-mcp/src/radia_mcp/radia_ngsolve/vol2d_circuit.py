"""Netgen dimension-2 ``.vol`` assembly for readable field/circuit studies.

This module deliberately stays separate from :mod:`netgen_vol`, whose public
contract is a three-dimensional tri/tet FEM/BEM mesh.  Netgen stores two-
dimensional triangles and quadrilaterals in ``surfaceelements``; here they are
domain cells and may be handed to standard NGSolve H1 or Radia H1Henrotte.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np

from .circuit_system import analyze_circuit_field


MESH_SCHEMA = "radia.netgen-vol-2d.v1"
RESULT_SCHEMA = "radia.vol2d-circuit-analysis.v1"
MU0 = 4.0e-7 * math.pi
_FAMILIES = {"P1", "Q1", "P2", "Q2", "P2_curved", "Q2_curved"}
_MAX_MCP_POINTS = 4096
_MAX_MCP_CELLS = 8192


def _package_version() -> str:
    try:
        return version("radia-mcp")
    except PackageNotFoundError:
        return "source-tree"


@dataclass(frozen=True)
class Netgen2DCell:
    material_number: int
    nodes: tuple[int, ...]

    @property
    def topology(self) -> str:
        return "triangle" if len(self.nodes) == 3 else "quadrilateral"


@dataclass(frozen=True)
class Netgen2DBoundaryEdge:
    boundary_number: int
    nodes: tuple[int, int]


@dataclass(frozen=True)
class Netgen2DVolMesh:
    points: tuple[tuple[float, float, float], ...]
    cells: tuple[Netgen2DCell, ...]
    boundary_edges: tuple[Netgen2DBoundaryEdge, ...]
    materials: dict[int, str]
    boundary_names: dict[int, str]
    content_sha256: str
    has_curved_geometry: bool
    source_name: str | None = None

    def material_name(self, number: int) -> str:
        return self.materials.get(number, f"material_{number}")

    def cell_area(self, cell: Netgen2DCell) -> float:
        xy = [(self.points[node - 1][0], self.points[node - 1][1]) for node in cell.nodes]
        twice = sum(
            x0 * y1 - x1 * y0
            for (x0, y0), (x1, y1) in zip(xy, (*xy[1:], xy[0]), strict=True)
        )
        return 0.5 * abs(twice)

    def axis_aligned_quads(self) -> bool:
        quads = [cell for cell in self.cells if len(cell.nodes) == 4]
        for cell in quads:
            xy = [(self.points[node - 1][0], self.points[node - 1][1]) for node in cell.nodes]
            scale = max(1.0, *(abs(value) for point in xy for value in point))
            tol = 1.0e-12 * scale
            xs = _clustered_values([point[0] for point in xy], tol)
            ys = _clustered_values([point[1] for point in xy], tol)
            if len(xs) != 2 or len(ys) != 2:
                return False
        return True

    def contract(self) -> dict[str, Any]:
        topology = Counter(cell.topology for cell in self.cells)
        material_counts: Counter[str] = Counter()
        material_areas: Counter[str] = Counter()
        for cell in self.cells:
            name = self.material_name(cell.material_number)
            material_counts[name] += 1
            material_areas[name] += self.cell_area(cell)
        contract = {
            "schema": MESH_SCHEMA,
            "sha256": self.content_sha256,
            "dimension": 2,
            "points": len(self.points),
            "cells": len(self.cells),
            "triangles": topology["triangle"],
            "quadrilaterals": topology["quadrilateral"],
            "boundary_edges": len(self.boundary_edges),
            "material_names": sorted(material_counts),
            "material_cell_counts": dict(sorted(material_counts.items())),
            "material_areas_m2": {
                key: float(value) for key, value in sorted(material_areas.items())
            },
            "boundary_names": sorted(set(self.boundary_names.values())),
            "has_curved_geometry": self.has_curved_geometry,
            "axis_aligned_quads": self.axis_aligned_quads(),
            "source_name": self.source_name,
        }
        canonical = {key: value for key, value in contract.items() if key != "source_name"}
        contract["contract_sha256"] = hashlib.sha256(
            json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        return contract


def _clustered_values(values: Sequence[float], tolerance: float) -> list[float]:
    groups: list[float] = []
    for value in sorted(values):
        if not groups or abs(value - groups[-1]) > tolerance:
            groups.append(value)
    return groups


def parse_netgen_2d_vol(text: str, *, source_name: str | None = None) -> Netgen2DVolMesh:
    """Parse Netgen's dimension-2 tri/quad topology without importing NGSolve."""

    if not isinstance(text, str) or not text.strip():
        raise ValueError("vol_text must be a non-empty string")
    encoded = text.encode("utf-8")
    if len(encoded) > 8 * 1024 * 1024:
        raise ValueError("vol_text exceeds the 8 MiB MCP assembly limit")
    lines = [line.strip() for line in text.splitlines()]
    index = 0
    dimension: int | None = None
    points: list[tuple[float, float, float]] = []
    cells: list[Netgen2DCell] = []
    edges: list[Netgen2DBoundaryEdge] = []
    materials: dict[int, str] = {}
    boundary_names: dict[int, str] = {}
    has_curved_geometry = False

    def next_data_line() -> str | None:
        nonlocal index
        while index < len(lines):
            line = lines[index]
            index += 1
            if line and not line.startswith("#"):
                return line
        return None

    def required(context: str) -> str:
        value = next_data_line()
        if value is None:
            raise ValueError(f"unexpected end of .vol while reading {context}")
        return value

    def count(context: str) -> int:
        value = int(required(context).split()[0])
        if value < 0:
            raise ValueError(f"{context} count must be nonnegative")
        return value

    def skip_counted(context: str) -> None:
        for _ in range(count(context)):
            required(context)

    while True:
        key = next_data_line()
        if key is None or key.lower() == "endmesh":
            break
        normalized = key.lower()
        if normalized == "mesh3d":
            continue
        if normalized == "dimension":
            dimension = int(required("dimension").split()[0])
        elif normalized == "geomtype":
            required("geomtype")
        elif normalized == "facedescriptors":
            skip_counted("facedescriptors")
        elif normalized in {"surfaceelements", "surfaceelementsuv"}:
            for cell_index in range(1, count("surfaceelements") + 1):
                record = required("surfaceelements").split()
                if len(record) < 6:
                    raise ValueError(f"surface element {cell_index} is too short")
                head = [int(value) for value in record[:5]]
                node_count = head[4]
                nodes = tuple(int(value) for value in record[5 : 5 + node_count])
                if len(nodes) != node_count:
                    raise ValueError(f"surface element {cell_index} node count is truncated")
                if node_count not in {3, 4}:
                    raise ValueError(
                        f"dimension-2 .vol supports triangle/quad cells; element "
                        f"{cell_index} has {node_count} nodes"
                    )
                # In a dimension-2 Netgen .vol, the second field is the domain
                # material number used by Mesh.GetMaterials().
                cells.append(Netgen2DCell(head[1], nodes))
        elif normalized == "volumeelements":
            volume_count = count("volumeelements")
            if volume_count:
                raise ValueError("dimension-2 .vol must not contain volume elements")
        elif normalized in {"edgesegmentsgi2", "edgesegmentsgi3"}:
            for edge_index in range(1, count(normalized) + 1):
                record = required(normalized).split()
                minimum_fields = 4 if normalized == "edgesegmentsgi2" else 7
                if len(record) < minimum_fields:
                    raise ValueError(f"boundary edge {edge_index} is too short")
                if normalized == "edgesegmentsgi2":
                    boundary_number = int(record[0])
                    nodes = (int(record[2]), int(record[3]))
                else:
                    # Netgen 6.2.2606 stores p1/p2 first and a zero-based
                    # EdgeDescriptor table index (edsi) in field seven.
                    boundary_number = int(record[6]) + 1
                    nodes = (int(record[0]), int(record[1]))
                edges.append(
                    Netgen2DBoundaryEdge(boundary_number, nodes)
                )
        elif normalized == "points":
            for point_index in range(1, count("points") + 1):
                values = [float(value) for value in required("points").split()]
                if len(values) < 3 or not all(math.isfinite(value) for value in values[:3]):
                    raise ValueError(f"point {point_index} must contain finite x y z")
                points.append((values[0], values[1], values[2]))
        elif normalized == "materials":
            for _ in range(count("materials")):
                record = required("materials").split(maxsplit=1)
                materials[int(record[0])] = (
                    record[1] if len(record) > 1 else f"material_{record[0]}"
                )
        elif normalized == "bcnames":
            for _ in range(count("bcnames")):
                record = required("bcnames").split(maxsplit=1)
                boundary_names[int(record[0])] = (
                    record[1] if len(record) > 1 else f"boundary_{record[0]}"
                )
        elif normalized in {
            "pointelements",
            "edgedescriptors",
            "face_colours",
            "face_transparencies",
        }:
            skip_counted(normalized)
        elif normalized == "curvedelements":
            has_curved_geometry = True
            break
        else:
            raise ValueError(f"unsupported or unexpected dimension-2 .vol section: {key}")

    if dimension != 2:
        raise ValueError(f"2-D assembly requires dimension 2, got {dimension!r}")
    if not points or not cells:
        raise ValueError("dimension-2 .vol requires non-empty points and cells")
    if any(abs(point[2]) > 1.0e-12 for point in points):
        raise ValueError("dimension-2 .vol points must lie in the z=0 storage plane")
    for label, records in (("cell", cells), ("boundary edge", edges)):
        for record_index, record in enumerate(records, start=1):
            for node in record.nodes:
                if node < 1 or node > len(points):
                    raise ValueError(
                        f"{label} {record_index} references node {node}, point count {len(points)}"
                    )
    mesh = Netgen2DVolMesh(
        points=tuple(points),
        cells=tuple(cells),
        boundary_edges=tuple(edges),
        materials=materials,
        boundary_names=boundary_names,
        content_sha256=hashlib.sha256(encoded).hexdigest(),
        has_curved_geometry=has_curved_geometry,
        source_name=source_name,
    )
    for cell_index, cell in enumerate(mesh.cells, start=1):
        if mesh.cell_area(cell) <= 1.0e-18:
            raise ValueError(f"cell {cell_index} has zero or negligible area")
    return mesh


def inspect_netgen_2d_vol(text: str, *, source_name: str | None = None) -> dict[str, Any]:
    """Return the neutral contract consumed by public and private MCP lanes."""

    return parse_netgen_2d_vol(text, source_name=source_name).contract()


def _family_contract(mesh: Netgen2DVolMesh, family: str, formulation: str) -> dict[str, Any]:
    if family not in _FAMILIES:
        raise ValueError(f"element_family must be one of {sorted(_FAMILIES)}")
    triangles = sum(len(cell.nodes) == 3 for cell in mesh.cells)
    quads = sum(len(cell.nodes) == 4 for cell in mesh.cells)
    if triangles and quads:
        raise ValueError("mixed triangle/quad dimension-2 assembly is not supported")
    if family.startswith("P") and not triangles:
        raise ValueError(f"{family} requires an all-triangle .vol")
    if family.startswith("Q") and not quads:
        raise ValueError(f"{family} requires an all-quad .vol")
    if family == "P2_curved" and not mesh.has_curved_geometry:
        raise ValueError("P2_curved requires curvedelements in the .vol")
    if formulation == "axisymmetric_henrotte" and family in {"Q1", "Q2"}:
        if not mesh.axis_aligned_quads():
            raise ValueError(f"{family} axisymmetric Henrotte requires axis-aligned quads")
    return {
        "order": 1 if family in {"P1", "Q1"} else 2,
        "curvedquad": family == "Q2_curved",
        "curved_geometry": family in {"P2_curved", "Q2_curved"},
    }


def _finite_positive(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be finite and positive") from exc
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{label} must be finite and positive")
    return result


def _material_permeability(mesh_contract: dict[str, Any], raw: Any) -> dict[str, float]:
    names = list(mesh_contract["material_names"])
    if isinstance(raw, Mapping):
        missing = sorted(set(names) - {str(key) for key in raw})
        unknown = sorted({str(key) for key in raw} - set(names))
        if missing or unknown:
            raise ValueError(
                f"permeability_h_per_m material mismatch: missing={missing}, unknown={unknown}"
            )
        return {
            name: _finite_positive(raw[name], f"permeability_h_per_m[{name}]")
            for name in names
        }
    value = _finite_positive(MU0 if raw is None else raw, "permeability_h_per_m")
    return {name: value for name in names}


def _runtime_vol_path(text: str, digest: str) -> Path:
    root = Path(os.environ.get("RADIA_MCP_TEMP", r"C:\temp")) / "radia_mcp_vol2d"
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{digest}.vol"
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    if not path.exists() or path.read_text(encoding="utf-8") != normalized:
        path.write_text(normalized, encoding="utf-8", newline="\n")
    return path


def _dense_matrix(matrix: Any) -> np.ndarray:
    rows, columns, values = matrix.COO()
    dense = np.zeros((matrix.height, matrix.width), dtype=float)
    np.add.at(
        dense,
        (np.asarray(rows, dtype=int), np.asarray(columns, dtype=int)),
        np.asarray(values, dtype=float),
    )
    return dense


def _validate_ngsolve_2d_execution_text(text: str) -> None:
    """Reject inventory-only shorthand before it reaches native NGSolve code.

    ``inspect_netgen_2d_vol`` intentionally accepts compact educational text,
    but native assembly needs the face descriptors and complete edge records
    written by Netgen.  Passing shorthand through can terminate the process in
    native code instead of raising a Python exception.
    """

    lines = [line.strip() for line in text.splitlines()]
    lowered = [line.lower() for line in lines]
    if "facedescriptors" not in lowered:
        raise ValueError(
            "solver-ready dimension-2 .vol requires a facedescriptors section; "
            "generate it with the Netgen writer"
        )
    edge_section = next(
        (name for name in ("edgesegmentsgi3", "edgesegmentsgi2") if name in lowered),
        None,
    )
    if edge_section is None:
        raise ValueError("solver-ready dimension-2 .vol has no edge segment section")
    try:
        edge_index = lowered.index(edge_section)
        edge_count = int(lines[edge_index + 1])
    except (ValueError, IndexError) as exc:
        raise ValueError(
            f"solver-ready dimension-2 .vol has invalid {edge_section}"
        ) from exc
    edge_rows: list[list[str]] = []
    for line in lines[edge_index + 2 :]:
        if not line or line.startswith("#"):
            continue
        edge_rows.append(line.split())
        if len(edge_rows) == edge_count:
            break
    minimum_fields = 7 if edge_section == "edgesegmentsgi3" else 12
    if len(edge_rows) != edge_count or any(
        len(row) < minimum_fields for row in edge_rows
    ):
        raise ValueError(
            "solver-ready dimension-2 .vol requires complete Netgen edge records"
        )
    if edge_section == "edgesegmentsgi3" and "edgedescriptors" not in lowered:
        raise ValueError(
            "solver-ready dimension-2 .vol requires edgedescriptors with edgesegmentsgi3"
        )


def assemble_vol2d_field(request: Mapping[str, Any]) -> dict[str, Any]:
    """Assemble one constrained global K and signed branch-source matrix."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    text = request.get("vol_text")
    if not isinstance(text, str):
        raise ValueError("vol_text must be text")
    source_name = str(request.get("source_name", "generated.vol"))
    mesh_view = parse_netgen_2d_vol(text, source_name=source_name)
    _validate_ngsolve_2d_execution_text(text)
    mesh_contract = mesh_view.contract()
    if mesh_contract["points"] > _MAX_MCP_POINTS or mesh_contract["cells"] > _MAX_MCP_CELLS:
        raise ValueError(
            "dimension-2 MCP assembly is limited to "
            f"{_MAX_MCP_POINTS} points and {_MAX_MCP_CELLS} cells"
        )
    family = str(request.get("element_family", ""))
    formulation = str(request.get("formulation", ""))
    if formulation not in {"planar", "axisymmetric_henrotte"}:
        raise ValueError("formulation must be planar or axisymmetric_henrotte")
    family_options = _family_contract(mesh_view, family, formulation)
    boundaries = request.get("dirichlet_boundaries")
    if not isinstance(boundaries, Sequence) or isinstance(boundaries, (str, bytes)):
        raise ValueError("dirichlet_boundaries must be a non-empty string sequence")
    boundary_names = [str(value) for value in boundaries]
    if not boundary_names or any(not value for value in boundary_names):
        raise ValueError("dirichlet_boundaries must be non-empty")
    missing_boundaries = sorted(set(boundary_names) - set(mesh_contract["boundary_names"]))
    if missing_boundaries:
        raise ValueError(f"Dirichlet boundaries are absent from .vol: {missing_boundaries}")

    permeability = _material_permeability(
        mesh_contract, request.get("permeability_h_per_m", MU0)
    )
    branches = request.get("branches")
    if not isinstance(branches, Sequence) or isinstance(branches, (str, bytes)) or not branches:
        raise ValueError("branches must be a non-empty sequence")
    parsed_branches: list[dict[str, Any]] = []
    branch_names: set[str] = set()
    for index, raw in enumerate(branches):
        if not isinstance(raw, Mapping):
            raise ValueError(f"branches[{index}] must be an object")
        name = str(raw.get("name", "")).strip()
        material = str(raw.get("material", "")).strip()
        turns = float(raw.get("turns"))
        if not name or name in branch_names:
            raise ValueError("branch names must be non-empty and unique")
        if material not in mesh_contract["material_names"]:
            raise ValueError(f"branch {name} material is absent from .vol: {material}")
        if not math.isfinite(turns) or turns == 0.0:
            raise ValueError(f"branch {name} turns must be finite and non-zero")
        branch_names.add(name)
        parsed_branches.append({"name": name, "material": material, "turns": turns})

    try:
        from ngsolve import (  # type: ignore
            BilinearForm,
            CoefficientFunction,
            H1,
            LinearForm,
            Mesh,
            dx,
            grad,
            x,
        )
    except ImportError as exc:
        raise RuntimeError("NGSolve is required for dimension-2 .vol assembly") from exc

    runtime_path = _runtime_vol_path(text, mesh_contract["sha256"])
    mesh = Mesh(str(runtime_path))
    dirichlet = "|".join(boundary_names)
    mu_cf = mesh.MaterialCF(permeability)
    if formulation == "planar":
        fes = H1(mesh, order=family_options["order"], dirichlet=dirichlet)
        trial, test = fes.TnT()
        stiffness = BilinearForm(fes, symmetric=True)
        stiffness += grad(trial) * grad(test) / mu_cf * dx
        backend = "ngsolve-h1"
    else:
        try:
            from radia.axifem import H1Henrotte, AxiHenrotteStiffnessBFI
        except ImportError as exc:
            raise RuntimeError("Radia axifem extension is required") from exc
        fes = H1Henrotte(
            mesh,
            order=family_options["order"],
            dirichlet=dirichlet,
            curvedquad=family_options["curvedquad"],
        )
        test = fes.TestFunction()
        stiffness = BilinearForm(fes, symmetric=True)
        stiffness += AxiHenrotteStiffnessBFI(mu_cf)
        backend = "radia-axifem-h1henrotte"

    linear_forms: list[Any] = []
    branch_areas: list[float] = []
    for branch in parsed_branches:
        area = float(mesh_contract["material_areas_m2"][branch["material"]])
        if area <= 0.0:
            raise ValueError(f"branch {branch['name']} material has zero area")
        current_density_per_amp = mesh.MaterialCF(
            {branch["material"]: branch["turns"] / area}, default=0.0
        )
        form = LinearForm(fes)
        source_integrand = current_density_per_amp * test
        if formulation == "axisymmetric_henrotte":
            source_integrand = 2.0 * math.pi * x * source_integrand
        form += source_integrand * dx
        linear_forms.append(form)
        branch_areas.append(area)

    # FastMCP invokes synchronous tools in a worker thread.  Opening NGSolve's
    # global TaskManager from that worker can deadlock the stdio request.  This
    # MCP path is deliberately capped at a small dense system, so serial
    # assembly is the deterministic protocol-safe choice.
    stiffness.Assemble()
    for form in linear_forms:
        form.Assemble()
    free = np.flatnonzero(np.asarray(fes.FreeDofs(), dtype=bool))
    if free.size == 0:
        raise ValueError("Dirichlet constraints leave no free field degrees of freedom")
    maximum_dense_dofs = int(request.get("maximum_dense_dofs", 512))
    if maximum_dense_dofs <= 0 or free.size > maximum_dense_dofs:
        raise ValueError(
            f"dense MCP solve is limited to {maximum_dense_dofs} free DOFs; got {free.size}"
        )
    full_matrix = _dense_matrix(stiffness.mat)
    field_matrix = full_matrix[np.ix_(free, free)]
    source_matrix = np.column_stack(
        [np.asarray(form.vec.FV().NumPy(), dtype=float)[free] for form in linear_forms]
    )
    symmetry_error = float(np.max(np.abs(field_matrix - field_matrix.T)))
    if symmetry_error > 1.0e-9 * max(1.0, float(np.max(np.abs(field_matrix)))):
        raise ValueError(f"assembled field matrix is not symmetric: {symmetry_error:.3e}")
    try:
        condition_number = float(np.linalg.cond(field_matrix))
        np.linalg.cholesky(0.5 * (field_matrix + field_matrix.T))
    except np.linalg.LinAlgError as exc:
        raise ValueError("constrained field matrix is not positive definite") from exc

    return {
        "schema": RESULT_SCHEMA,
        "status": "assembled",
        "backend": backend,
        "formulation": formulation,
        "element_family": family,
        "family_options": family_options,
        "mesh_contract": mesh_contract,
        "field_dofs": int(fes.ndof),
        "free_field_dofs": int(free.size),
        "free_dof_indices_0based": free.tolist(),
        "branch_order": [branch["name"] for branch in parsed_branches],
        "branch_materials": [branch["material"] for branch in parsed_branches],
        "branch_turns": [float(branch["turns"]) for branch in parsed_branches],
        "branch_area_m2": branch_areas,
        "field_matrix": field_matrix.tolist(),
        "source_matrix": source_matrix.tolist(),
        "field_rhs": [0.0] * int(free.size),
        "matrix_symmetry_error": symmetry_error,
        "matrix_condition_number": condition_number,
        "source_column_norms": [
            float(np.linalg.norm(source_matrix[:, index]))
            for index in range(source_matrix.shape[1])
        ],
        "generated_mesh_git_required": False,
    }


def analyze_vol2d_circuit(request: Mapping[str, Any]) -> dict[str, Any]:
    """Compose dimension-2 assembly with the v2 field-circuit solver."""

    assembly = assemble_vol2d_field(request)
    circuit = request.get("circuit")
    if not isinstance(circuit, Mapping):
        raise ValueError("circuit must be an object")
    circuit_payload = dict(circuit)
    circuit_payload.update(
        {
            "field_matrix": assembly["field_matrix"],
            "source_matrix": assembly["source_matrix"],
            "field_rhs": assembly["field_rhs"],
        }
    )
    solution = analyze_circuit_field(circuit_payload)
    return {
        "execution_version": {"radia_mcp": _package_version()},
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema": RESULT_SCHEMA,
        "status": "solved",
        "assembly": assembly,
        "solution": solution,
    }


def write_structured_rect_vol(
    path: str | Path,
    *,
    x0: float = 0.0,
    x1: float = 1.0,
    y0: float = 0.0,
    y1: float = 1.0,
    nx: int = 2,
    ny: int = 2,
    quads: bool = False,
    material: str = "domain",
) -> dict[str, Any]:
    """Generate a deterministic replay mesh; the ``.vol`` need not enter git."""

    if nx < 2 or ny < 2:
        raise ValueError("nx and ny must be at least two")
    if not (x1 > x0 and y1 > y0):
        raise ValueError("rectangle bounds must be strictly increasing")
    from netgen.meshing import (  # type: ignore
        EdgeDescriptor,
        Element1D,
        Element2D,
        FaceDescriptor,
        Mesh as NgMesh,
        MeshPoint,
        Pnt,
    )

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    mesh = NgMesh()
    mesh.dim = 2
    mesh.SetMaterial(1, material)
    for boundary, name in enumerate(("bottom", "right", "top", "left"), start=1):
        mesh.Add(FaceDescriptor(surfnr=boundary, domin=0, bc=boundary))
        mesh.SetBCName(boundary - 1, name)
        edge = EdgeDescriptor()
        edge.edgenr = boundary
        edge.surfnr = (boundary, -1)
        edge.domin = 1
        edge.domout = 0
        edge.name = name
        mesh.Add(edge)
    point_ids: list[list[Any]] = []
    for j in range(ny + 1):
        row = []
        y_value = y0 + (y1 - y0) * j / ny
        for i in range(nx + 1):
            x_value = x0 + (x1 - x0) * i / nx
            row.append(mesh.Add(MeshPoint(Pnt(x_value, y_value, 0.0))))
        point_ids.append(row)
    for j in range(ny):
        for i in range(nx):
            p00, p10 = point_ids[j][i], point_ids[j][i + 1]
            p01, p11 = point_ids[j + 1][i], point_ids[j + 1][i + 1]
            if quads:
                mesh.Add(Element2D(1, [p00, p10, p11, p01]))
            else:
                mesh.Add(Element2D(1, [p00, p10, p11]))
                mesh.Add(Element2D(1, [p00, p11, p01]))
    for i in range(nx):
        mesh.Add(Element1D([point_ids[0][i], point_ids[0][i + 1]], index=1))
        mesh.Add(Element1D([point_ids[ny][i], point_ids[ny][i + 1]], index=3))
    for j in range(ny):
        mesh.Add(Element1D([point_ids[j][nx], point_ids[j + 1][nx]], index=2))
        mesh.Add(Element1D([point_ids[j][0], point_ids[j + 1][0]], index=4))
    mesh.Save(str(target))
    return inspect_netgen_2d_vol(target.read_text(encoding="utf-8"), source_name=target.name)


def write_structured_material_rect_vol(
    path: str | Path,
    *,
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    nx: int,
    ny: int,
    rectangles: Sequence[Mapping[str, Any]],
    background_material: str = "air",
    quads: bool = False,
) -> dict[str, Any]:
    """Generate a deterministic multi-material replay mesh.

    Rectangles are assigned by cell centroid and therefore must align with the
    structured grid when exact interfaces are required.  The generated
    ``.vol`` is an execution artifact; callers need only version the script.
    """

    if nx < 2 or ny < 2:
        raise ValueError("nx and ny must be at least two")
    if not (x1 > x0 and y1 > y0):
        raise ValueError("rectangle bounds must be strictly increasing")
    background = str(background_material).strip()
    if not background:
        raise ValueError("background_material must be non-empty")
    parsed: list[dict[str, Any]] = []
    names = {background}
    for index, raw in enumerate(rectangles):
        if not isinstance(raw, Mapping):
            raise ValueError(f"rectangles[{index}] must be an object")
        name = str(raw.get("name", "")).strip()
        if not name or name in names:
            raise ValueError("rectangle material names must be non-empty and unique")
        bounds = tuple(
            float(raw.get(key)) for key in ("x0", "x1", "y0", "y1")
        )
        if not all(math.isfinite(value) for value in bounds):
            raise ValueError(f"rectangles[{index}] bounds must be finite")
        rx0, rx1, ry0, ry1 = bounds
        if not (x0 <= rx0 < rx1 <= x1 and y0 <= ry0 < ry1 <= y1):
            raise ValueError(f"rectangles[{index}] must lie inside the domain")
        for other in parsed:
            overlap_x = min(rx1, other["x1"]) - max(rx0, other["x0"])
            overlap_y = min(ry1, other["y1"]) - max(ry0, other["y0"])
            if overlap_x > 0.0 and overlap_y > 0.0:
                raise ValueError("material rectangles must not overlap")
        parsed.append(
            {"name": name, "x0": rx0, "x1": rx1, "y0": ry0, "y1": ry1}
        )
        names.add(name)

    from netgen.meshing import (  # type: ignore
        EdgeDescriptor,
        Element1D,
        Element2D,
        FaceDescriptor,
        Mesh as NgMesh,
        MeshPoint,
        Pnt,
    )

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    mesh = NgMesh()
    mesh.dim = 2
    material_names = [background, *(row["name"] for row in parsed)]
    for material_index, name in enumerate(material_names, start=1):
        mesh.SetMaterial(material_index, name)
    for boundary, name in enumerate(("bottom", "right", "top", "left"), start=1):
        mesh.Add(FaceDescriptor(surfnr=boundary, domin=0, bc=boundary))
        mesh.SetBCName(boundary - 1, name)
        edge = EdgeDescriptor()
        edge.edgenr = boundary
        edge.surfnr = (boundary, -1)
        edge.domin = 1
        edge.domout = 0
        edge.name = name
        mesh.Add(edge)
    point_ids: list[list[Any]] = []
    for j in range(ny + 1):
        row = []
        y_value = y0 + (y1 - y0) * j / ny
        for i in range(nx + 1):
            x_value = x0 + (x1 - x0) * i / nx
            row.append(mesh.Add(MeshPoint(Pnt(x_value, y_value, 0.0))))
        point_ids.append(row)

    def material_at(x_value: float, y_value: float) -> int:
        for material_index, rectangle in enumerate(parsed, start=2):
            if (
                rectangle["x0"] <= x_value <= rectangle["x1"]
                and rectangle["y0"] <= y_value <= rectangle["y1"]
            ):
                return material_index
        return 1

    for j in range(ny):
        for i in range(nx):
            xa = x0 + (x1 - x0) * i / nx
            xb = x0 + (x1 - x0) * (i + 1) / nx
            ya = y0 + (y1 - y0) * j / ny
            yb = y0 + (y1 - y0) * (j + 1) / ny
            material = material_at(0.5 * (xa + xb), 0.5 * (ya + yb))
            p00, p10 = point_ids[j][i], point_ids[j][i + 1]
            p01, p11 = point_ids[j + 1][i], point_ids[j + 1][i + 1]
            if quads:
                mesh.Add(Element2D(material, [p00, p10, p11, p01]))
            else:
                mesh.Add(Element2D(material, [p00, p10, p11]))
                mesh.Add(Element2D(material, [p00, p11, p01]))
    for i in range(nx):
        mesh.Add(Element1D([point_ids[0][i], point_ids[0][i + 1]], index=1))
        mesh.Add(Element1D([point_ids[ny][i], point_ids[ny][i + 1]], index=3))
    for j in range(ny):
        mesh.Add(Element1D([point_ids[j][nx], point_ids[j + 1][nx]], index=2))
        mesh.Add(Element1D([point_ids[j][0], point_ids[j + 1][0]], index=4))
    mesh.Save(str(target))
    return inspect_netgen_2d_vol(
        target.read_text(encoding="utf-8"), source_name=target.name
    )
