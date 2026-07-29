"""Portable scalar-PDE analysis for dimension-2 Netgen ``.vol`` meshes.

Electrostatics, current flow, and steady heat conduction share the operator
``-div(C grad(u)) = s``.  This module gives that operator one closed-world,
SI-unit result contract for planar-depth and full-revolution axisymmetric
models.  It deliberately accepts mesh text rather than arbitrary paths and
writes visualization files only below the owned temporary directory.
"""

from __future__ import annotations

import csv
import io
import json
import math
import time
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .vol2d_circuit import (
    _family_contract,
    _runtime_vol_path,
    parse_netgen_2d_vol,
)
from .vol2d_postprocess import (
    _SAFE_NAME,
    _canonical,
    _export_entry,
    _package_version,
    _sha,
)


SCALAR_SCHEMA = "radia.vol2d-scalar-analysis.v1"
REPLAY_SCHEMA = "radia.vol2d-scalar-replay.v1"
_PHYSICS = {"electrostatic", "current_flow", "steady_heat", "transient_heat"}
_FORMULATIONS = {"planar", "axisymmetric"}
_EPS0 = 8.8541878128e-12


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be finite")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _positive(value: Any, label: str) -> float:
    result = _finite(value, label)
    if result <= 0.0:
        raise ValueError(f"{label} must be positive")
    return result


def _complex(value: Any, label: str) -> complex:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 2:
            raise ValueError(f"{label} complex pair must contain [real, imag]")
        result = complex(_finite(value[0], label), _finite(value[1], label))
    else:
        result = complex(_finite(value, label), 0.0)
    if not (math.isfinite(result.real) and math.isfinite(result.imag)):
        raise ValueError(f"{label} must be finite")
    return result


def _pair(value: complex | float) -> list[float]:
    parsed = complex(value)
    return [float(parsed.real), float(parsed.imag)]


def _real(value: Any, label: str) -> float:
    parsed = complex(value)
    scale = max(1.0, abs(parsed.real))
    if abs(parsed.imag) > 1.0e-10 * scale:
        raise ValueError(f"{label} unexpectedly contains an imaginary part")
    return float(parsed.real)


def _tensor(value: Any, label: str) -> tuple[float, float]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 2:
            raise ValueError(f"{label} diagonal tensor must contain [xx, yy]")
        xx = _positive(value[0], f"{label}[0]")
        yy = _positive(value[1], f"{label}[1]")
        return xx, yy
    scalar = _positive(value, label)
    return scalar, scalar


def _nonnegative_tensor(value: Any, label: str) -> tuple[float, float]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 2:
            raise ValueError(f"{label} diagonal tensor must contain [xx, yy]")
        xx = _finite(value[0], f"{label}[0]")
        yy = _finite(value[1], f"{label}[1]")
    else:
        xx = yy = _finite(value, label)
    if xx < 0.0 or yy < 0.0:
        raise ValueError(f"{label} must be nonnegative")
    return xx, yy


def _material_contract(
    raw: Any,
    material_names: Sequence[str],
    physics: str,
    frequency_hz: float,
) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ValueError("materials must be an object keyed by .vol material name")
    provided = {str(key) for key in raw}
    missing = sorted(set(material_names) - provided)
    unknown = sorted(provided - set(material_names))
    if missing or unknown:
        raise ValueError(f"materials mismatch: missing={missing}, unknown={unknown}")
    units = {
        "electrostatic": ("F/m", "C/m^3"),
        "current_flow": ("S/m", "A/m^3"),
        "steady_heat": ("W/(m K)", "W/m^3"),
        "transient_heat": ("W/(m K)", "W/m^3"),
    }[physics]
    normalized: dict[str, Any] = {}
    for name in material_names:
        row = raw[name]
        if not isinstance(row, Mapping):
            raise ValueError(f"materials[{name}] must be an object")
        allowed = {"coefficient_si", "volumetric_source_si"}
        if physics == "transient_heat":
            allowed.add("volumetric_heat_capacity_j_per_m3_k")
        if physics == "current_flow":
            allowed.add("relative_permittivity")
        extra = sorted({str(key) for key in row} - allowed)
        if extra:
            raise ValueError(f"materials[{name}] contains unsupported or mixed-unit keys: {extra}")
        if "coefficient_si" not in row:
            raise ValueError(f"materials[{name}].coefficient_si is required")
        if physics == "current_flow":
            xx, yy = _nonnegative_tensor(
                row["coefficient_si"], f"materials[{name}].coefficient_si"
            )
        else:
            xx, yy = _tensor(row["coefficient_si"], f"materials[{name}].coefficient_si")
        source = _finite(row.get("volumetric_source_si", 0.0), f"materials[{name}].volumetric_source_si")
        relative_permittivity = (0.0, 0.0)
        if physics == "current_flow":
            relative_permittivity = _nonnegative_tensor(
                row.get("relative_permittivity", 0.0),
                f"materials[{name}].relative_permittivity",
            )
            if frequency_hz == 0.0 and relative_permittivity != (0.0, 0.0):
                raise ValueError("relative_permittivity requires positive frequency_hz")
            if frequency_hz == 0.0 and (xx <= 0.0 or yy <= 0.0):
                raise ValueError("DC current_flow conductivity must be positive in both directions")
            if frequency_hz > 0.0 and (
                (xx <= 0.0 and relative_permittivity[0] <= 0.0)
                or (yy <= 0.0 and relative_permittivity[1] <= 0.0)
            ):
                raise ValueError(
                    "AC current_flow needs conductivity or permittivity in both directions"
                )
        normalized[name] = {
            "coefficient_si": [xx, yy],
            "coefficient_unit": units[0],
            "volumetric_source_si": source,
            "source_unit": units[1],
            "relative_permittivity": list(relative_permittivity),
        }
        if physics == "transient_heat":
            normalized[name]["volumetric_heat_capacity_j_per_m3_k"] = _positive(
                row.get("volumetric_heat_capacity_j_per_m3_k"),
                f"materials[{name}].volumetric_heat_capacity_j_per_m3_k",
            )
    contract = {
        "schema": "radia.vol2d-scalar-materials.v1",
        "physics": physics,
        "frequency_hz": frequency_hz,
        "materials": normalized,
    }
    contract["contract_sha256"] = _sha(contract)
    return contract


def _boundary_contract(
    request: Mapping[str, Any],
    boundary_names: Sequence[str],
    physics: str,
) -> dict[str, Any]:
    raw_dirichlet = request.get("dirichlet_values")
    if not isinstance(raw_dirichlet, Mapping):
        raise ValueError("dirichlet_values must be an object keyed by boundary name")
    unknown = sorted({str(key) for key in raw_dirichlet} - set(boundary_names))
    if unknown:
        raise ValueError(f"unknown Dirichlet boundaries: {unknown}")
    dirichlet = {
        str(name): _complex(value, f"dirichlet_values[{name}]")
        for name, value in raw_dirichlet.items()
    }
    if physics != "current_flow" and any(abs(value.imag) > 0.0 for value in dirichlet.values()):
        raise ValueError("only current_flow accepts complex Dirichlet values")

    raw_robin = request.get("robin_boundaries", {})
    if not isinstance(raw_robin, Mapping):
        raise ValueError("robin_boundaries must be an object")
    if physics not in {"steady_heat", "transient_heat"} and raw_robin:
        raise ValueError("Robin convection is available only for heat studies")
    robin: dict[str, dict[str, float]] = {}
    for name, raw in raw_robin.items():
        key = str(name)
        if key not in boundary_names:
            raise ValueError(f"unknown Robin boundary: {key}")
        if key in dirichlet:
            raise ValueError(f"boundary {key} cannot be both Dirichlet and Robin")
        if not isinstance(raw, Mapping):
            raise ValueError(f"robin_boundaries[{key}] must be an object")
        allowed = {"transfer_w_per_m2_k", "ambient_k"}
        extra = sorted({str(item) for item in raw} - allowed)
        if extra:
            raise ValueError(f"robin_boundaries[{key}] contains unsupported keys: {extra}")
        robin[key] = {
            "transfer_w_per_m2_k": _positive(
                raw.get("transfer_w_per_m2_k"),
                f"robin_boundaries[{key}].transfer_w_per_m2_k",
            ),
            "ambient_k": _finite(raw.get("ambient_k"), f"robin_boundaries[{key}].ambient_k"),
        }
    if physics != "transient_heat" and not dirichlet and not robin:
        raise ValueError("unconstrained scalar nullspace: specify Dirichlet or positive Robin data")
    unit = "K" if physics in {"steady_heat", "transient_heat"} else "V"
    contract = {
        "schema": "radia.vol2d-scalar-boundaries.v1",
        "dirichlet_values": {name: _pair(value) for name, value in sorted(dirichlet.items())},
        "dirichlet_unit": unit,
        "robin_boundaries": dict(sorted(robin.items())),
        "natural_boundaries": sorted(set(boundary_names) - set(dirichlet) - set(robin)),
    }
    contract["contract_sha256"] = _sha(contract)
    return contract


def _dense_matrix(matrix: Any, *, complex_matrix: bool) -> np.ndarray:
    rows, columns, values = matrix.COO()
    dtype = complex if complex_matrix else float
    dense = np.zeros((matrix.height, matrix.width), dtype=dtype)
    np.add.at(
        dense,
        (np.asarray(rows, dtype=int), np.asarray(columns, dtype=int)),
        np.asarray(values, dtype=dtype),
    )
    return dense


def _dof_indices(bits: Any) -> np.ndarray:
    return np.asarray([index for index, flag in enumerate(bits) if flag], dtype=int)


def _gmsh_export(
    mesh: Any,
    gfu: Any,
    field: Any,
    *,
    basename: str,
    request_sha256: str,
    complex_field: bool,
) -> dict[str, str]:
    from ngsolve import CoefficientFunction, Conj, InnerProduct, sqrt  # type: ignore
    from radia.gmsh_post_export import GmshPostExport

    root = Path(r"C:\temp") / "radia_mcp_vol2d_scalar" / request_sha256[:20]
    root.mkdir(parents=True, exist_ok=True)
    msh = root / f"{basename}.msh"
    post = GmshPostExport(mesh)
    if complex_field:
        post.add_scalar_field("primary_real", gfu.real)
        post.add_scalar_field("primary_magnitude", sqrt((gfu * Conj(gfu)).real))
        field_real = field.real
    else:
        post.add_scalar_field("primary", gfu)
        field_real = field
    post.add_scalar_field("field_magnitude", sqrt(InnerProduct(field_real, field_real)))
    post.add_vector_field(
        "field_vector",
        CoefficientFunction((field_real[0], field_real[1], 0.0)),
    )
    post.write(str(msh))
    geo = msh.with_suffix(".geo")
    paths = {
        "gmsh_msh": msh,
        "gmsh_geo": geo,
        "gmsh_geo_opt": Path(str(geo) + ".opt"),
        "gmsh_msh_opt": Path(str(msh) + ".opt"),
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"Gmsh exporter did not create required companions: {missing}")
    return {name: path.read_text(encoding="utf-8") for name, path in paths.items()}


def _csv_export(observables: Mapping[str, Any]) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["observable", "value_json"])
    for name, value in sorted(observables.items()):
        writer.writerow([name, json.dumps(value, sort_keys=True, separators=(",", ":"))])
    return output.getvalue()


def _electrostatic_boundary_force(
    mesh_view: Any,
    mesh: Any,
    electric_field: Any,
    materials: Mapping[str, Any],
    boundary_name: str,
    *,
    formulation: str,
    model_depth_m: float | None,
) -> list[float]:
    """Integrate dielectric-side Maxwell traction on one boundary.

    NGSolve's boundary trace of ``grad(H1)`` is tangential and therefore loses
    the conductor-normal electric field.  Sample from the adjacent volume
    element instead, with a three-point Gauss rule on every boundary edge.
    """

    boundary_numbers = {
        number
        for number, name in mesh_view.boundary_names.items()
        if name == boundary_name
    }
    gauss = (
        (0.5 - math.sqrt(15.0) / 10.0, 5.0 / 18.0),
        (0.5, 8.0 / 18.0),
        (0.5 + math.sqrt(15.0) / 10.0, 5.0 / 18.0),
    )
    force = np.zeros(2, dtype=float)
    matched = 0
    for edge in mesh_view.boundary_edges:
        if edge.boundary_number not in boundary_numbers:
            continue
        edge_nodes = set(edge.nodes)
        adjacent_cells = [
            cell for cell in mesh_view.cells if edge_nodes.issubset(cell.nodes)
        ]
        if len(adjacent_cells) != 1:
            raise ValueError(
                f"force boundary {boundary_name} edge must have one adjacent material cell"
            )
        cell = adjacent_cells[0]
        material = materials["materials"][mesh_view.material_name(cell.material_number)]
        eps_x, eps_y = material["coefficient_si"]
        p0 = np.asarray(mesh_view.points[edge.nodes[0] - 1][:2], dtype=float)
        p1 = np.asarray(mesh_view.points[edge.nodes[1] - 1][:2], dtype=float)
        tangent = p1 - p0
        length = float(np.linalg.norm(tangent))
        if length <= 0.0:
            raise ValueError(f"force boundary {boundary_name} contains a zero-length edge")
        midpoint = 0.5 * (p0 + p1)
        centroid = np.mean(
            [np.asarray(mesh_view.points[node - 1][:2], dtype=float) for node in cell.nodes],
            axis=0,
        )
        normal = np.asarray((tangent[1], -tangent[0]), dtype=float) / length
        if float(np.dot(normal, midpoint - centroid)) < 0.0:
            normal = -normal
        inward = centroid - midpoint
        inward_norm = float(np.linalg.norm(inward))
        if inward_norm <= 0.0:
            raise ValueError(f"force boundary {boundary_name} has an invalid adjacent cell")
        inward /= inward_norm
        offset = 1.0e-8 * max(length, math.sqrt(mesh_view.cell_area(cell)))
        for coordinate, weight in gauss:
            point = (1.0 - coordinate) * p0 + coordinate * p1 + offset * inward
            value = electric_field(mesh(float(point[0]), float(point[1])))
            electric = np.asarray((float(value[0]), float(value[1])), dtype=float)
            displacement = np.asarray(
                (eps_x * electric[0], eps_y * electric[1]), dtype=float
            )
            traction = (
                displacement * float(np.dot(electric, normal))
                - 0.5 * float(np.dot(electric, displacement)) * normal
            )
            measure = length * weight
            if formulation == "planar":
                measure *= float(model_depth_m)
            else:
                measure *= 2.0 * math.pi * float(point[0])
            force += traction * measure
        matched += 1
    if matched == 0:
        raise ValueError(f"force boundary {boundary_name} has no boundary edges")
    if formulation == "axisymmetric":
        return [0.0, float(force[1])]
    return [float(force[0]), float(force[1])]


def _prepare(request: Mapping[str, Any]) -> tuple[dict[str, Any], Any, dict[str, Any], dict[str, Any]]:
    physics = str(request.get("physics", ""))
    if physics not in _PHYSICS:
        raise ValueError(f"physics must be one of {sorted(_PHYSICS)}")
    formulation = str(request.get("formulation", ""))
    if formulation not in _FORMULATIONS:
        raise ValueError(f"formulation must be one of {sorted(_FORMULATIONS)}")
    text = request.get("vol_text")
    if not isinstance(text, str) or not text.strip():
        raise ValueError("vol_text must be non-empty")
    mesh_view = parse_netgen_2d_vol(text, source_name=str(request.get("source_name", "input.vol")))
    family = str(request.get("element_family", "P1"))
    family_contract = _family_contract(mesh_view, family, "planar")
    if formulation == "axisymmetric":
        radii = [point[0] for point in mesh_view.points]
        scale = max(1.0, *(abs(value) for value in radii))
        if min(radii) < -1.0e-12 * scale or max(radii) <= 0.0:
            raise ValueError("axisymmetric .vol requires a nonnegative radius coordinate and r > 0 extent")
        if request.get("model_depth_m") is not None:
            raise ValueError("axisymmetric full-revolution analysis must not specify model_depth_m")
        model_depth_m = None
    else:
        model_depth_m = _positive(request.get("model_depth_m"), "model_depth_m")
    frequency_hz = _finite(request.get("frequency_hz", 0.0), "frequency_hz")
    if frequency_hz < 0.0:
        raise ValueError("frequency_hz must be nonnegative")
    if physics != "current_flow" and frequency_hz != 0.0:
        raise ValueError("frequency_hz is available only for current_flow")
    material_contract = _material_contract(
        request.get("materials"),
        list(mesh_view.contract()["material_names"]),
        physics,
        frequency_hz,
    )
    boundary_contract = _boundary_contract(
        request,
        list(mesh_view.contract()["boundary_names"]),
        physics,
    )
    prepared = {
        "physics": physics,
        "formulation": formulation,
        "element_family": family,
        "family_contract": family_contract,
        "model_depth_m": model_depth_m,
        "frequency_hz": frequency_hz,
    }
    return prepared, mesh_view, material_contract, boundary_contract


def solve_vol2d_scalar(request: Mapping[str, Any]) -> dict[str, Any]:
    """Solve one scalar field and return a replayable result artifact."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    prepare_started = time.perf_counter()
    prepared, mesh_view, materials, boundaries = _prepare(request)
    if prepared["physics"] == "transient_heat":
        raise ValueError("transient_heat requires operation='transient_heat'")
    basename = str(request.get("export_basename", "vol2d_scalar")).strip()
    if not _SAFE_NAME.fullmatch(basename):
        raise ValueError("export_basename must be a portable filename stem")
    mesh_contract = mesh_view.contract()
    request_contract = {
        "schema": "radia.vol2d-scalar-request.v1",
        **prepared,
        "mesh_contract": mesh_contract,
        "mesh_contract_sha256": mesh_contract["contract_sha256"],
        "material_contract": materials,
        "material_contract_sha256": materials["contract_sha256"],
        "boundary_contract": boundaries,
        "boundary_contract_sha256": boundaries["contract_sha256"],
        "terminal_pair": request.get("terminal_pair"),
        "force_boundaries": request.get("force_boundaries", []),
        "export_basename": basename,
    }
    request_sha = _sha(request_contract)
    for key in ("mesh_contract_sha256", "material_contract_sha256", "boundary_contract_sha256"):
        expected = request.get(f"expected_{key}")
        if expected is not None and str(expected) != str(request_contract[key]):
            raise ValueError(f"expected_{key} does not match the prepared request")

    from ngsolve import (  # type: ignore
        BND,
        H1,
        BilinearForm,
        GridFunction,
        Conj,
        Integrate,
        LinearForm,
        Mesh,
        grad,
        ds,
        dx,
        x,
    )

    runtime = _runtime_vol_path(str(request["vol_text"]), mesh_view.content_sha256)
    mesh = Mesh(str(runtime))
    if prepared["family_contract"]["curved_geometry"]:
        mesh.Curve(prepared["family_contract"]["order"])
    physics = prepared["physics"]
    complex_field = physics == "current_flow" and prepared["frequency_hz"] > 0.0
    dirichlet_pairs = boundaries["dirichlet_values"]
    dirichlet_names = list(dirichlet_pairs)
    fes = H1(
        mesh,
        order=prepared["family_contract"]["order"],
        complex=complex_field,
        dirichlet="|".join(dirichlet_names),
    )
    trial, test = fes.TnT()
    omega = 2.0 * math.pi * prepared["frequency_hz"]
    coeff_x: dict[str, complex | float] = {}
    coeff_y: dict[str, complex | float] = {}
    source_values: dict[str, float] = {}
    for name, row in materials["materials"].items():
        xx, yy = row["coefficient_si"]
        if complex_field:
            erx, ery = row["relative_permittivity"]
            coeff_x[name] = complex(xx, omega * _EPS0 * erx)
            coeff_y[name] = complex(yy, omega * _EPS0 * ery)
        else:
            coeff_x[name] = xx
            coeff_y[name] = yy
        source_values[name] = row["volumetric_source_si"]
    cxx = mesh.MaterialCF(coeff_x)
    cyy = mesh.MaterialCF(coeff_y)
    source = mesh.MaterialCF(source_values)
    if prepared["formulation"] == "planar":
        measure = prepared["model_depth_m"]
        measure_contract = {"kind": "planar_depth", "model_depth_m": measure}
    else:
        measure = 2.0 * math.pi * x
        measure_contract = {"kind": "axisymmetric_full_revolution", "weight": "2*pi*r"}

    form = BilinearForm(fes, symmetric=True)
    form += (cxx * grad(trial)[0] * grad(test)[0] + cyy * grad(trial)[1] * grad(test)[1]) * measure * dx
    rhs = LinearForm(fes)
    rhs += source * test * measure * dx
    for name, row in boundaries["robin_boundaries"].items():
        h = row["transfer_w_per_m2_k"]
        ambient = row["ambient_k"]
        region = mesh.Boundaries(name)
        form += h * trial * test * measure * ds(definedon=region)
        rhs += h * ambient * test * measure * ds(definedon=region)
    form.Assemble()
    rhs.Assemble()
    gfu = GridFunction(fes)
    if dirichlet_names:
        values = {name: complex(*pair) if complex_field else pair[0] for name, pair in dirichlet_pairs.items()}
        gfu.Set(mesh.BoundaryCF(values, default=0.0), BND)
    matrix = _dense_matrix(form.mat, complex_matrix=complex_field)
    load = np.asarray(rhs.vec.FV().NumPy(), dtype=complex if complex_field else float).copy()
    initial = np.asarray(gfu.vec.FV().NumPy(), dtype=complex if complex_field else float).copy()
    free = _dof_indices(fes.FreeDofs())
    if len(free) == 0:
        raise ValueError("scalar solve has no free degrees of freedom")
    reduced = matrix[np.ix_(free, free)]
    reduced_rhs = (load - matrix @ initial)[free]
    operator_contract = {
        "matrix_shape": list(matrix.shape),
        "free_dofs": free.tolist(),
        "measure": measure_contract,
        "complex": complex_field,
        "matrix_sha256": _sha([[_pair(value) for value in row] for row in matrix]),
    }
    operator_sha = _sha(operator_contract)
    if request.get("expected_operator_sha256") is not None:
        if str(request["expected_operator_sha256"]) != operator_sha:
            raise ValueError("expected_operator_sha256 does not match the assembled operator")
    prepare_assembly_s = time.perf_counter() - prepare_started

    solve_started = time.perf_counter()
    try:
        if complex_field:
            q, r = np.linalg.qr(reduced)
            solved = np.linalg.solve(r, q.conj().T @ reduced_rhs)
            factor_payload = {"kind": "qr", "q": [[_pair(v) for v in row] for row in q], "r": [[_pair(v) for v in row] for row in r]}
        else:
            factor = np.linalg.cholesky(0.5 * (reduced + reduced.T))
            solved = np.linalg.solve(factor.T, np.linalg.solve(factor, reduced_rhs))
            factor_payload = {"kind": "cholesky", "factor": factor.tolist()}
    except np.linalg.LinAlgError as exc:
        raise ValueError("scalar operator is singular or not factorable") from exc
    state = initial.copy()
    state[free] += solved
    gfu.vec.FV().NumPy()[:] = state
    residual = matrix @ state - load
    factorization_contract = {
        "kind": factor_payload["kind"],
        "operator_sha256": operator_sha,
        "factors_sha256": _sha(factor_payload),
    }
    factorization_sha = _sha(factorization_contract)
    if request.get("expected_factorization_sha256") is not None:
        if str(request["expected_factorization_sha256"]) != factorization_sha:
            raise ValueError("expected_factorization_sha256 does not match the solved operator")
    factor_solve_s = time.perf_counter() - solve_started

    observable_started = time.perf_counter()
    field = -grad(gfu)
    boundary_reactions: dict[str, list[float]] = {}
    boundary_reaction_values: dict[str, complex] = {}
    for name in dirichlet_names:
        indicator = GridFunction(fes)
        indicator.Set(1.0, definedon=mesh.Boundaries(name))
        indicator_vector = np.asarray(
            indicator.vec.FV().NumPy(), dtype=complex if complex_field else float
        )
        reaction = complex(np.dot(indicator_vector, residual))
        boundary_reaction_values[name] = reaction
        boundary_reactions[name] = _pair(reaction)

    terminal_pair = request.get("terminal_pair")
    terminal: dict[str, Any] | None = None
    if terminal_pair is not None:
        if not isinstance(terminal_pair, Mapping):
            raise ValueError("terminal_pair must be an object")
        positive = str(terminal_pair.get("positive_boundary", ""))
        negative = str(terminal_pair.get("negative_boundary", ""))
        if positive == negative or positive not in dirichlet_pairs or negative not in dirichlet_pairs:
            raise ValueError("terminal_pair must name two distinct Dirichlet boundaries")
        response_positive = boundary_reaction_values[positive]
        response_negative = boundary_reaction_values[negative]
        voltage = complex(*dirichlet_pairs[positive]) - complex(*dirichlet_pairs[negative])
        if abs(voltage) == 0.0:
            raise ValueError("terminal_pair requires a nonzero value difference")
        terminal = {
            "positive_boundary": positive,
            "negative_boundary": negative,
            "value_difference": _pair(voltage),
            "positive_reaction": _pair(response_positive),
            "negative_reaction": _pair(response_negative),
            "reaction_closure": _pair(response_positive + response_negative),
        }

    if complex_field:
        operator_quadratic = 0.5 * float(np.real(np.vdot(state, matrix @ state)))
        gradient_quadratic = 0.5 * _real(
            Integrate(
                (
                    mesh.MaterialCF({name: row["coefficient_si"][0] for name, row in materials["materials"].items()})
                    * grad(gfu)[0] * Conj(grad(gfu)[0])
                    + mesh.MaterialCF({name: row["coefficient_si"][1] for name, row in materials["materials"].items()})
                    * grad(gfu)[1] * Conj(grad(gfu)[1])
                ).real * measure,
                mesh,
            ),
            "gradient quadratic",
        )
    else:
        operator_quadratic = 0.5 * float(state @ matrix @ state)
        gradient_quadratic = 0.5 * _real(
            Integrate(
                (cxx * grad(gfu)[0] * grad(gfu)[0] + cyy * grad(gfu)[1] * grad(gfu)[1])
                * measure,
                mesh,
            ),
            "gradient quadratic",
        )
    source_total = _real(Integrate(source * measure, mesh), "source integral")
    observables: dict[str, Any] = {
        "operator_quadratic_half": operator_quadratic,
        "field_gradient_quadratic_half": gradient_quadratic,
        "volumetric_source_total": source_total,
        "residual_inf": float(np.linalg.norm(residual[free], ord=np.inf)),
        "boundary_reactions": boundary_reactions,
    }
    if physics == "electrostatic":
        observables["electric_energy_j"] = gradient_quadratic
        if terminal:
            dv = abs(complex(*terminal["value_difference"]))
            observables["capacitance_f"] = 2.0 * gradient_quadratic / (dv * dv)
            observables["terminal_charge_c"] = terminal["positive_reaction"]
        force_names = request.get("force_boundaries", [])
        if isinstance(force_names, str):
            force_names = [force_names]
        if not isinstance(force_names, Sequence) or isinstance(force_names, (bytes, str)):
            raise ValueError("force_boundaries must be a sequence of boundary names")
        unknown_force = sorted({str(name) for name in force_names} - set(mesh_view.contract()["boundary_names"]))
        if unknown_force:
            raise ValueError(f"unknown force boundaries: {unknown_force}")
        if force_names:
            electric = field
            force_rows: dict[str, list[float]] = {}
            for name in force_names:
                force_rows[str(name)] = _electrostatic_boundary_force(
                    mesh_view,
                    mesh,
                    electric,
                    materials,
                    str(name),
                    formulation=prepared["formulation"],
                    model_depth_m=prepared["model_depth_m"],
                )
            observables["electrostatic_force_by_boundary_n"] = force_rows
            observables["electrostatic_force_convention"] = (
                "integral of dielectric Maxwell traction T*n_domain; negate for force on an excluded conductor"
            )
    elif physics == "current_flow":
        if terminal:
            voltage = complex(*terminal["value_difference"])
            current = complex(*terminal["positive_reaction"])
            admittance = current / voltage
            observables["terminal_current_a"] = _pair(current)
            observables["admittance_s"] = _pair(admittance)
            observables["complex_power_va_rms"] = _pair(voltage * current.conjugate())
            if prepared["frequency_hz"] == 0.0:
                power = float((voltage * current.conjugate()).real)
                observables["conduction_power_w"] = power
                observables["resistance_ohm"] = float(abs(voltage) ** 2 / power)
    else:
        convection_total = 0.0
        convection_rows: dict[str, float] = {}
        for name, row in boundaries["robin_boundaries"].items():
            flux = _real(
                Integrate(
                    row["transfer_w_per_m2_k"] * (gfu - row["ambient_k"])
                    * measure * ds(definedon=mesh.Boundaries(name)),
                    mesh,
                ),
                f"convection[{name}]",
            )
            convection_rows[name] = flux
            convection_total += flux
        fixed_outflow = 0.0
        if dirichlet_names:
            fixed_dofs = sorted(set().union(*(
                set(_dof_indices(fes.GetDofs(mesh.Boundaries(name))).tolist())
                for name in dirichlet_names
            )))
            fixed_outflow = -float(np.sum(np.real(residual[fixed_dofs])))
        observables.update({
            "generated_heat_w": source_total,
            "convection_outflow_w": convection_total,
            "convection_by_boundary_w": convection_rows,
            "dirichlet_outflow_w": fixed_outflow,
            "heat_balance_residual_w": source_total - convection_total - fixed_outflow,
        })
    observables_sha = _sha(observables)
    if request.get("expected_observables_sha256") is not None:
        if str(request["expected_observables_sha256"]) != observables_sha:
            raise ValueError("expected_observables_sha256 does not match recomputed observables")
    observable_s = time.perf_counter() - observable_started

    export_started = time.perf_counter()
    gmsh = _gmsh_export(
        mesh,
        gfu,
        field,
        basename=basename,
        request_sha256=request_sha,
        complex_field=complex_field,
    )
    csv_content = _csv_export(observables)
    state_rows = [_pair(value) for value in state]
    state_sha = _sha(state_rows)
    if request.get("expected_field_state_sha256") is not None:
        if str(request["expected_field_state_sha256"]) != state_sha:
            raise ValueError("expected_field_state_sha256 does not match the solved field")
    result_contract = {
        "schema": SCALAR_SCHEMA,
        "status": "solved",
        "operation": "solve",
        "request_contract": request_contract,
        "request_contract_sha256": request_sha,
        "mesh_contract_sha256": request_contract["mesh_contract_sha256"],
        "material_contract_sha256": request_contract["material_contract_sha256"],
        "boundary_contract_sha256": request_contract["boundary_contract_sha256"],
        "operator_contract": operator_contract,
        "operator_sha256": operator_sha,
        "factorization_contract": factorization_contract,
        "factorization_sha256": factorization_sha,
        "field_state": state_rows,
        "field_state_sha256": state_sha,
        "terminal": terminal,
        "observables": observables,
        "observables_sha256": observables_sha,
        "factorization_count": 1,
        "solve_count": 1,
        "generated_vol_git_required": False,
    }
    json_content = _canonical(result_contract)
    exports = {
        "json": _export_entry(json_content, f"{basename}.json", "application/json"),
        "csv": _export_entry(csv_content, f"{basename}.csv", "text/csv"),
        "gmsh_msh": _export_entry(gmsh["gmsh_msh"], f"{basename}.msh", "model/mesh"),
        "gmsh_geo": _export_entry(gmsh["gmsh_geo"], f"{basename}.geo", "text/plain"),
        "gmsh_geo_opt": _export_entry(gmsh["gmsh_geo_opt"], f"{basename}.geo.opt", "text/plain"),
        "gmsh_msh_opt": _export_entry(gmsh["gmsh_msh_opt"], f"{basename}.msh.opt", "text/plain"),
    }
    result_contract["export_content_sha256"] = {name: row["sha256"] for name, row in exports.items()}
    expected_exports = request.get("expected_export_content_sha256")
    if expected_exports is not None:
        if not isinstance(expected_exports, Mapping) or dict(expected_exports) != result_contract["export_content_sha256"]:
            raise ValueError("expected_export_content_sha256 does not match generated exports")
    result_contract["canonical_json_sha256"] = exports["json"]["sha256"]
    export_s = time.perf_counter() - export_started
    total = prepare_assembly_s + factor_solve_s + observable_s + export_s
    return {
        "schema": SCALAR_SCHEMA,
        "status": "solved",
        "operation": "solve",
        "execution_version": {
            "radia_mcp": _package_version(),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "timing_s": {
            "prepare_and_assemble": prepare_assembly_s,
            "factorize_and_solve": factor_solve_s,
            "observables": observable_s,
            "export": export_s,
            "total": total,
        },
        "result_contract": result_contract,
        "exports": exports,
    }


def scalar_replay_gate(artifact: Mapping[str, Any]) -> dict[str, Any]:
    """Reject stale or internally inconsistent scalar artifacts."""

    if not isinstance(artifact, Mapping):
        raise ValueError("replay_artifact must be an object")
    contract = artifact.get("result_contract")
    exports = artifact.get("exports")
    if not isinstance(contract, Mapping) or not isinstance(exports, Mapping):
        raise ValueError("replay artifact needs result_contract and exports")
    checks = {
        "schema": contract.get("schema") == SCALAR_SCHEMA,
        "status": contract.get("status") == "solved",
        "request_contract_sha256": _sha(contract.get("request_contract")) == contract.get("request_contract_sha256"),
        "mesh_contract_sha256": contract.get("mesh_contract_sha256") == contract.get("request_contract", {}).get("mesh_contract_sha256"),
        "material_contract_sha256": contract.get("material_contract_sha256") == contract.get("request_contract", {}).get("material_contract_sha256"),
        "boundary_contract_sha256": contract.get("boundary_contract_sha256") == contract.get("request_contract", {}).get("boundary_contract_sha256"),
        "operator_sha256": _sha(contract.get("operator_contract")) == contract.get("operator_sha256"),
        "factorization_sha256": _sha(contract.get("factorization_contract")) == contract.get("factorization_sha256"),
        "factorization_operator": contract.get("factorization_contract", {}).get("operator_sha256") == contract.get("operator_sha256"),
        "field_state_sha256": _sha(contract.get("field_state")) == contract.get("field_state_sha256"),
        "observables_sha256": _sha(contract.get("observables")) == contract.get("observables_sha256"),
        "single_factorization": contract.get("factorization_count") == 1,
        "single_solve": contract.get("solve_count") == 1,
    }
    hashes = contract.get("export_content_sha256", {})
    for name in ("csv", "gmsh_msh", "gmsh_geo", "gmsh_geo_opt", "gmsh_msh_opt"):
        row = exports.get(name)
        checks[f"{name}_sha256"] = (
            isinstance(row, Mapping)
            and _sha(str(row.get("content", ""))) == row.get("sha256")
            and row.get("sha256") == hashes.get(name)
        )
    json_row = exports.get("json")
    stored = dict(contract)
    canonical_sha = stored.pop("canonical_json_sha256", None)
    stored.pop("export_content_sha256", None)
    expected_json = _canonical(stored)
    checks["canonical_json"] = (
        isinstance(json_row, Mapping)
        and json_row.get("content") == expected_json
        and _sha(expected_json) == json_row.get("sha256")
        and canonical_sha == json_row.get("sha256")
    )
    checks["gmsh_v41"] = "$MeshFormat\n4.1 " in str(exports.get("gmsh_msh", {}).get("content", ""))
    passed = all(checks.values())
    return {
        "schema": REPLAY_SCHEMA,
        "status": "accepted" if passed else "rejected",
        "pass": passed,
        "checks": checks,
        "request_contract_sha256": contract.get("request_contract_sha256"),
    }


def analyze_vol2d_scalar(request: Mapping[str, Any]) -> dict[str, Any]:
    """Dispatch the closed-world solve or replay operation."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    operation = str(request.get("operation", "solve"))
    if operation == "solve":
        return solve_vol2d_scalar(request)
    if operation == "transient_heat":
        from .vol2d_thermal import solve_vol2d_transient_heat

        return solve_vol2d_transient_heat(request)
    if operation == "electrostatic_system":
        from .vol2d_electrostatic import solve_vol2d_electrostatic_system

        return solve_vol2d_electrostatic_system(request)
    if operation == "replay_gate":
        return scalar_replay_gate(request.get("replay_artifact"))
    raise ValueError(
        "operation must be solve, transient_heat, electrostatic_system, or replay_gate"
    )
