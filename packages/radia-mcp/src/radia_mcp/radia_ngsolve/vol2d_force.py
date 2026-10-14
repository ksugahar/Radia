"""Boundary-aware force extraction from solved dimension-2 Netgen ``.vol`` models.

The readable contract is deliberately small: solve the same field operators as
``vol2d_dynamics``, reconstruct the finite-element field, and extract force in
an air-only weighted-stress band.  A conductor target can additionally use a
Lorentz integral.  A passive magnetic target cannot: its force must come from
weighted stress and an independent virtual-work sweep.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from typing import Any

import numpy as np

from .force import eggshell_force_2d, eggshell_force_axi, lorentz_force_2d
from .vol2d_dynamics import (
    _open_space,
    _prepare_request,
    assemble_vol2d_dynamics,
    solve_vol2d_nonlinear_static,
)


FORCE_SCHEMA = "radia.vol2d-force-analysis.v1"
FORCE_TARGET_SCHEMA = "radia.vol2d-force-target.v1"
FORCE_REFINEMENT_SCHEMA = "radia.vol2d-force-refinement.v1"
VIRTUAL_WORK_SCHEMA = "radia.vol2d-force-virtual-work.v1"


def _package_version() -> str:
    try:
        return version("radia-mcp")
    except PackageNotFoundError:
        return "source-tree"


def _sha(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _finite(value: Any, label: str) -> float:
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


def _vector2(value: Any, label: str) -> tuple[float, float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 2:
        raise ValueError(f"{label} must have exactly two components")
    return _finite(value[0], f"{label}[0]"), _finite(value[1], f"{label}[1]")


def _boundary_contract(mesh_view: Any) -> dict[str, Any]:
    counts: dict[str, int] = {}
    rows: list[dict[str, Any]] = []
    for edge in mesh_view.boundary_edges:
        name = mesh_view.boundary_names.get(
            edge.boundary_number, f"boundary_{edge.boundary_number}"
        )
        counts[name] = counts.get(name, 0) + 1
        rows.append({"name": name, "nodes_1based": list(edge.nodes)})
    payload = {
        "schema": "radia.vol2d-boundary-contract.v1",
        "boundary_edge_counts": dict(sorted(counts.items())),
        "edges": rows,
    }
    payload["contract_sha256"] = _sha(payload)
    return payload


def _cell_radius_range(mesh_view: Any, cell: Any, center: tuple[float, float]) -> tuple[float, float]:
    radii = [
        math.hypot(
            mesh_view.points[node - 1][0] - center[0],
            mesh_view.points[node - 1][1] - center[1],
        )
        for node in cell.nodes
    ]
    return min(radii), max(radii)


def normalize_vol2d_force_target(
    mesh_view: Any,
    material_contract: Mapping[str, Any],
    operator_sha256: str,
    raw: Any,
    *,
    formulation: str,
    dirichlet_boundaries: Sequence[str],
) -> dict[str, Any]:
    """Validate force geometry, units, and lineage before field extraction."""

    if not isinstance(raw, Mapping):
        raise ValueError("force_target must be an object")
    mesh_contract = mesh_view.contract()
    target = str(raw.get("target_material", "")).strip()
    air = str(raw.get("air_material", "")).strip()
    materials = set(mesh_contract["material_names"])
    if target not in materials:
        raise ValueError(f"target_material is absent from .vol: {target or '<empty>'}")
    if air not in materials:
        raise ValueError(f"air_material is absent from .vol: {air or '<empty>'}")
    if target == air:
        raise ValueError("target_material and air_material must differ")

    center = _vector2(raw.get("center_m"), "center_m")
    inner = _positive(raw.get("inner_radius_m"), "inner_radius_m")
    outer = _positive(raw.get("outer_radius_m"), "outer_radius_m")
    if not outer > inner:
        raise ValueError("outer_radius_m must exceed inner_radius_m")

    target_cells = [
        cell for cell in mesh_view.cells if mesh_view.material_name(cell.material_number) == target
    ]
    if not target_cells:
        raise ValueError("target_material has no cells")
    target_max = max(_cell_radius_range(mesh_view, cell, center)[1] for cell in target_cells)
    if target_max >= inner * (1.0 - 1.0e-12):
        raise ValueError("inner_radius_m must strictly enclose every target cell")

    band_air_cells = 0
    crossing: set[str] = set()
    for cell in mesh_view.cells:
        name = mesh_view.material_name(cell.material_number)
        radius_min, radius_max = _cell_radius_range(mesh_view, cell, center)
        overlaps = radius_max > inner and radius_min < outer
        if overlaps and name == air:
            band_air_cells += 1
        elif overlaps and name != target:
            crossing.add(name)
    if crossing:
        raise ValueError(f"weighted-stress band crosses non-air materials: {sorted(crossing)}")
    if not band_air_cells:
        raise ValueError("weighted-stress band contains no air cells")

    boundary_nodes = {
        node for edge in mesh_view.boundary_edges for node in edge.nodes
    }
    minimum_boundary_radius = min(
        math.hypot(
            mesh_view.points[node - 1][0] - center[0],
            mesh_view.points[node - 1][1] - center[1],
        )
        for node in boundary_nodes
    )
    if minimum_boundary_radius <= outer:
        raise ValueError("weighted-stress band reaches the exterior boundary")

    outer_names_raw = raw.get("outer_boundary_names")
    if (
        not isinstance(outer_names_raw, Sequence)
        or isinstance(outer_names_raw, (str, bytes))
        or not outer_names_raw
    ):
        raise ValueError("outer_boundary_names must be a non-empty sequence")
    outer_names = sorted({str(value) for value in outer_names_raw})
    available_boundaries = set(mesh_contract["boundary_names"])
    missing = sorted(set(outer_names) - available_boundaries)
    if missing:
        raise ValueError(f"outer boundary names are absent from .vol: {missing}")
    if set(outer_names) != {str(value) for value in dirichlet_boundaries}:
        raise ValueError("outer_boundary_names must match the solved Dirichlet boundary set")

    method = str(raw.get("method", "")).strip()
    allowed_methods = {"weighted_stress", "dual_lorentz_weighted_stress"}
    if method not in allowed_methods:
        raise ValueError(f"force method must be one of {sorted(allowed_methods)}")
    if formulation == "axisymmetric_henrotte" and method != "weighted_stress":
        raise ValueError("axisymmetric force uses full-revolution weighted stress only")

    branch = raw.get("target_branch")
    if branch is not None:
        branch = str(branch).strip()
        if not branch:
            raise ValueError("target_branch must be non-empty when provided")
    if method == "dual_lorentz_weighted_stress" and not branch:
        raise ValueError("dual force extraction requires target_branch")

    agreement = _positive(raw.get("agreement_relative_tolerance", 0.2), "agreement_relative_tolerance")
    if agreement > 1.0:
        raise ValueError("agreement_relative_tolerance must not exceed one")

    expected = {
        "mesh_contract_sha256": mesh_contract["contract_sha256"],
        "material_contract_sha256": material_contract["contract_sha256"],
        "operator_sha256": operator_sha256,
    }
    for key, actual in expected.items():
        if str(raw.get(key, "")) != actual:
            raise ValueError(f"{key} does not match the solved .vol operators")

    if formulation == "planar":
        depth = _positive(raw.get("model_depth_m"), "model_depth_m")
        frame = "global_cartesian_xy"
        unit_basis = "N_per_m_out_of_plane_and_depth_integrated_N"
    elif formulation == "axisymmetric_henrotte":
        if raw.get("model_depth_m") is not None:
            raise ValueError("axisymmetric full-revolution force must not specify model_depth_m")
        depth = None
        frame = "meridional_rz_axial_resultant_only"
        unit_basis = "full_revolution_N"
    else:
        raise ValueError("unsupported force formulation")

    declared_frame = str(raw.get("force_component_frame", frame))
    if declared_frame != frame:
        raise ValueError(f"force_component_frame must be {frame}")

    boundary = _boundary_contract(mesh_view)
    normalized = {
        "schema": FORCE_TARGET_SCHEMA,
        "target_material": target,
        "air_material": air,
        "center_m": list(center),
        "inner_radius_m": inner,
        "outer_radius_m": outer,
        "outer_boundary_names": outer_names,
        "boundary_contract_sha256": boundary["contract_sha256"],
        "method": method,
        "target_branch": branch,
        "agreement_relative_tolerance": agreement,
        "formulation": formulation,
        "force_component_frame": frame,
        "force_unit_basis": unit_basis,
        "model_depth_m": depth,
        "mesh_contract_sha256": expected["mesh_contract_sha256"],
        "material_contract_sha256": expected["material_contract_sha256"],
        "operator_sha256": expected["operator_sha256"],
        "band_air_cell_count": band_air_cells,
        "target_max_radius_m": target_max,
        "minimum_exterior_boundary_radius_m": minimum_boundary_radius,
    }
    normalized["force_request_sha256"] = _sha(normalized)
    return {"target": normalized, "boundary_contract": boundary}


def _field_state(request: Mapping[str, Any], operators: Mapping[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    currents = np.asarray(request.get("branch_current_a"), dtype=float)
    source = np.asarray(operators["assembly"]["source_matrix"], dtype=float)
    if currents.ndim != 1 or currents.size != source.shape[1] or not np.all(np.isfinite(currents)):
        raise ValueError("branch_current_a must contain one finite value per branch")
    nonlinear = any(
        row["kind"] == "nonlinear_bh"
        for row in operators["material_contract"]["materials"].values()
    )
    if nonlinear:
        solved = solve_vol2d_nonlinear_static(request)
        state = np.asarray(solved["field_state"], dtype=float)
        metadata = {
            "field_solver": "nonlinear_picard",
            "field_residual_inf": solved["residual"]["field_inf"],
            "nonlinear_iterations": solved["iterations"],
            "coenergy": solved["magnetic_energy_j"],
        }
    else:
        stiffness = np.asarray(operators["assembly"]["field_matrix"], dtype=float)
        rhs = source @ currents
        try:
            state = np.linalg.solve(stiffness, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError("linear force field matrix is singular") from exc
        residual = stiffness @ state - rhs
        metadata = {
            "field_solver": "linear_static",
            "field_residual_inf": float(np.linalg.norm(residual, ord=np.inf)),
            "nonlinear_iterations": 0,
            "coenergy": float(0.5 * state @ rhs),
        }
    metadata["field_state_sha256"] = _sha(state.tolist())
    return state, metadata


def solve_vol2d_force(request: Mapping[str, Any]) -> dict[str, Any]:
    """Solve one static field and extract a boundary-aware force resultant."""

    started = time.perf_counter()
    operators = assemble_vol2d_dynamics(request)
    assembly_s = time.perf_counter() - started
    prepared, mesh_view, materials = _prepare_request(request)
    mesh, fes, _ = _open_space(prepared, mesh_view)
    target_bundle = normalize_vol2d_force_target(
        mesh_view,
        materials,
        operators["operator_sha256"],
        request.get("force_target"),
        formulation=prepared["formulation"],
        dirichlet_boundaries=prepared.get("dirichlet_boundaries", []),
    )
    target = target_bundle["target"]

    solve_started = time.perf_counter()
    state, state_meta = _field_state(request, operators)
    from ngsolve import CoefficientFunction, GridFunction, grad, x  # type: ignore

    gfu = GridFunction(fes)
    vector = gfu.vec.FV().NumPy()
    vector[:] = 0.0
    free = np.asarray(operators["assembly"]["free_dof_indices_0based"], dtype=int)
    vector[free] = state
    if prepared["formulation"] == "planar":
        field = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    else:
        field = CoefficientFunction((-grad(gfu)[1], grad(gfu)[0] + gfu / x))
    solve_s = time.perf_counter() - solve_started

    post_started = time.perf_counter()
    if prepared["formulation"] == "planar":
        weighted = np.asarray(
            eggshell_force_2d(
                field,
                mesh,
                target["center_m"],
                target["inner_radius_m"],
                target["outer_radius_m"],
                target["air_material"],
            ),
            dtype=float,
        )
        depth = float(target["model_depth_m"])
        force_payload: dict[str, Any] = {
            "weighted_stress_n_per_m": weighted.tolist(),
            "weighted_stress_integrated_n": (weighted * depth).tolist(),
            "lorentz_scope": "not_requested",
        }
        if target["method"] == "dual_lorentz_weighted_stress":
            order = operators["assembly"]["branch_order"]
            try:
                branch_index = order.index(target["target_branch"])
            except ValueError as exc:
                raise ValueError("target_branch is absent from assembled branches") from exc
            if operators["assembly"]["branch_materials"][branch_index] != target["target_material"]:
                raise ValueError("target_branch material does not match target_material")
            currents = np.asarray(request["branch_current_a"], dtype=float)
            area = float(operators["assembly"]["branch_area_m2"][branch_index])
            turns = float(operators["assembly"]["branch_turns"][branch_index])
            jz = mesh.MaterialCF(
                {target["target_material"]: turns * currents[branch_index] / area},
                default=0.0,
            )
            lorentz = np.asarray(
                lorentz_force_2d(jz, field, mesh, target["target_material"]), dtype=float
            )
            scale = max(float(np.linalg.norm(weighted)), float(np.linalg.norm(lorentz)))
            if scale <= 1.0e-14:
                raise ValueError("dual force comparison is indeterminate for near-zero forces")
            disagreement = float(np.linalg.norm(weighted - lorentz) / scale)
            if disagreement > target["agreement_relative_tolerance"]:
                raise ValueError(
                    "Lorentz and weighted-stress force disagree: "
                    f"relative={disagreement:.6g}"
                )
            force_payload.update(
                {
                    "lorentz_n_per_m": lorentz.tolist(),
                    "lorentz_integrated_n": (lorentz * depth).tolist(),
                    "dual_method_relative_disagreement": disagreement,
                    "lorentz_scope": "total_conductor_force_only",
                }
            )
        elif target["target_material"] not in operators["assembly"]["branch_materials"]:
            force_payload["lorentz_scope"] = "not_applicable_to_passive_magnetic_target"
    else:
        axial = float(
            eggshell_force_axi(
                field,
                mesh,
                target["center_m"],
                target["inner_radius_m"],
                target["outer_radius_m"],
                target["air_material"],
            )
        )
        force_payload = {
            "weighted_stress_full_revolution_axial_n": axial,
            "net_radial_force_n": 0.0,
            "radial_force_claim_forbidden": True,
            "toroidal_weight": "2*pi*r",
            "lorentz_scope": "not_used_for_total_axisymmetric_target_force",
        }
    post_s = time.perf_counter() - post_started

    result = {
        "execution_version": {
            "radia_mcp": _package_version(),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema": FORCE_SCHEMA,
        "status": "solved",
        "operation": "solve",
        "formulation": prepared["formulation"],
        "element_family": prepared["element_family"],
        "mesh_contract_sha256": operators["assembly"]["mesh_contract"]["contract_sha256"],
        "material_contract_sha256": materials["contract_sha256"],
        "operator_sha256": operators["operator_sha256"],
        "field_state_sha256": state_meta["field_state_sha256"],
        "boundary_contract_sha256": target["boundary_contract_sha256"],
        "force_request_sha256": target["force_request_sha256"],
        "force_target": target,
        "boundary_contract": target_bundle["boundary_contract"],
        "field": state_meta,
        "force": force_payload,
        "result_output_schema_id": "radia.vol2d-force-table.v1",
        "result_output_columns": ["force_component", "value"],
        "result_output_units": target["force_unit_basis"],
        "timing_breakdown_s": {
            "assemble": assembly_s,
            "solve": solve_s,
            "force_postprocess": post_s,
            "total": time.perf_counter() - started,
        },
    }
    return result


def vol2d_force_virtual_work_gate(request: Mapping[str, Any]) -> dict[str, Any]:
    """Compare a centered fixed-current coenergy derivative with weighted stress."""

    rows = request.get("rows")
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)) or len(rows) != 3:
        raise ValueError("virtual-work rows must contain exactly three displacement samples")
    parsed = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"rows[{index}] must be an object")
        parsed.append(
            {
                "displacement_m": _finite(row.get("displacement_m"), "displacement_m"),
                "coenergy": _finite(row.get("coenergy"), "coenergy"),
                "physics_identity_sha256": str(row.get("physics_identity_sha256", "")),
                "fixed_current_sha256": str(row.get("fixed_current_sha256", "")),
                "coenergy_unit": str(row.get("coenergy_unit", "")),
            }
        )
    parsed.sort(key=lambda row: row["displacement_m"])
    positions = np.asarray([row["displacement_m"] for row in parsed])
    if not positions[0] < positions[1] < positions[2]:
        raise ValueError("virtual-work displacements must be strictly increasing")
    left = positions[1] - positions[0]
    right = positions[2] - positions[1]
    if not math.isclose(left, right, rel_tol=1.0e-8, abs_tol=1.0e-15):
        raise ValueError("virtual-work samples must be centered and equally spaced")
    for key in ("physics_identity_sha256", "fixed_current_sha256", "coenergy_unit"):
        values = {row[key] for row in parsed}
        if len(values) != 1 or not next(iter(values)):
            raise ValueError(f"virtual-work rows require one non-empty {key}")
    unit = parsed[0]["coenergy_unit"]
    if unit not in {"J", "J_per_m"}:
        raise ValueError("coenergy_unit must be J or J_per_m")
    derivative = (parsed[2]["coenergy"] - parsed[0]["coenergy"]) / (positions[2] - positions[0])
    weighted = _finite(request.get("weighted_stress_force"), "weighted_stress_force")
    scale = max(abs(derivative), abs(weighted))
    if scale <= 1.0e-14:
        raise ValueError("virtual-work comparison is indeterminate for near-zero force")
    relative = abs(derivative - weighted) / scale
    tolerance = _positive(request.get("relative_tolerance", 0.2), "relative_tolerance")
    if relative > tolerance:
        raise ValueError(f"virtual work and weighted stress disagree: relative={relative:.6g}")
    return {
        "schema": VIRTUAL_WORK_SCHEMA,
        "status": "ok",
        "force_from_coenergy_derivative": derivative,
        "weighted_stress_force": weighted,
        "force_unit": "N" if unit == "J" else "N_per_m",
        "relative_disagreement": relative,
        "center_displacement_m": float(positions[1]),
        "spacing_m": float(left),
        "fixed_current": True,
        "physics_identity_sha256": parsed[0]["physics_identity_sha256"],
        "fixed_current_sha256": parsed[0]["fixed_current_sha256"],
    }


def vol2d_force_refinement_gate(request: Mapping[str, Any]) -> dict[str, Any]:
    """Require a nonzero, convergent force sequence with one physics identity."""

    rows = request.get("rows")
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)) or len(rows) < 2:
        raise ValueError("force refinement requires at least two rows")
    parsed = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"rows[{index}] must be an object")
        vector = _vector2(row.get("force_vector"), "force_vector")
        parsed.append(
            {
                "mesh_cells": int(row.get("mesh_cells", 0)),
                "force_vector": vector,
                "force_unit": str(row.get("force_unit", "")),
                "physics_identity_sha256": str(row.get("physics_identity_sha256", "")),
            }
        )
    if any(row["mesh_cells"] <= 0 for row in parsed):
        raise ValueError("mesh_cells must be positive")
    if any(right["mesh_cells"] <= left["mesh_cells"] for left, right in zip(parsed, parsed[1:])):
        raise ValueError("mesh refinement rows must have increasing cell counts")
    for key in ("force_unit", "physics_identity_sha256"):
        values = {row[key] for row in parsed}
        if len(values) != 1 or not next(iter(values)):
            raise ValueError(f"force refinement rows require one non-empty {key}")
    norms = [float(np.linalg.norm(row["force_vector"])) for row in parsed]
    minimum = _positive(request.get("minimum_force_magnitude", 1.0e-12), "minimum_force_magnitude")
    if max(norms) <= minimum:
        raise ValueError("force refinement cannot use a near-zero reference")
    relative_changes = []
    for left, right in zip(parsed, parsed[1:]):
        scale = max(np.linalg.norm(left["force_vector"]), np.linalg.norm(right["force_vector"]), minimum)
        relative_changes.append(
            float(np.linalg.norm(np.subtract(right["force_vector"], left["force_vector"])) / scale)
        )
    tolerance = _positive(request.get("terminal_relative_tolerance", 0.1), "terminal_relative_tolerance")
    if relative_changes[-1] > tolerance:
        raise ValueError(
            f"terminal force refinement is not converged: relative={relative_changes[-1]:.6g}"
        )
    return {
        "schema": FORCE_REFINEMENT_SCHEMA,
        "status": "ok",
        "row_count": len(parsed),
        "relative_changes": relative_changes,
        "terminal_relative_change": relative_changes[-1],
        "force_unit": parsed[0]["force_unit"],
        "physics_identity_sha256": parsed[0]["physics_identity_sha256"],
    }


def analyze_vol2d_force(request: Mapping[str, Any]) -> dict[str, Any]:
    """Closed operation router used by the owned MCP worker."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    operation = str(request.get("operation", "solve"))
    if operation == "solve":
        return solve_vol2d_force(request)
    if operation == "virtual_work_gate":
        return vol2d_force_virtual_work_gate(request)
    if operation == "refinement_gate":
        return vol2d_force_refinement_gate(request)
    raise ValueError("operation must be solve, virtual_work_gate, or refinement_gate")
