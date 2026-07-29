"""Multi-conductor electrostatics and Maxwell-force artifacts for 2-D ``.vol``."""

from __future__ import annotations

import csv
import io
import time
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from typing import Any

import numpy as np

from .vol2d_circuit import parse_netgen_2d_vol
from .vol2d_postprocess import _canonical, _export_entry, _package_version, _sha
from .vol2d_scalar import _SAFE_NAME, _finite, solve_vol2d_scalar


ELECTROSTATIC_SYSTEM_SCHEMA = "radia.vol2d-electrostatic-system.v1"


def _names(raw: Any, label: str, *, minimum: int = 1) -> list[str]:
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
        raise ValueError(f"{label} must be a sequence of boundary names")
    names = [str(value) for value in raw]
    if len(names) < minimum or any(not name for name in names):
        raise ValueError(f"{label} must contain at least {minimum} non-empty names")
    if len(set(names)) != len(names):
        raise ValueError(f"{label} must not contain duplicates")
    return names


def _matrix_csv(names: Sequence[str], matrix: np.ndarray) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["conductor", *names])
    for name, row in zip(names, matrix):
        writer.writerow([name, *row.tolist()])
    return output.getvalue()


def solve_vol2d_electrostatic_system(request: Mapping[str, Any]) -> dict[str, Any]:
    """Build a Maxwell capacitance matrix and one loaded electrostatic state.

    Every conductor is excited to one volt while all peers are grounded.  The
    boundary reaction vectors form ``C`` directly from the common FEM operator.
    A final requested voltage vector supplies charge, energy, and optional
    Maxwell traction on named conductor boundaries.
    """

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    started = time.perf_counter()
    text = request.get("vol_text")
    if not isinstance(text, str) or not text.strip():
        raise ValueError("vol_text must be non-empty")
    mesh_view = parse_netgen_2d_vol(
        text, source_name=str(request.get("source_name", "electrostatic_system.vol"))
    )
    boundary_names = set(mesh_view.contract()["boundary_names"])
    conductors = _names(request.get("conductors"), "conductors", minimum=2)
    unknown = sorted(set(conductors) - boundary_names)
    if unknown:
        raise ValueError(f"unknown conductor boundaries: {unknown}")
    raw_force_boundaries = request.get("force_boundaries")
    if raw_force_boundaries is None or (
        isinstance(raw_force_boundaries, Sequence)
        and not isinstance(raw_force_boundaries, (str, bytes))
        and len(raw_force_boundaries) == 0
    ):
        raw_force_boundaries = conductors[:1]
    force_boundaries = _names(raw_force_boundaries, "force_boundaries")
    force_unknown = sorted(set(force_boundaries) - set(conductors))
    if force_unknown:
        raise ValueError(
            f"force_boundaries must be selected from conductors: {force_unknown}"
        )
    raw_voltages = request.get("applied_voltages_v")
    if not isinstance(raw_voltages, Sequence) or isinstance(raw_voltages, (str, bytes)):
        raise ValueError("applied_voltages_v must contain one value per conductor")
    voltages = np.asarray(
        [_finite(value, "applied_voltages_v") for value in raw_voltages], dtype=float
    )
    if voltages.shape != (len(conductors),):
        raise ValueError("applied_voltages_v must contain one value per conductor")
    if np.ptp(voltages) <= 0.0:
        raise ValueError("applied_voltages_v must contain a nonzero voltage difference")
    basename = str(request.get("export_basename", "vol2d_electrostatic_system")).strip()
    if not _SAFE_NAME.fullmatch(basename):
        raise ValueError("export_basename must be a portable filename stem")
    materials = request.get("materials")
    if not isinstance(materials, Mapping):
        raise ValueError("materials must be an object")
    for name, row in materials.items():
        if not isinstance(row, Mapping):
            raise ValueError(f"materials[{name}] must be an object")
        if abs(_finite(row.get("volumetric_source_si", 0.0), f"materials[{name}].volumetric_source_si")) > 0.0:
            raise ValueError("capacitance matrices require zero volumetric charge source")

    base = {
        key: value
        for key, value in request.items()
        if key
        not in {
            "operation",
            "conductors",
            "force_boundaries",
            "applied_voltages_v",
        }
    }
    base.update({"physics": "electrostatic", "frequency_hz": 0.0})
    columns: list[list[float]] = []
    excitation_hashes: list[str] = []
    for column, excited in enumerate(conductors):
        one_hot = dict(base)
        one_hot["dirichlet_values"] = {
            name: float(index == column) for index, name in enumerate(conductors)
        }
        one_hot["force_boundaries"] = []
        one_hot["terminal_pair"] = None
        one_hot["export_basename"] = f"{basename}_basis_{column + 1}"
        solved = solve_vol2d_scalar(one_hot)
        reactions = solved["result_contract"]["observables"]["boundary_reactions"]
        columns.append([float(reactions[name][0]) for name in conductors])
        excitation_hashes.append(
            str(solved["result_contract"]["field_state_sha256"])
        )
    matrix = np.asarray(columns, dtype=float).T
    scale = max(1.0e-30, float(np.max(np.abs(matrix))))
    symmetry_error = float(np.max(np.abs(matrix - matrix.T))) / scale
    symmetric = 0.5 * (matrix + matrix.T)
    eigenvalues = np.linalg.eigvalsh(symmetric)
    off_diagonal = matrix - np.diag(np.diag(matrix))
    row_sum_error = float(np.max(np.abs(np.sum(matrix, axis=1)))) / scale
    checks = {
        "reciprocal": symmetry_error <= 1.0e-8,
        "positive_semidefinite": float(eigenvalues[0]) >= -1.0e-9 * scale,
        "positive_self_capacitance": bool(np.all(np.diag(matrix) > 0.0)),
        "nonpositive_mutual_capacitance": bool(np.max(off_diagonal) <= 1.0e-9 * scale),
        "charge_neutral_under_common_mode": row_sum_error <= 1.0e-7,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ValueError(f"electrostatic capacitance gate failed: {failed}")

    applied = dict(base)
    applied["dirichlet_values"] = {
        name: float(value) for name, value in zip(conductors, voltages)
    }
    applied["force_boundaries"] = force_boundaries
    applied["terminal_pair"] = None
    applied["export_basename"] = basename
    loaded = solve_vol2d_scalar(applied)
    loaded_observables = loaded["result_contract"]["observables"]
    charges = matrix @ voltages
    matrix_energy = 0.5 * float(voltages @ matrix @ voltages)
    field_energy = float(loaded_observables["electric_energy_j"])
    energy_relative_error = abs(matrix_energy - field_energy) / max(
        1.0e-30, abs(matrix_energy), abs(field_energy)
    )
    if energy_relative_error > 1.0e-8:
        raise ValueError(
            "capacitance energy does not match the loaded field energy: "
            f"relative_error={energy_relative_error:.3e}"
        )
    raw_forces = loaded_observables.get("electrostatic_force_by_boundary_n", {})
    conductor_forces = {
        name: [-float(value[0]), -float(value[1])]
        for name, value in raw_forces.items()
    }
    prepare_solve_s = time.perf_counter() - started
    request_contract = {
        "schema": "radia.vol2d-electrostatic-system-request.v1",
        "mesh_contract_sha256": mesh_view.contract()["contract_sha256"],
        "conductors": conductors,
        "force_boundaries": force_boundaries,
        "applied_voltages_v": voltages.tolist(),
        "element_family": str(request.get("element_family", "P1")),
        "formulation": str(request.get("formulation", "planar")),
        "material_contract_sha256": loaded["result_contract"][
            "material_contract_sha256"
        ],
        "operator_sha256": loaded["result_contract"]["operator_sha256"],
    }
    result_contract = {
        "schema": ELECTROSTATIC_SYSTEM_SCHEMA,
        "status": "solved",
        "request_contract": request_contract,
        "request_contract_sha256": _sha(request_contract),
        "conductor_order": conductors,
        "capacitance_matrix_f": matrix.tolist(),
        "capacitance_matrix_sha256": _sha(matrix.tolist()),
        "basis_field_state_sha256": excitation_hashes,
        "applied_voltage_v": voltages.tolist(),
        "terminal_charge_c": charges.tolist(),
        "stored_energy_j": field_energy,
        "matrix_energy_j": matrix_energy,
        "energy_relative_error": energy_relative_error,
        "electrostatic_force_on_conductor_n": conductor_forces,
        "force_convention": "negative dielectric-domain Maxwell traction",
        "checks": checks,
        "metrics": {
            "symmetry_relative_error": symmetry_error,
            "minimum_eigenvalue_f": float(eigenvalues[0]),
            "common_mode_row_sum_relative_error": row_sum_error,
        },
        "loaded_field_state_sha256": loaded["result_contract"]["field_state_sha256"],
        "operator_sha256": loaded["result_contract"]["operator_sha256"],
        "material_contract_sha256": loaded["result_contract"][
            "material_contract_sha256"
        ],
        "generated_vol_git_required": False,
    }
    json_content = _canonical(result_contract)
    csv_content = _matrix_csv(conductors, matrix)
    exports = {
        "json": _export_entry(json_content, f"{basename}.json", "application/json"),
        "csv": _export_entry(csv_content, f"{basename}.csv", "text/csv"),
    }
    for name in ("gmsh_msh", "gmsh_geo", "gmsh_geo_opt", "gmsh_msh_opt"):
        exports[name] = loaded["exports"][name]
    result_contract["export_content_sha256"] = {
        name: row["sha256"] for name, row in exports.items()
    }
    total = time.perf_counter() - started
    return {
        "schema": ELECTROSTATIC_SYSTEM_SCHEMA,
        "status": "solved",
        "operation": "electrostatic_system",
        "execution_version": {
            "radia_mcp": _package_version(),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "timing_s": {
            "basis_and_loaded_solves": prepare_solve_s,
            "contract_and_export": total - prepare_solve_s,
            "total": total,
        },
        "result_contract": result_contract,
        "exports": exports,
    }
