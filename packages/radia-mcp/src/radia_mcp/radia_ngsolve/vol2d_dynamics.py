"""Nonlinear and transient dimension-2 ``.vol`` field execution.

The linear mesh/circuit bridge lives in :mod:`vol2d_circuit`.  This module
adds an explicit material contract, conductivity mass assembly, a readable
Picard solve, and a theta-method transient.  Generated meshes remain replay
artifacts identified by digest; they are not package fixtures.
"""

from __future__ import annotations

import hashlib
from importlib import metadata
import json
import math
import time
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from typing import Any

import numpy as np

from .circuit_system import compile_circuit_state_space
from .vol2d_circuit import (
    MU0,
    _dense_matrix,
    _family_contract,
    _runtime_vol_path,
    assemble_vol2d_field,
    parse_netgen_2d_vol,
)


MATERIAL_SCHEMA = "radia.vol2d-material-contract.v1"
DYNAMIC_SCHEMA = "radia.vol2d-dynamic-analysis.v1"
_OPERATIONS = {
    "assemble",
    "nonlinear_static",
    "harmonic",
    "transient",
    "state_space",
}


def _package_version() -> str:
    try:
        return metadata.version("radia-mcp")
    except metadata.PackageNotFoundError:
        return "unknown"


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


def _complex_value(value: Any, label: str) -> complex:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 2:
            raise ValueError(f"{label} complex pair must contain [real, imag]")
        result = complex(_finite(value[0], label), _finite(value[1], label))
    else:
        result = complex(_finite(value, label), 0.0)
    if not (math.isfinite(result.real) and math.isfinite(result.imag)):
        raise ValueError(f"{label} must be finite")
    return result


def _complex_pair(value: complex | float) -> list[float]:
    parsed = complex(value)
    return [float(parsed.real), float(parsed.imag)]


def _sha256_json(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _matrix_sha256(*matrices: np.ndarray, metadata: Mapping[str, Any]) -> str:
    digest = hashlib.sha256(
        json.dumps(dict(metadata), sort_keys=True, separators=(",", ":")).encode("utf-8")
    )
    for matrix in matrices:
        array = np.ascontiguousarray(np.asarray(matrix, dtype="<f8"))
        digest.update(str(array.shape).encode("ascii"))
        digest.update(array.tobytes())
    return digest.hexdigest()


def _normalize_bh_curve(raw: Any, label: str) -> list[dict[str, float]]:
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or len(raw) < 3:
        raise ValueError(f"{label} must contain at least three B-H rows")
    rows: list[dict[str, float]] = []
    for index, row in enumerate(raw):
        if not isinstance(row, Mapping):
            raise ValueError(f"{label}[{index}] must be an object")
        rows.append(
            {
                "b_t": _finite(row.get("b_t"), f"{label}[{index}].b_t"),
                "h_a_per_m": _finite(
                    row.get("h_a_per_m"), f"{label}[{index}].h_a_per_m"
                ),
            }
        )
    b_values = [row["b_t"] for row in rows]
    h_values = [row["h_a_per_m"] for row in rows]
    if abs(b_values[0]) > 1.0e-15 or abs(h_values[0]) > 1.0e-12:
        raise ValueError(f"{label} must start at the B=0, H=0 origin")
    if any(right <= left for left, right in zip(b_values, b_values[1:])):
        raise ValueError(f"{label} B values must be strictly increasing")
    if any(right <= left for left, right in zip(h_values, h_values[1:])):
        raise ValueError(f"{label} H values must be strictly increasing")
    return rows


def normalize_vol2d_materials(
    mesh_contract: Mapping[str, Any], raw: Any
) -> dict[str, Any]:
    """Validate one complete isotropic material law per mesh material."""

    names = mesh_contract.get("material_names")
    if not isinstance(names, list) or not names:
        raise ValueError("mesh_contract material_names must be a non-empty list")
    if not isinstance(raw, Mapping):
        raise ValueError("materials must be an object keyed by mesh material")
    supplied = {str(key) for key in raw}
    expected = set(names)
    if supplied != expected:
        raise ValueError(
            "materials must cover mesh materials exactly: "
            f"missing={sorted(expected - supplied)}, unknown={sorted(supplied - expected)}"
        )

    normalized: dict[str, dict[str, Any]] = {}
    for name in names:
        law = raw[name]
        if not isinstance(law, Mapping):
            raise ValueError(f"materials[{name}] must be an object")
        conductivity = _finite(
            law.get("conductivity_s_per_m", 0.0),
            f"materials[{name}].conductivity_s_per_m",
        )
        if conductivity < 0.0:
            raise ValueError(f"materials[{name}].conductivity_s_per_m must be nonnegative")
        has_bh = "bh_curve" in law
        has_mu = "permeability_h_per_m" in law
        if has_bh == has_mu:
            raise ValueError(
                f"materials[{name}] must define exactly one of bh_curve or permeability_h_per_m"
            )
        if has_bh:
            rows = _normalize_bh_curve(law["bh_curve"], f"materials[{name}].bh_curve")
            first_db = rows[1]["b_t"] - rows[0]["b_t"]
            first_dh = rows[1]["h_a_per_m"] - rows[0]["h_a_per_m"]
            entry = {
                "kind": "nonlinear_bh",
                "interpolation": "piecewise_linear_H_of_B",
                "bh_curve": rows,
                "initial_permeability_h_per_m": first_db / first_dh,
                "conductivity_s_per_m": conductivity,
            }
        else:
            entry = {
                "kind": "linear",
                "permeability_h_per_m": _positive(
                    law["permeability_h_per_m"],
                    f"materials[{name}].permeability_h_per_m",
                ),
                "conductivity_s_per_m": conductivity,
            }
        entry["material_sha256"] = _sha256_json(entry)
        normalized[name] = entry

    contract = {
        "schema": MATERIAL_SCHEMA,
        "units": {
            "b": "T",
            "h": "A/m",
            "permeability": "H/m",
            "conductivity": "S/m",
        },
        "mesh_contract_sha256": mesh_contract.get("contract_sha256"),
        "materials": normalized,
    }
    contract["contract_sha256"] = _sha256_json(contract)
    return contract


def _initial_permeability(material: Mapping[str, Any]) -> float:
    if material["kind"] == "linear":
        return float(material["permeability_h_per_m"])
    return float(material["initial_permeability_h_per_m"])


def _prepare_request(request: Mapping[str, Any]) -> tuple[dict[str, Any], Any, dict[str, Any]]:
    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    prepared = dict(request)
    mesh_view = parse_netgen_2d_vol(
        prepared.get("vol_text"), source_name=str(prepared.get("source_name", "generated.vol"))
    )
    mesh_contract = mesh_view.contract()
    material_contract = normalize_vol2d_materials(mesh_contract, prepared.get("materials"))
    prepared["permeability_h_per_m"] = {
        name: _initial_permeability(material)
        for name, material in material_contract["materials"].items()
    }
    return prepared, mesh_view, material_contract


def _open_space(
    prepared: Mapping[str, Any],
    mesh_view: Any,
    *,
    complex_space: bool = False,
) -> tuple[Any, Any, dict[str, Any]]:
    try:
        from ngsolve import H1, Mesh  # type: ignore
    except ImportError as exc:
        raise RuntimeError("NGSolve is required for dimension-2 dynamics") from exc
    family = str(prepared.get("element_family", ""))
    formulation = str(prepared.get("formulation", ""))
    options = _family_contract(mesh_view, family, formulation)
    boundaries = prepared.get("dirichlet_boundaries")
    if not isinstance(boundaries, Sequence) or isinstance(boundaries, (str, bytes)):
        raise ValueError("dirichlet_boundaries must be a non-empty sequence")
    dirichlet = "|".join(str(value) for value in boundaries)
    path = _runtime_vol_path(prepared["vol_text"], mesh_view.content_sha256)
    mesh = Mesh(str(path))
    if formulation == "planar":
        fes = H1(
            mesh,
            order=options["order"],
            dirichlet=dirichlet,
            complex=complex_space,
        )
    elif formulation == "axisymmetric_henrotte":
        try:
            from radia.axifem import H1Henrotte
        except ImportError as exc:
            raise RuntimeError("Radia axifem extension is required") from exc
        fes = H1Henrotte(
            mesh,
            order=options["order"],
            dirichlet=dirichlet,
            curvedquad=options["curvedquad"],
            complex=complex_space,
        )
    else:
        raise ValueError("formulation must be planar or axisymmetric_henrotte")
    return mesh, fes, options


def assemble_vol2d_dynamics(request: Mapping[str, Any]) -> dict[str, Any]:
    """Assemble constrained stiffness, conductivity mass, and branch sources."""

    prepared, mesh_view, materials = _prepare_request(request)
    field = assemble_vol2d_field(prepared)
    mesh, fes, _ = _open_space(prepared, mesh_view)
    conductivity = {
        name: material["conductivity_s_per_m"]
        for name, material in materials["materials"].items()
    }
    sigma_cf = mesh.MaterialCF(conductivity)
    try:
        from ngsolve import BilinearForm, dx  # type: ignore
    except ImportError as exc:
        raise RuntimeError("NGSolve is required for dimension-2 dynamics") from exc
    mass = BilinearForm(fes, symmetric=True)
    if prepared["formulation"] == "planar":
        trial, test = fes.TnT()
        mass += sigma_cf * trial * test * dx
    else:
        try:
            from radia.axifem import AxiHenrotteSigmaMassBFI
        except ImportError as exc:
            raise RuntimeError("Radia axifem extension is required") from exc
        mass += AxiHenrotteSigmaMassBFI(sigma_cf)
    mass.Assemble()
    free = np.asarray(field["free_dof_indices_0based"], dtype=int)
    mass_matrix = _dense_matrix(mass.mat)[np.ix_(free, free)]
    mass_matrix = 0.5 * (mass_matrix + mass_matrix.T)
    scale = max(1.0, float(np.max(np.abs(mass_matrix))))
    symmetry_error = float(np.max(np.abs(mass_matrix - mass_matrix.T)))
    eigenvalues = np.linalg.eigvalsh(mass_matrix)
    if eigenvalues[0] < -1.0e-10 * scale:
        raise ValueError("conductivity mass matrix is not positive semidefinite")
    stiffness = np.asarray(field["field_matrix"], dtype=float)
    sources = np.asarray(field["source_matrix"], dtype=float)
    operator_sha = _matrix_sha256(
        stiffness,
        mass_matrix,
        sources,
        metadata={
            "mesh": field["mesh_contract"]["contract_sha256"],
            "materials": materials["contract_sha256"],
            "family": field["element_family"],
            "formulation": field["formulation"],
        },
    )
    expected = request.get("expected_operator_sha256")
    if expected is not None and str(expected) != operator_sha:
        raise ValueError("expected_operator_sha256 does not match assembled operators")
    return {
        "schema": DYNAMIC_SCHEMA,
        "status": "assembled",
        "assembly": field,
        "material_contract": materials,
        "conductivity_mass_matrix": mass_matrix.tolist(),
        "mass_matrix_symmetry_error": symmetry_error,
        "mass_matrix_minimum_eigenvalue": float(eigenvalues[0]),
        "mass_matrix_rank": int(np.linalg.matrix_rank(mass_matrix)),
        "operator_sha256": operator_sha,
    }


def _piecewise_h_of_b(bmag: Any, rows: Sequence[Mapping[str, float]]) -> Any:
    from ngsolve import IfPos  # type: ignore

    def segment(left: int, right: int) -> Any:
        b0, h0 = rows[left]["b_t"], rows[left]["h_a_per_m"]
        b1, h1 = rows[right]["b_t"], rows[right]["h_a_per_m"]
        return h0 + (h1 - h0) * (bmag - b0) / (b1 - b0)

    value = segment(0, 1)
    for left in range(1, len(rows) - 1):
        value = IfPos(bmag - rows[left]["b_t"], segment(left, left + 1), value)
    last = len(rows) - 1
    tail = segment(last - 1, last)
    return IfPos(bmag - rows[last]["b_t"], tail, value)


def _constitutive_cf(mesh: Any, gfu: Any, materials: Mapping[str, Any], formulation: str) -> Any:
    from ngsolve import CoefficientFunction, InnerProduct, IfPos, grad, sqrt, x  # type: ignore

    if formulation == "planar":
        field = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    else:
        field = CoefficientFunction((grad(gfu)[0] + gfu / x, -grad(gfu)[1]))
    bmag = sqrt(InnerProduct(field, field) + 1.0e-30)
    coefficient = 0.0
    for name, material in materials["materials"].items():
        indicator = mesh.MaterialCF({name: 1.0}, default=0.0)
        if material["kind"] == "linear":
            mu = material["permeability_h_per_m"]
            local = 1.0 / mu if formulation == "planar" else mu
        else:
            rows = material["bh_curve"]
            hmag = _piecewise_h_of_b(bmag, rows)
            nu0 = 1.0 / material["initial_permeability_h_per_m"]
            nu = IfPos(bmag - 1.0e-12, hmag / (bmag + 1.0e-30), nu0)
            local = nu if formulation == "planar" else 1.0 / nu
        coefficient += indicator * local
    return coefficient


def _assemble_nonlinear_stiffness(
    prepared: Mapping[str, Any], mesh: Any, fes: Any, gfu: Any, materials: Mapping[str, Any]
) -> np.ndarray:
    from ngsolve import BilinearForm, dx, grad  # type: ignore

    coefficient = _constitutive_cf(mesh, gfu, materials, prepared["formulation"])
    stiffness = BilinearForm(fes, symmetric=True)
    if prepared["formulation"] == "planar":
        trial, test = fes.TnT()
        stiffness += coefficient * grad(trial) * grad(test) * dx
    else:
        from radia.axifem import AxiHenrotteStiffnessBFI

        stiffness += AxiHenrotteStiffnessBFI(coefficient)
    stiffness.Assemble()
    free = np.flatnonzero(np.asarray(fes.FreeDofs(), dtype=bool))
    matrix = _dense_matrix(stiffness.mat)[np.ix_(free, free)]
    return 0.5 * (matrix + matrix.T)


def solve_vol2d_nonlinear_static(request: Mapping[str, Any]) -> dict[str, Any]:
    """Solve one prescribed-current nonlinear magnetostatic operating point."""

    operators = assemble_vol2d_dynamics(request)
    prepared, mesh_view, materials = _prepare_request(request)
    mesh, fes, _ = _open_space(prepared, mesh_view)
    currents = np.asarray(request.get("branch_current_a"), dtype=float)
    sources = np.asarray(operators["assembly"]["source_matrix"], dtype=float)
    if currents.ndim != 1 or currents.size != sources.shape[1] or not np.all(np.isfinite(currents)):
        raise ValueError("branch_current_a must contain one finite real value per branch")
    rhs = sources @ currents
    relaxation = _positive(request.get("relaxation", 0.5), "relaxation")
    if relaxation > 1.0:
        raise ValueError("relaxation must not exceed one")
    tolerance = _positive(request.get("relative_tolerance", 1.0e-7), "relative_tolerance")
    maximum_iterations = int(request.get("maximum_iterations", 80))
    if maximum_iterations < 2 or maximum_iterations > 500:
        raise ValueError("maximum_iterations must be in [2, 500]")

    from ngsolve import GridFunction  # type: ignore

    gfu = GridFunction(fes)
    state = np.zeros(sources.shape[0], dtype=float)
    free = np.asarray(operators["assembly"]["free_dof_indices_0based"], dtype=int)
    history: list[dict[str, float]] = []
    converged = False
    for iteration in range(1, maximum_iterations + 1):
        stiffness = _assemble_nonlinear_stiffness(prepared, mesh, fes, gfu, materials)
        try:
            candidate = np.linalg.solve(stiffness, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError("nonlinear stiffness matrix is singular") from exc
        updated = (1.0 - relaxation) * state + relaxation * candidate
        relative_change = float(
            np.linalg.norm(updated - state) / max(np.linalg.norm(updated), 1.0e-30)
        )
        state = updated
        vector = gfu.vec.FV().NumPy()
        vector[:] = 0.0
        vector[free] = state
        residual = float(np.linalg.norm(stiffness @ state - rhs, ord=np.inf))
        history.append(
            {
                "iteration": float(iteration),
                "relative_state_change": relative_change,
                "field_residual_inf": residual,
            }
        )
        if iteration >= 3 and relative_change <= tolerance:
            converged = True
            break
    if not converged:
        raise ValueError(
            f"nonlinear Picard iteration did not converge in {maximum_iterations} iterations"
        )
    stiffness = _assemble_nonlinear_stiffness(prepared, mesh, fes, gfu, materials)
    residual = stiffness @ state - rhs
    result = {
        "schema": DYNAMIC_SCHEMA,
        "status": "solved",
        "operation": "nonlinear_static",
        "mesh_contract_sha256": operators["assembly"]["mesh_contract"]["contract_sha256"],
        "material_contract_sha256": materials["contract_sha256"],
        "linear_operator_sha256": operators["operator_sha256"],
        "nonlinear_operator_sha256": _matrix_sha256(
            stiffness,
            metadata={"materials": materials["contract_sha256"], "state": state.tolist()},
        ),
        "branch_current_a": currents.tolist(),
        "field_state": state.tolist(),
        "field_state_l2": float(np.linalg.norm(state)),
        "magnetic_energy_j": float(0.5 * state @ stiffness @ state),
        "iterations": len(history),
        "iteration_history": history,
        "residual": {"field_inf": float(np.linalg.norm(residual, ord=np.inf))},
        "converged": True,
    }
    return result


def _solve_harmonic_matrices(
    stiffness: np.ndarray,
    mass: np.ndarray,
    sources: np.ndarray,
    *,
    frequency_hz: float,
    branch_current_a: Sequence[Any],
) -> dict[str, Any]:
    """Pure matrix kernel for one ``exp(+j*omega*t)`` eddy-current solve."""

    frequency_hz = _positive(frequency_hz, "frequency_hz")
    stiffness = np.asarray(stiffness, dtype=float)
    mass = np.asarray(mass, dtype=float)
    sources = np.asarray(sources, dtype=float)
    if stiffness.ndim != 2 or stiffness.shape[0] != stiffness.shape[1]:
        raise ValueError("stiffness must be square")
    if mass.shape != stiffness.shape:
        raise ValueError("mass must match stiffness")
    if sources.ndim != 2 or sources.shape[0] != stiffness.shape[0]:
        raise ValueError("sources row count must match stiffness")
    omega = 2.0 * math.pi * frequency_hz
    raw_currents = branch_current_a
    if (
        not isinstance(raw_currents, Sequence)
        or isinstance(raw_currents, (str, bytes))
    ):
        raise ValueError("branch_current_a must contain one phasor per branch")
    currents = np.asarray(
        [
            _complex_value(value, f"branch_current_a[{index}]")
            for index, value in enumerate(raw_currents)
        ],
        dtype=complex,
    )
    if currents.shape != (sources.shape[1],):
        raise ValueError("branch_current_a must contain one phasor per branch")
    if not np.any(np.abs(currents) > 0.0):
        raise ValueError("harmonic branch_current_a must contain a nonzero excitation")

    operator = stiffness.astype(complex) + 1j * omega * mass
    rhs = sources @ currents
    solve_started = time.perf_counter()
    try:
        state = np.linalg.solve(operator, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("harmonic eddy operator is singular") from exc
    solve_s = time.perf_counter() - solve_started
    residual = operator @ state - rhs
    flux_linkage = sources.T @ state
    branch_voltage = 1j * omega * flux_linkage
    apparent_power = 0.5 * np.vdot(currents, branch_voltage)
    magnetic_energy = 0.25 * float(np.real(np.vdot(state, stiffness @ state)))
    eddy_loss = 0.5 * omega * omega * float(np.real(np.vdot(state, mass @ state)))
    power_error = float(np.real(apparent_power) - eddy_loss)
    power_scale = max(1.0, abs(eddy_loss), abs(float(np.real(apparent_power))))
    if abs(power_error) > 1.0e-9 * power_scale:
        raise ValueError(
            "harmonic branch power does not close the conductivity loss: "
            f"error={power_error:.6e} W"
        )
    return {
        "frequency_hz": frequency_hz,
        "angular_frequency_rad_s": omega,
        "currents": currents,
        "state": state,
        "flux_linkage": flux_linkage,
        "branch_voltage": branch_voltage,
        "apparent_power": apparent_power,
        "magnetic_energy": magnetic_energy,
        "eddy_loss": eddy_loss,
        "power_error": power_error,
        "residual": residual,
        "solve_s": solve_s,
    }


def solve_vol2d_harmonic(request: Mapping[str, Any]) -> dict[str, Any]:
    """Solve ``(K + j*omega*M) a = S i`` at one steady AC frequency.

    ``M`` is the conductivity mass matrix, so ``0.5*omega^2*a^H*M*a``
    is the time-average eddy-current loss.  The branch voltage convention is
    ``v = j*omega*S^T*a`` for ``exp(+j*omega*t)`` phasors.
    """

    started = time.perf_counter()
    operators = assemble_vol2d_dynamics(request)
    assembly_s = time.perf_counter() - started
    nonlinear = [
        name
        for name, material in operators["material_contract"]["materials"].items()
        if material["kind"] != "linear"
    ]
    if nonlinear:
        raise ValueError(
            "harmonic eddy currently requires linear permeability; nonlinear "
            f"materials={nonlinear} need a separately validated harmonic iteration"
        )
    stiffness = np.asarray(operators["assembly"]["field_matrix"], dtype=float)
    mass = np.asarray(operators["conductivity_mass_matrix"], dtype=float)
    sources = np.asarray(operators["assembly"]["source_matrix"], dtype=float)
    solved = _solve_harmonic_matrices(
        stiffness,
        mass,
        sources,
        frequency_hz=request.get("frequency_hz"),
        branch_current_a=request.get("branch_current_a"),
    )
    frequency_hz = solved["frequency_hz"]
    omega = solved["angular_frequency_rad_s"]
    currents = solved["currents"]
    state = solved["state"]
    flux_linkage = solved["flux_linkage"]
    branch_voltage = solved["branch_voltage"]
    apparent_power = solved["apparent_power"]
    magnetic_energy = solved["magnetic_energy"]
    eddy_loss = solved["eddy_loss"]
    power_error = solved["power_error"]
    residual = solved["residual"]
    solve_s = solved["solve_s"]

    export_started = time.perf_counter()
    prepared, mesh_view, _materials = _prepare_request(request)
    mesh, fes, _options = _open_space(prepared, mesh_view, complex_space=True)
    from ngsolve import CoefficientFunction, GridFunction, grad, x  # type: ignore
    from .vol2d_scalar import _export_entry, _gmsh_export

    gfu = GridFunction(fes)
    free = np.asarray(operators["assembly"]["free_dof_indices_0based"], dtype=int)
    vector = gfu.vec.FV().NumPy()
    vector[:] = 0.0
    vector[free] = state
    if prepared["formulation"] == "planar":
        field = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    else:
        field = CoefficientFunction((grad(gfu)[0] + gfu / x, -grad(gfu)[1]))
    basename = str(request.get("export_basename", "vol2d_harmonic_eddy")).strip()
    gmsh = _gmsh_export(
        mesh,
        gfu,
        field,
        basename=basename,
        request_sha256=operators["operator_sha256"],
        complex_field=True,
    )
    exports = {
        "gmsh_msh": _export_entry(gmsh["gmsh_msh"], f"{basename}.msh", "model/mesh"),
        "gmsh_geo": _export_entry(gmsh["gmsh_geo"], f"{basename}.geo", "text/plain"),
        "gmsh_geo_opt": _export_entry(
            gmsh["gmsh_geo_opt"], f"{basename}.geo.opt", "text/plain"
        ),
        "gmsh_msh_opt": _export_entry(
            gmsh["gmsh_msh_opt"], f"{basename}.msh.opt", "text/plain"
        ),
    }
    export_s = time.perf_counter() - export_started
    result = {
        "schema": DYNAMIC_SCHEMA,
        "status": "solved",
        "operation": "harmonic",
        "phasor_convention": "exp(+j*omega*t), RMS branch current",
        "frequency_hz": frequency_hz,
        "angular_frequency_rad_s": omega,
        "branch_current_a": [_complex_pair(value) for value in currents],
        "field_state": [_complex_pair(value) for value in state],
        "flux_linkage_wb_turn": [_complex_pair(value) for value in flux_linkage],
        "branch_voltage_v": [_complex_pair(value) for value in branch_voltage],
        "apparent_power_va": _complex_pair(apparent_power),
        "magnetic_energy_j": magnetic_energy,
        "eddy_loss_w": eddy_loss,
        "power_closure_error_w": power_error,
        "residual_inf": float(np.linalg.norm(residual, ord=np.inf)),
        "mesh_contract_sha256": operators["assembly"]["mesh_contract"]["contract_sha256"],
        "material_contract_sha256": operators["material_contract"]["contract_sha256"],
        "operator_sha256": operators["operator_sha256"],
        "timing_s": {
            "prepare_and_assemble": assembly_s,
            "factorize_and_solve": solve_s,
            "export": export_s,
            "total": assembly_s + solve_s + export_s,
        },
        "exports": exports,
    }
    return result


def _time_axis(value: Any) -> np.ndarray:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) < 2:
        raise ValueError("time_s must contain at least two samples")
    times = np.asarray([_finite(item, "time_s") for item in value], dtype=float)
    if times[0] < 0.0 or np.any(np.diff(times) <= 0.0):
        raise ValueError("time_s must be nonnegative and strictly increasing")
    if times.size > 2000:
        raise ValueError("time_s is limited to 2000 samples")
    return times


def solve_vol2d_transient(request: Mapping[str, Any]) -> dict[str, Any]:
    """Advance ``M da/dt + K a = S i(t)`` with a theta method."""

    operators = assemble_vol2d_dynamics(request)
    stiffness = np.asarray(operators["assembly"]["field_matrix"], dtype=float)
    mass = np.asarray(operators["conductivity_mass_matrix"], dtype=float)
    sources = np.asarray(operators["assembly"]["source_matrix"], dtype=float)
    times = _time_axis(request.get("time_s"))
    currents = np.asarray(request.get("branch_current_history_a"), dtype=float)
    if currents.shape != (times.size, sources.shape[1]) or not np.all(np.isfinite(currents)):
        raise ValueError(
            "branch_current_history_a must have shape [len(time_s), branch_count]"
        )
    theta = _finite(request.get("theta", 1.0), "theta")
    if theta < 0.5 or theta > 1.0:
        raise ValueError("theta must be in [0.5, 1]")
    initial = request.get("initial_state", "zero")
    if initial == "zero":
        state = np.zeros(stiffness.shape[0], dtype=float)
        initial_policy = "zero"
    elif initial == "magnetostatic_equilibrium":
        state = np.linalg.solve(stiffness, sources @ currents[0])
        initial_policy = "magnetostatic_equilibrium"
    elif isinstance(initial, Sequence) and not isinstance(initial, (str, bytes)):
        state = np.asarray(initial, dtype=float)
        if state.shape != (stiffness.shape[0],) or not np.all(np.isfinite(state)):
            raise ValueError("initial_state vector must match free field DOFs")
        initial_policy = "provided"
    else:
        raise ValueError("initial_state must be zero, magnetostatic_equilibrium, or a vector")

    states = [state.copy()]
    energies = [float(0.5 * state @ stiffness @ state)]
    residuals = [0.0]
    dissipation = [0.0]
    for index, dt in enumerate(np.diff(times), start=1):
        previous = state
        lhs = mass / dt + theta * stiffness
        weighted_current = theta * currents[index] + (1.0 - theta) * currents[index - 1]
        rhs = (
            (mass / dt - (1.0 - theta) * stiffness) @ previous
            + sources @ weighted_current
        )
        try:
            state = np.linalg.solve(lhs, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f"transient step {index} matrix is singular") from exc
        derivative = (state - previous) / dt
        weighted_state = theta * state + (1.0 - theta) * previous
        residual = mass @ derivative + stiffness @ weighted_state - sources @ weighted_current
        states.append(state.copy())
        energies.append(float(0.5 * state @ stiffness @ state))
        residuals.append(float(np.linalg.norm(residual, ord=np.inf)))
        dissipation.append(float(derivative @ mass @ derivative * dt))

    unforced = bool(np.max(np.abs(currents)) == 0.0)
    passive_decay = all(
        right <= left + 1.0e-10 * max(1.0, abs(left))
        for left, right in zip(energies, energies[1:])
    ) if unforced else None
    waveform_sha = _sha256_json(
        {"time_s": times.tolist(), "branch_current_history_a": currents.tolist()}
    )
    return {
        "schema": DYNAMIC_SCHEMA,
        "status": "solved",
        "operation": "transient",
        "method": "theta",
        "theta": theta,
        "initial_state_policy": initial_policy,
        "time_s": times.tolist(),
        "branch_current_history_a": currents.tolist(),
        "field_state_history": [value.tolist() for value in states],
        "magnetic_energy_history_j": energies,
        "step_dissipation_j": dissipation,
        "step_residual_inf": residuals,
        "maximum_step_residual_inf": max(residuals),
        "unforced": unforced,
        "passive_energy_decay": passive_decay,
        "waveform_sha256": waveform_sha,
        "mesh_contract_sha256": operators["assembly"]["mesh_contract"]["contract_sha256"],
        "material_contract_sha256": operators["material_contract"]["contract_sha256"],
        "operator_sha256": operators["operator_sha256"],
    }


def compile_vol2d_state_space(request: Mapping[str, Any]) -> dict[str, Any]:
    """Compile the v3 field/source operators into the native MEX contract."""

    operators = assemble_vol2d_dynamics(request)
    circuit = request.get("state_space")
    if not isinstance(circuit, Mapping):
        raise ValueError("state_space must be an object")
    payload = dict(circuit)
    payload.update(
        {
            "field_matrix": operators["assembly"]["field_matrix"],
            "source_matrix": operators["assembly"]["source_matrix"],
            "field_rhs": operators["assembly"]["field_rhs"],
        }
    )
    result = compile_circuit_state_space(payload)
    result["artifact_identity"] = {
        "mesh_contract_sha256": operators["assembly"]["mesh_contract"]["contract_sha256"],
        "material_contract_sha256": operators["material_contract"]["contract_sha256"],
        "operator_sha256": operators["operator_sha256"],
    }
    return result


def analyze_vol2d_dynamics(request: Mapping[str, Any]) -> dict[str, Any]:
    """Dispatch a closed-world dimension-2 nonlinear/transient operation."""

    if not isinstance(request, Mapping):
        raise ValueError("request must be an object")
    operation = str(request.get("operation", ""))
    if operation not in _OPERATIONS:
        raise ValueError(f"operation must be one of {sorted(_OPERATIONS)}")
    if operation == "assemble":
        result = assemble_vol2d_dynamics(request)
    elif operation == "nonlinear_static":
        result = solve_vol2d_nonlinear_static(request)
    elif operation == "harmonic":
        result = solve_vol2d_harmonic(request)
    elif operation == "transient":
        result = solve_vol2d_transient(request)
    else:
        result = compile_vol2d_state_space(request)
    return {
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_version": {
            "radia_mcp": _package_version(),
            "ngsolve": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "schema": DYNAMIC_SCHEMA,
        "status": result["status"],
        "operation": operation,
        "result": result,
    }
