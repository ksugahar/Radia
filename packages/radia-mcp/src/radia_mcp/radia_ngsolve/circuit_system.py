"""Executable linear field-circuit and AGE sweep kernels.

The finite-element matrix and coil source vectors are supplied explicitly so
the same algebra serves planar 2D and axisymmetric formulations.  Parallel
circuits are solved as one augmented system; prescribed total current is never
split before the field equations are assembled.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from .airgap_element import annular_rotation_phase, planar_translation_phase


SCHEMA = "radia.circuit-field-analysis.v1"


def _complex_scalar(value: Any, label: str) -> complex:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a finite scalar")
    if isinstance(value, (int, float, complex)):
        result = complex(value)
    elif isinstance(value, Mapping) and set(value) <= {"re", "im"}:
        result = complex(value.get("re", 0.0), value.get("im", 0.0))
    else:
        raise ValueError(f"{label} must be a number or {{re, im}} object")
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError(f"{label} must be finite")
    return result


def _vector(value: Any, label: str) -> np.ndarray:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{label} must be a sequence")
    result = np.asarray(
        [_complex_scalar(item, f"{label}[{index}]") for index, item in enumerate(value)],
        dtype=complex,
    )
    if result.ndim != 1 or result.size == 0:
        raise ValueError(f"{label} must be a non-empty vector")
    return result


def _matrix(value: Any, label: str) -> np.ndarray:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{label} must be a non-empty row sequence")
    rows: list[list[complex]] = []
    width: int | None = None
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            raise ValueError(f"{label}[{row_index}] must be a sequence")
        parsed = [
            _complex_scalar(item, f"{label}[{row_index}][{column_index}]")
            for column_index, item in enumerate(row)
        ]
        if width is None:
            width = len(parsed)
        if not parsed or len(parsed) != width:
            raise ValueError(f"{label} rows must have one consistent non-zero width")
        rows.append(parsed)
    return np.asarray(rows, dtype=complex)


def _pair(value: complex) -> list[float]:
    return [float(value.real), float(value.imag)]


def _pairs(value: np.ndarray) -> list[Any]:
    array = np.asarray(value, dtype=complex)
    if array.ndim == 1:
        return [_pair(item) for item in array]
    return [[_pair(item) for item in row] for row in array]


def _field_inputs(payload: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    field = _matrix(payload.get("field_matrix"), "field_matrix")
    sources = _matrix(payload.get("source_matrix"), "source_matrix")
    rhs = _vector(payload.get("field_rhs", [0.0] * field.shape[0]), "field_rhs")
    if field.shape[0] != field.shape[1]:
        raise ValueError("field_matrix must be square")
    if sources.shape[0] != field.shape[0]:
        raise ValueError("source_matrix row count must match field_matrix")
    if rhs.size != field.shape[0]:
        raise ValueError("field_rhs length must match field_matrix")
    return field, sources, rhs


def _frequency(payload: Mapping[str, Any]) -> tuple[float, complex, str]:
    frequency = float(payload.get("frequency_hz", 0.0))
    if not math.isfinite(frequency) or frequency < 0.0:
        raise ValueError("frequency_hz must be finite and nonnegative")
    convention = str(payload.get("phase_convention", "exp(+jwt)"))
    if convention not in {"exp(+jwt)", "exp(-jwt)"}:
        raise ValueError("phase_convention must be exp(+jwt) or exp(-jwt)")
    sign = 1.0 if convention == "exp(+jwt)" else -1.0
    return frequency, sign * 1j * 2.0 * math.pi * frequency, convention


def solve_parallel_field_circuit(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Solve a current-driven parallel field-circuit augmented system."""

    field, sources, rhs = _field_inputs(payload)
    impedance = _vector(payload.get("branch_impedance_ohm"), "branch_impedance_ohm")
    if impedance.size != sources.shape[1]:
        raise ValueError("one branch impedance is required per source_matrix column")
    if np.any(np.real(impedance) < 0.0) or np.any(np.abs(impedance) == 0.0):
        raise ValueError("branch impedances must be non-zero with nonnegative resistance")
    total_current = _complex_scalar(payload.get("total_current_a", 0.0), "total_current_a")
    frequency, derivative, convention = _frequency(payload)

    nfield, nbranch = sources.shape
    size = nfield + nbranch + 1
    augmented = np.zeros((size, size), dtype=complex)
    augmented[:nfield, :nfield] = field
    augmented[:nfield, nfield : nfield + nbranch] = -sources
    augmented[nfield : nfield + nbranch, :nfield] = derivative * sources.T
    augmented[nfield : nfield + nbranch, nfield : nfield + nbranch] = np.diag(impedance)
    augmented[nfield : nfield + nbranch, -1] = -1.0
    augmented[-1, nfield : nfield + nbranch] = 1.0
    augmented_rhs = np.concatenate((rhs, np.zeros(nbranch, dtype=complex), [total_current]))
    try:
        solution = np.linalg.solve(augmented, augmented_rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("parallel field-circuit augmented system is singular") from exc

    field_state = solution[:nfield]
    branch_current = solution[nfield : nfield + nbranch]
    common_voltage = solution[-1]
    flux = sources.T @ field_state
    field_residual = field @ field_state - sources @ branch_current - rhs
    branch_residual = impedance * branch_current + derivative * flux - common_voltage
    total_residual = np.sum(branch_current) - total_current
    residual_inf = max(
        float(np.max(np.abs(field_residual))),
        float(np.max(np.abs(branch_residual))),
        float(abs(total_residual)),
    )
    return {
        "schema": SCHEMA,
        "status": "solved",
        "operation": "parallel",
        "frequency_hz": frequency,
        "phase_convention": convention,
        "field_state": _pairs(field_state),
        "branch_current_a": _pairs(branch_current),
        "common_terminal_voltage_v": _pair(common_voltage),
        "flux_linkage_wb_turn": _pairs(flux),
        "total_current_a": _pair(total_current),
        "equal_current_split_assumed": False,
        "residual": {
            "field_inf": float(np.max(np.abs(field_residual))),
            "branch_voltage_inf": float(np.max(np.abs(branch_residual))),
            "total_current_abs": float(abs(total_residual)),
            "maximum": residual_inf,
        },
        "augmented_shape": list(augmented.shape),
        "condition_number": float(np.linalg.cond(augmented)),
    }


def solve_series_field_circuit(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Solve the field for one prescribed current through every series branch."""

    field, sources, rhs = _field_inputs(payload)
    impedance = _vector(payload.get("branch_impedance_ohm"), "branch_impedance_ohm")
    if impedance.size != sources.shape[1]:
        raise ValueError("one branch impedance is required per source_matrix column")
    if np.any(np.real(impedance) < 0.0):
        raise ValueError("branch resistance must be nonnegative")
    current = _complex_scalar(payload.get("circuit_current_a", 0.0), "circuit_current_a")
    frequency, derivative, convention = _frequency(payload)
    branch_current = np.full(sources.shape[1], current, dtype=complex)
    try:
        field_state = np.linalg.solve(field, rhs + sources @ branch_current)
    except np.linalg.LinAlgError as exc:
        raise ValueError("series field system is singular") from exc
    flux = sources.T @ field_state
    branch_voltage = impedance * branch_current + derivative * flux
    field_residual = field @ field_state - sources @ branch_current - rhs
    return {
        "schema": SCHEMA,
        "status": "solved",
        "operation": "series",
        "frequency_hz": frequency,
        "phase_convention": convention,
        "field_state": _pairs(field_state),
        "branch_current_a": _pairs(branch_current),
        "branch_voltage_v": _pairs(branch_voltage),
        "circuit_terminal_voltage_v": _pair(np.sum(branch_voltage)),
        "flux_linkage_wb_turn": _pairs(flux),
        "residual": {"field_inf": float(np.max(np.abs(field_residual)))},
        "condition_number": float(np.linalg.cond(field)),
    }


def compile_age_sweep(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Compile an annular or planar AGE position sweep with one fixed operator."""

    kind = str(payload.get("kind", ""))
    positions = payload.get("positions", [])
    if not isinstance(positions, Sequence) or isinstance(positions, (str, bytes)) or not positions:
        raise ValueError("positions must be a non-empty sequence")
    position_values = [float(value) for value in positions]
    if not all(math.isfinite(value) for value in position_values):
        raise ValueError("positions must be finite")

    if kind == "annular_age":
        modes = payload.get("harmonics", [])
        if not isinstance(modes, Sequence) or isinstance(modes, (str, bytes)) or not modes:
            raise ValueError("harmonics must be a non-empty sequence")
        normalized = [int(value) for value in modes]
        if any(
            isinstance(value, bool) or int(value) != value or int(value) <= 0
            for value in modes
        ):
            raise ValueError("harmonics must be positive integers")
        factors = [
            {str(mode): _pair(annular_rotation_phase(mode, position)) for mode in normalized}
            for position in position_values
        ]
        observable = "mesh-independent harmonic torque"
        coordinate = "angle_rad"
    elif kind == "planar_age":
        modes = payload.get("wavenumbers_per_m", [])
        if not isinstance(modes, Sequence) or isinstance(modes, (str, bytes)) or not modes:
            raise ValueError("wavenumbers_per_m must be a non-empty sequence")
        normalized = [float(value) for value in modes]
        if any(not math.isfinite(value) or value <= 0.0 for value in normalized):
            raise ValueError("wavenumbers_per_m must be positive and finite")
        factors = [
            {
                format(mode, ".17g"): _pair(planar_translation_phase(mode, position))
                for mode in normalized
            }
            for position in position_values
        ]
        observable = "mesh-independent harmonic thrust"
        coordinate = "displacement_m"
    else:
        raise ValueError("kind must be annular_age or planar_age")

    closure = {
        key: float(abs(complex(*factors[-1][key]) - complex(*factors[0][key])))
        for key in factors[0]
    }
    return {
        "schema": SCHEMA,
        "status": "compiled",
        "operation": "age_sweep",
        "kind": kind,
        "coordinate": coordinate,
        "positions": position_values,
        "phase_factors": factors,
        "endpoint_closure_error": closure,
        "factorization_count": 1,
        "operator_rebuild_count": 0,
        "mesh_rebuild_count": 0,
        "position_solve_count": len(position_values),
        "mechanical_observable": observable,
    }


def compile_circuit_state_space(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Reduce a linear magnetic field matrix to an RL state-space MEX model."""

    field, sources, _ = _field_inputs(payload)
    resistance = _vector(payload.get("branch_resistance_ohm"), "branch_resistance_ohm")
    if resistance.size != sources.shape[1]:
        raise ValueError("one branch resistance is required per source_matrix column")
    if np.max(np.abs(field.imag)) > 1.0e-13 or np.max(np.abs(sources.imag)) > 1.0e-13:
        raise ValueError("time-domain state-space reduction requires real field/source matrices")
    if np.max(np.abs(resistance.imag)) > 1.0e-13 or np.any(resistance.real <= 0.0):
        raise ValueError("branch resistances must be real and positive")
    sample_time = float(payload.get("sample_time_s", 0.0))
    if not math.isfinite(sample_time) or sample_time <= 0.0:
        raise ValueError("sample_time_s must be positive and finite")
    mode = str(payload.get("voltage_input_mode", "common"))
    if mode not in {"common", "per_branch"}:
        raise ValueError("voltage_input_mode must be common or per_branch")

    real_field = field.real
    real_sources = sources.real
    try:
        inductance = real_sources.T @ np.linalg.solve(real_field, real_sources)
    except np.linalg.LinAlgError as exc:
        raise ValueError("field matrix is singular") from exc
    inductance = 0.5 * (inductance + inductance.T)
    eigenvalues = np.linalg.eigvalsh(inductance)
    if eigenvalues[0] <= max(1.0, eigenvalues[-1]) * 1.0e-12:
        raise ValueError("reduced inductance matrix must be positive definite")

    resistance_matrix = np.diag(resistance.real)
    input_map = np.ones((resistance.size, 1)) if mode == "common" else np.eye(resistance.size)
    continuous_a = -np.linalg.solve(inductance, resistance_matrix)
    continuous_b = np.linalg.solve(inductance, input_map)
    eigvals, eigvecs = np.linalg.eig(continuous_a)
    discrete_a = np.real_if_close(
        eigvecs @ np.diag(np.exp(eigvals * sample_time)) @ np.linalg.inv(eigvecs),
        tol=1000,
    ).real
    discrete_b = np.linalg.solve(
        continuous_a, (discrete_a - np.eye(resistance.size)) @ continuous_b
    )
    output_c = np.vstack((np.eye(resistance.size), inductance))
    output_d = np.zeros((2 * resistance.size, input_map.shape[1]))
    return {
        "schema": "radia.circuit-field.state-space.v1",
        "status": "compiled",
        "operation": "state_space",
        "backend": "native-mex-sfunction",
        "mex_s_function": "radia_state_space_mex_sfunction",
        "mex_commands": [
            "simulink.state_space.create",
            "simulink.state_space.info",
            "simulink.state_space.step",
            "simulink.state_space.reset",
            "simulink.state_space.destroy",
        ],
        "sample_time_s": sample_time,
        "voltage_input_mode": mode,
        "state_order": int(resistance.size),
        "input_count": int(input_map.shape[1]),
        "output_count": int(output_c.shape[0]),
        "output_order": [
            *[f"branch_current_a:{index}" for index in range(resistance.size)],
            *[f"flux_linkage_wb_turn:{index}" for index in range(resistance.size)],
        ],
        "inductance_matrix_h": inductance.tolist(),
        "A": continuous_a.tolist(),
        "B": continuous_b.tolist(),
        "C": output_c.tolist(),
        "D": output_d.tolist(),
        "Ad": discrete_a.tolist(),
        "Bd": discrete_b.tolist(),
        "Cd": output_c.tolist(),
        "Dd": output_d.tolist(),
        "x0": [0.0] * resistance.size,
        "stable": bool(np.all(np.real(eigvals) < 0.0)),
        "python_per_step": False,
        "field_factorization_runtime_count": 0,
    }


def analyze_circuit_field(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Dispatch one closed-world circuit/AGE analysis request."""

    if not isinstance(payload, Mapping):
        raise ValueError("payload must be an object")
    operation = str(payload.get("operation", ""))
    if operation == "parallel":
        return solve_parallel_field_circuit(payload)
    if operation == "series":
        return solve_series_field_circuit(payload)
    if operation == "age_sweep":
        return compile_age_sweep(payload)
    if operation == "state_space":
        return compile_circuit_state_space(payload)
    raise ValueError("operation must be parallel, series, age_sweep, or state_space")
