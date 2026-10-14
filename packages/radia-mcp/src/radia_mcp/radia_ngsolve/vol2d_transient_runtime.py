"""Deterministic fixed-width transient runtime for portable 2-D operators.

The runtime is deliberately handle-free.  Every mutating call consumes a
digest-bound state token and returns the next generation, so a MATLAB MEX,
Simulink S-function, Python worker, or an MCP caller can own the token without
hidden server state.  The field lane advances ``M da/dt + K a = B u`` and an
optional thermal lane receives the resulting conductivity loss.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from typing import Any

import numpy as np


MODEL_SCHEMA = "radia.vol2d-transient-model.v1"
STATE_SCHEMA = "radia.vol2d-transient-state.v1"
RESULT_SCHEMA = "radia.vol2d-transient-runtime.v1"
_OPERATIONS = {"initialize", "step", "reset", "terminate"}


def _sha(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


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


def _matrix(value: Any, label: str) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite matrix") from exc
    if result.ndim != 2 or not result.size or not np.all(np.isfinite(result)):
        raise ValueError(f"{label} must be a finite non-empty matrix")
    return result


def _vector(value: Any, size: int, label: str) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must contain {size} finite values") from exc
    if result.shape != (size,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{label} must contain {size} finite values")
    return result


def _symmetric_psd(value: np.ndarray, label: str) -> np.ndarray:
    if value.shape[0] != value.shape[1]:
        raise ValueError(f"{label} must be square")
    scale = max(1.0, float(np.max(np.abs(value))))
    if float(np.max(np.abs(value - value.T))) > 1.0e-12 * scale:
        raise ValueError(f"{label} must be symmetric")
    result = 0.5 * (value + value.T)
    if float(np.linalg.eigvalsh(result)[0]) < -1.0e-12 * scale:
        raise ValueError(f"{label} must be positive semidefinite")
    return result


def _unit_list(raw: Any, size: int, label: str) -> list[str]:
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or len(raw) != size:
        raise ValueError(f"{label} must contain one unit per channel")
    units = [str(item).strip() for item in raw]
    if any(not item for item in units):
        raise ValueError(f"{label} entries must be non-empty")
    return units


def normalize_transient_model(raw: Any) -> dict[str, Any]:
    """Validate and digest a fixed-width field plus optional thermal model."""

    if not isinstance(raw, Mapping):
        raise ValueError("model must be an object")
    field = raw.get("field")
    if not isinstance(field, Mapping):
        raise ValueError("model.field must be an object")
    mass = _symmetric_psd(_matrix(field.get("mass_matrix"), "field.mass_matrix"), "field.mass_matrix")
    stiffness = _symmetric_psd(
        _matrix(field.get("stiffness_matrix"), "field.stiffness_matrix"),
        "field.stiffness_matrix",
    )
    if mass.shape != stiffness.shape:
        raise ValueError("field mass and stiffness matrices must have equal shape")
    state_width = mass.shape[0]
    if state_width > 512:
        raise ValueError("field state width is limited to 512")
    inputs = _matrix(field.get("input_matrix"), "field.input_matrix")
    if inputs.shape[0] != state_width:
        raise ValueError("field.input_matrix row count must match field state width")
    input_width = inputs.shape[1]
    if input_width > 128:
        raise ValueError("input width is limited to 128")
    field_state_unit = str(field.get("state_unit", "")).strip()
    if not field_state_unit:
        raise ValueError("field.state_unit is required")
    input_units = _unit_list(raw.get("input_units"), input_width, "input_units")

    normalized: dict[str, Any] = {
        "schema": MODEL_SCHEMA,
        "field": {
            "mass_matrix": mass.tolist(),
            "stiffness_matrix": stiffness.tolist(),
            "input_matrix": inputs.tolist(),
            "state_unit": field_state_unit,
        },
        "input_units": input_units,
    }

    thermal = raw.get("thermal")
    if thermal is not None:
        if not isinstance(thermal, Mapping):
            raise ValueError("model.thermal must be an object")
        capacity = _symmetric_psd(
            _matrix(thermal.get("capacity_matrix"), "thermal.capacity_matrix"),
            "thermal.capacity_matrix",
        )
        conductance = _symmetric_psd(
            _matrix(thermal.get("conductance_matrix"), "thermal.conductance_matrix"),
            "thermal.conductance_matrix",
        )
        if capacity.shape != conductance.shape:
            raise ValueError("thermal capacity and conductance matrices must have equal shape")
        thermal_width = capacity.shape[0]
        distribution = _vector(
            thermal.get("loss_distribution"), thermal_width, "thermal.loss_distribution"
        )
        if np.any(distribution < 0.0) or not math.isclose(
            float(np.sum(distribution)), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12
        ):
            raise ValueError("thermal.loss_distribution must be nonnegative and sum to one")
        direct = thermal.get("input_matrix")
        if direct is None:
            direct_matrix = np.zeros((thermal_width, input_width), dtype=float)
        else:
            direct_matrix = _matrix(direct, "thermal.input_matrix")
            if direct_matrix.shape != (thermal_width, input_width):
                raise ValueError(
                    "thermal.input_matrix must have shape [thermal width, input width]"
                )
        thermal_unit = str(thermal.get("state_unit", "")).strip()
        if not thermal_unit:
            raise ValueError("thermal.state_unit is required")
        normalized["thermal"] = {
            "capacity_matrix": capacity.tolist(),
            "conductance_matrix": conductance.tolist(),
            "loss_distribution": distribution.tolist(),
            "input_matrix": direct_matrix.tolist(),
            "state_unit": thermal_unit,
        }

    normalized["fixed_width_abi"] = {
        "input_width": input_width,
        "field_state_width": state_width,
        "thermal_state_width": len(normalized.get("thermal", {}).get("loss_distribution", [])),
        "dynamic_paths": False,
        "per_step_python_required": False,
    }
    normalized["model_sha256"] = _sha(normalized)
    expected = raw.get("expected_model_sha256")
    if expected is not None and str(expected) != normalized["model_sha256"]:
        raise ValueError("expected_model_sha256 does not match normalized model")
    return normalized


def _state_payload(token: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in token.items() if key != "state_sha256"}


def _seal_state(token: dict[str, Any]) -> dict[str, Any]:
    token["state_sha256"] = _sha(_state_payload(token))
    return token


def _new_state(model: Mapping[str, Any], initial: Mapping[str, Any] | None = None) -> dict[str, Any]:
    raw = dict(initial or {})
    field_width = model["fixed_width_abi"]["field_state_width"]
    thermal_width = model["fixed_width_abi"]["thermal_state_width"]
    input_width = model["fixed_width_abi"]["input_width"]
    field = _vector(raw.get("field", [0.0] * field_width), field_width, "initial_state.field")
    thermal = _vector(
        raw.get("thermal", [0.0] * thermal_width), thermal_width, "initial_state.thermal"
    )
    previous_input = _vector(
        raw.get("previous_input", [0.0] * input_width), input_width, "initial_state.previous_input"
    )
    token = {
        "schema": STATE_SCHEMA,
        "status": "active",
        "model_sha256": model["model_sha256"],
        "generation": 0,
        "step_count": 0,
        "time_s": 0.0,
        "field_state": field.tolist(),
        "thermal_state": thermal.tolist(),
        "previous_input": previous_input.tolist(),
        "initial_state": {
            "field": field.tolist(),
            "thermal": thermal.tolist(),
            "previous_input": previous_input.tolist(),
        },
    }
    return _seal_state(token)


def _validate_state(
    raw: Any, model: Mapping[str, Any], expected_generation: Any, *, active: bool = True
) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ValueError("state_token must be an object")
    token = json.loads(json.dumps(raw))
    if token.get("schema") != STATE_SCHEMA:
        raise ValueError("state_token schema is unsupported")
    if token.get("model_sha256") != model["model_sha256"]:
        raise ValueError("state_token belongs to a different model")
    if token.get("state_sha256") != _sha(_state_payload(token)):
        raise ValueError("state_token digest mismatch")
    try:
        generation = int(token.get("generation"))
        expected = int(expected_generation)
    except (TypeError, ValueError) as exc:
        raise ValueError("expected_generation must be an integer") from exc
    if generation != expected:
        raise ValueError("expected_generation does not match state_token generation")
    if active and token.get("status") != "active":
        raise ValueError("state_token is not active")
    _vector(token.get("field_state"), model["fixed_width_abi"]["field_state_width"], "state_token.field_state")
    _vector(token.get("thermal_state"), model["fixed_width_abi"]["thermal_state_width"], "state_token.thermal_state")
    _vector(token.get("previous_input"), model["fixed_width_abi"]["input_width"], "state_token.previous_input")
    return token


def _theta_step(
    mass: np.ndarray,
    stiffness: np.ndarray,
    previous: np.ndarray,
    old_rhs: np.ndarray,
    new_rhs: np.ndarray,
    dt: float,
    theta: float,
    label: str,
) -> tuple[np.ndarray, np.ndarray, float]:
    lhs = mass / dt + theta * stiffness
    rhs = (mass / dt - (1.0 - theta) * stiffness) @ previous
    rhs += theta * new_rhs + (1.0 - theta) * old_rhs
    try:
        current = np.linalg.solve(lhs, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError(f"{label} effective matrix is singular") from exc
    derivative = (current - previous) / dt
    weighted = theta * current + (1.0 - theta) * previous
    residual = mass @ derivative + stiffness @ weighted
    residual -= theta * new_rhs + (1.0 - theta) * old_rhs
    return current, derivative, float(np.linalg.norm(residual, ord=np.inf))


def execute_transient_runtime(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Execute one deterministic lifecycle operation on a fixed-width token."""

    if not isinstance(packet, Mapping):
        raise ValueError("packet must be an object")
    operation = str(packet.get("operation", ""))
    if operation not in _OPERATIONS:
        raise ValueError(f"operation must be one of {sorted(_OPERATIONS)}")
    model = normalize_transient_model(packet.get("model"))
    if operation == "initialize":
        warm = packet.get("warm_start")
        if warm is not None:
            source = _validate_state(
                warm, model, packet.get("expected_generation"), active=True
            )
            initial = {
                "field": source["field_state"],
                "thermal": source["thermal_state"],
                "previous_input": source["previous_input"],
            }
            state = _new_state(model, initial)
            state["warm_start_source_sha256"] = source["state_sha256"]
            state = _seal_state(state)
        else:
            state = _new_state(model, packet.get("initial_state"))
        return {
            "schema": RESULT_SCHEMA,
            "status": "initialized",
            "operation": operation,
            "model": model,
            "state_token": state,
            "fixed_width_abi": model["fixed_width_abi"],
        }

    token = _validate_state(
        packet.get("state_token"), model, packet.get("expected_generation"), active=True
    )
    if operation == "reset":
        reset = _new_state(model, token["initial_state"])
        reset["generation"] = token["generation"] + 1
        reset = _seal_state(reset)
        return {
            "schema": RESULT_SCHEMA,
            "status": "reset",
            "operation": operation,
            "model_sha256": model["model_sha256"],
            "state_token": reset,
        }
    if operation == "terminate":
        token["status"] = "terminated"
        token["generation"] += 1
        token = _seal_state(token)
        return {
            "schema": RESULT_SCHEMA,
            "status": "terminated",
            "operation": operation,
            "model_sha256": model["model_sha256"],
            "state_token": token,
        }

    control = packet.get("control", {})
    if not isinstance(control, Mapping):
        raise ValueError("control must be an object")
    if control.get("cancel_requested") is True:
        raise ValueError("step cancelled before execution")
    timeout_s = _positive(control.get("timeout_s", 30.0), "control.timeout_s")
    if timeout_s > 3600.0:
        raise ValueError("control.timeout_s must not exceed 3600 seconds")
    deadline = control.get("deadline_utc")
    if deadline is not None:
        try:
            parsed = datetime.fromisoformat(str(deadline).replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError("control.deadline_utc must be ISO-8601") from exc
        if parsed.tzinfo is None or parsed <= datetime.now(timezone.utc):
            raise ValueError("step deadline has expired")

    dt = _positive(packet.get("dt_s"), "dt_s")
    theta = _finite(packet.get("theta", 1.0), "theta")
    if theta < 0.5 or theta > 1.0:
        raise ValueError("theta must be in [0.5, 1]")
    input_width = model["fixed_width_abi"]["input_width"]
    new_input = _vector(packet.get("input"), input_width, "input")
    old_input = _vector(token["previous_input"], input_width, "state_token.previous_input")

    field_model = model["field"]
    field_mass = np.asarray(field_model["mass_matrix"], dtype=float)
    field_stiffness = np.asarray(field_model["stiffness_matrix"], dtype=float)
    field_input = np.asarray(field_model["input_matrix"], dtype=float)
    old_field = np.asarray(token["field_state"], dtype=float)
    field, derivative, field_residual = _theta_step(
        field_mass,
        field_stiffness,
        old_field,
        field_input @ old_input,
        field_input @ new_input,
        dt,
        theta,
        "field",
    )
    joule_power = float(derivative @ field_mass @ derivative)

    thermal = np.asarray(token["thermal_state"], dtype=float)
    thermal_residual = 0.0
    if "thermal" in model:
        thermal_model = model["thermal"]
        capacity = np.asarray(thermal_model["capacity_matrix"], dtype=float)
        conductance = np.asarray(thermal_model["conductance_matrix"], dtype=float)
        distribution = np.asarray(thermal_model["loss_distribution"], dtype=float)
        direct = np.asarray(thermal_model["input_matrix"], dtype=float)
        old_rhs = direct @ old_input
        new_rhs = direct @ new_input + distribution * joule_power
        thermal, _, thermal_residual = _theta_step(
            capacity,
            conductance,
            thermal,
            old_rhs,
            new_rhs,
            dt,
            theta,
            "thermal",
        )

    token.update(
        {
            "generation": token["generation"] + 1,
            "step_count": token["step_count"] + 1,
            "time_s": token["time_s"] + dt,
            "field_state": field.tolist(),
            "thermal_state": thermal.tolist(),
            "previous_input": new_input.tolist(),
        }
    )
    token = _seal_state(token)
    return {
        "schema": RESULT_SCHEMA,
        "status": "stepped",
        "operation": operation,
        "model_sha256": model["model_sha256"],
        "state_token": token,
        "observables": {
            "joule_power_w": joule_power,
            "field_quadratic_energy_j": float(0.5 * field @ field_stiffness @ field),
            "field_residual_inf": field_residual,
            "thermal_residual_inf": thermal_residual,
        },
        "control": {
            "cancel_checked_before_solve": True,
            "timeout_s": timeout_s,
            "owned_processes_started": 0,
        },
    }
