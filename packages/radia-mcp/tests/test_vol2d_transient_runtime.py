import copy
import hashlib
import json
import math

import pytest

from radia_mcp.fem.server import fem_vol2d_transient_runtime
from radia_mcp.radia_ngsolve.vol2d_transient_runtime import execute_transient_runtime


def _model(*, thermal: bool = False) -> dict:
    model = {
        "field": {
            "mass_matrix": [[1.0]],
            "stiffness_matrix": [[2.0]],
            "input_matrix": [[1.0]],
            "state_unit": "Wb/m",
        },
        "input_units": ["A"],
    }
    if thermal:
        model["thermal"] = {
            "capacity_matrix": [[4.0]],
            "conductance_matrix": [[1.0]],
            "loss_distribution": [1.0],
            "input_matrix": [[0.0]],
            "state_unit": "K",
        }
    return model


def _initialize(*, thermal: bool = False, field: float = 0.0) -> dict:
    return execute_transient_runtime(
        {
            "operation": "initialize",
            "model": _model(thermal=thermal),
            "initial_state": {
                "field": [field],
                "thermal": [0.0] if thermal else [],
                "previous_input": [0.0],
            },
        }
    )


def _step(initialized: dict, *, thermal: bool = False, current: float = 0.0) -> dict:
    return execute_transient_runtime(
        {
            "operation": "step",
            "model": _model(thermal=thermal),
            "state_token": initialized["state_token"],
            "expected_generation": initialized["state_token"]["generation"],
            "dt_s": 0.1,
            "theta": 1.0,
            "input": [current],
            "control": {"timeout_s": 1.0, "cancel_requested": False},
        }
    )


def _reseal(token: dict) -> None:
    payload = {key: value for key, value in token.items() if key != "state_sha256"}
    token["state_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def test_rl_decay_matches_backward_euler_and_is_passive() -> None:
    result = _step(_initialize(field=1.0))
    state = result["state_token"]["field_state"][0]
    assert state == pytest.approx(1.0 / 1.2)
    assert result["observables"]["field_quadratic_energy_j"] < 1.0
    assert result["observables"]["field_residual_inf"] < 1.0e-14


def test_forced_step_is_fixed_width_and_digest_bound() -> None:
    initialized = _initialize()
    result = _step(initialized, current=3.0)
    assert result["state_token"]["field_state"] == pytest.approx([0.25])
    assert result["state_token"]["generation"] == 1
    assert result["control"]["owned_processes_started"] == 0
    assert len(result["state_token"]["state_sha256"]) == 64


def test_field_loss_drives_optional_thermal_state() -> None:
    result = _step(_initialize(thermal=True), thermal=True, current=3.0)
    assert result["observables"]["joule_power_w"] == pytest.approx(6.25)
    assert result["state_token"]["thermal_state"][0] == pytest.approx(6.25 / 41.0)
    assert result["observables"]["thermal_residual_inf"] < 1.0e-14


def test_reset_repeats_and_terminate_closes_lifecycle() -> None:
    initialized = _initialize(field=1.0)
    first = _step(initialized)
    reset = execute_transient_runtime(
        {
            "operation": "reset",
            "model": _model(),
            "state_token": first["state_token"],
            "expected_generation": 1,
        }
    )
    second = _step(reset)
    assert second["state_token"]["field_state"] == pytest.approx(
        first["state_token"]["field_state"]
    )
    terminated = execute_transient_runtime(
        {
            "operation": "terminate",
            "model": _model(),
            "state_token": second["state_token"],
            "expected_generation": 3,
        }
    )
    assert terminated["status"] == "terminated"
    with pytest.raises(ValueError, match="not active"):
        execute_transient_runtime(
            {
                "operation": "step",
                "model": _model(),
                "state_token": terminated["state_token"],
                "expected_generation": 4,
                "dt_s": 0.1,
                "input": [0.0],
            }
        )


def test_warm_start_requires_exact_model_and_generation() -> None:
    first = _step(_initialize(), current=3.0)
    warm = execute_transient_runtime(
        {
            "operation": "initialize",
            "model": _model(),
            "warm_start": first["state_token"],
            "expected_generation": 1,
        }
    )
    assert warm["state_token"]["field_state"] == first["state_token"]["field_state"]
    assert warm["state_token"]["warm_start_source_sha256"] == first["state_token"]["state_sha256"]
    stale = copy.deepcopy(first["state_token"])
    stale["generation"] = 99
    with pytest.raises(ValueError, match="digest mismatch"):
        execute_transient_runtime(
            {
                "operation": "initialize",
                "model": _model(),
                "warm_start": stale,
                "expected_generation": 99,
            }
        )


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        (
            lambda request: (
                request["state_token"].update({"field_state": [0.0, 1.0]}),
                _reseal(request["state_token"]),
            ),
            "must contain 1",
        ),
        (lambda request: request.update({"expected_generation": 4}), "expected_generation"),
        (lambda request: request["control"].update({"cancel_requested": True}), "cancelled"),
    ],
)
def test_step_rejects_width_generation_and_cancel_mutations(mutation, match: str) -> None:
    initialized = _initialize()
    request = {
        "operation": "step",
        "model": _model(),
        "state_token": copy.deepcopy(initialized["state_token"]),
        "expected_generation": 0,
        "dt_s": 0.1,
        "theta": 1.0,
        "input": [0.0],
        "control": {"timeout_s": 1.0, "cancel_requested": False},
    }
    mutation(request)
    with pytest.raises(ValueError, match=match):
        execute_transient_runtime(request)


def test_rejects_stale_operator_identity() -> None:
    initialized = _initialize()
    changed = _model()
    changed["field"]["stiffness_matrix"] = [[3.0]]
    with pytest.raises(ValueError, match="different model"):
        execute_transient_runtime(
            {
                "operation": "step",
                "model": changed,
                "state_token": initialized["state_token"],
                "expected_generation": 0,
                "dt_s": 0.1,
                "input": [0.0],
            }
        )


def test_mcp_surface_returns_structured_invalid_input() -> None:
    result = json.loads(
        fem_vol2d_transient_runtime(
            json.dumps(
                {
                    "operation": "step",
                    "model": _model(),
                    "state_token": {},
                    "expected_generation": 0,
                    "dt_s": 0.1,
                    "input": [0.0],
                }
            )
        )
    )
    assert result["status"] == "invalid_input"
    assert result["pass"] is False
