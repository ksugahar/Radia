from __future__ import annotations

import copy
import json

from radia_mcp.motor.server import motor_electrothermal_result_chain_gate
from radia_mcp.motor.thermal_handoff import evaluate_motor_electrothermal_result_chain


def _chain() -> dict:
    ids = {
        "electromagnetic": "em-result",
        "stator_core_loss": "stator-loss-result",
        "rotor_core_loss": "rotor-loss-result",
        "thermal": "thermal-result",
    }
    digests = {
        "electromagnetic": "a" * 64,
        "stator_core_loss": "b" * 64,
        "rotor_core_loss": "c" * 64,
        "thermal": "d" * 64,
    }
    stages = [
        {
            "stage": "electromagnetic",
            "artifact_id": ids["electromagnetic"],
            "result_digest": digests["electromagnetic"],
            "input_result_digests": {},
            "completed": True,
            "fresh": True,
            "solve_s": 12.0,
        },
        {
            "stage": "stator_core_loss",
            "artifact_id": ids["stator_core_loss"],
            "result_digest": digests["stator_core_loss"],
            "input_result_digests": {
                ids["electromagnetic"]: digests["electromagnetic"]
            },
            "completed": True,
            "fresh": True,
            "solve_s": 2.0,
        },
        {
            "stage": "rotor_core_loss",
            "artifact_id": ids["rotor_core_loss"],
            "result_digest": digests["rotor_core_loss"],
            "input_result_digests": {
                ids["electromagnetic"]: digests["electromagnetic"]
            },
            "completed": True,
            "fresh": True,
            "solve_s": 2.5,
        },
        {
            "stage": "thermal",
            "artifact_id": ids["thermal"],
            "result_digest": digests["thermal"],
            "input_result_digests": {
                ids["electromagnetic"]: digests["electromagnetic"],
                ids["stator_core_loss"]: digests["stator_core_loss"],
                ids["rotor_core_loss"]: digests["rotor_core_loss"],
            },
            "completed": True,
            "fresh": True,
            "solve_s": 3.0,
        },
    ]
    powers = {
        "rotor_conductor_joule": 20.0,
        "phase_u_joule": 10.0,
        "phase_v_joule": 10.0,
        "phase_w_joule": 10.0,
        "stator_core_loss": 4.0,
        "rotor_core_loss": 2.0,
    }
    owners = {
        "rotor_conductor_joule": "electromagnetic",
        "phase_u_joule": "electromagnetic",
        "phase_v_joule": "electromagnetic",
        "phase_w_joule": "electromagnetic",
        "stator_core_loss": "stator_core_loss",
        "rotor_core_loss": "rotor_core_loss",
    }
    return {
        "schema": "motor-electrothermal-result-chain/v1",
        "stages": stages,
        "source_buckets": [
            {
                "channel": channel,
                "upstream_stage": owner,
                "upstream_artifact_id": ids[owner],
                "upstream_result_digest": digests[owner],
                "power_W": power,
            }
            for channel, power in powers.items()
            for owner in [owners[channel]]
        ],
        "symmetry_fraction": 0.5,
        "thermal_summary": {
            "steady_state": True,
            "input_power_W": 28.0,
            "ambient_temperature_C": 20.0,
            "maximum_temperature_C": 82.0,
        },
    }


def test_electrothermal_result_chain_accepts_digest_pinned_symmetry_handoff() -> None:
    result = evaluate_motor_electrothermal_result_chain(_chain())
    assert result["status"] == "ok"
    assert result["metrics"]["expected_thermal_input_W"] == 28.0
    assert result["metrics"]["temperature_rise_C"] == 62.0
    assert all(result["checks"].values())

    wrapped = json.loads(motor_electrothermal_result_chain_gate(json.dumps(_chain())))
    assert wrapped["status"] == "ok"
    assert wrapped["policy"] == "fresh_digest_pinned_symmetry_scaled_power_handoff"


def test_electrothermal_result_chain_rejects_stale_upstream_digest() -> None:
    chain = _chain()
    chain["stages"][3]["input_result_digests"]["em-result"] = "e" * 64
    result = evaluate_motor_electrothermal_result_chain(chain)
    assert result["status"] == "needs_attention"
    assert result["checks"]["upstream_result_digests_pinned"] is False


def test_electrothermal_result_chain_rejects_missing_channel_and_wrong_symmetry() -> None:
    chain = _chain()
    chain["source_buckets"].pop()
    chain["thermal_summary"]["input_power_W"] = 56.0
    result = evaluate_motor_electrothermal_result_chain(chain)
    assert result["status"] == "needs_attention"
    assert result["checks"]["six_loss_channels_owned_once"] is False
    assert result["checks"]["symmetry_scaled_power_closure"] is False


def test_electrothermal_result_chain_rejects_reused_result_and_no_temperature_rise() -> None:
    chain = copy.deepcopy(_chain())
    chain["stages"][2]["result_digest"] = chain["stages"][1]["result_digest"]
    chain["thermal_summary"]["maximum_temperature_C"] = 20.0
    result = evaluate_motor_electrothermal_result_chain(chain)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fresh_unique_stage_results"] is False
    assert result["checks"]["steady_temperature_rise_present"] is False


def test_electrothermal_result_chain_rejects_invalid_json() -> None:
    try:
        motor_electrothermal_result_chain_gate("not-json")
    except ValueError as exc:
        assert "chain_json must be valid JSON" in str(exc)
    else:
        raise AssertionError("invalid JSON was accepted")
