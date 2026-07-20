from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.motor_map_induction_identity_v56 import INDUCTION, MAP, validate_public_identity

CASE_IDS = {
    "v56_public_motormap_speed_torque_input_output_loss_efficiency_owner_mismatch",
    "v56_public_inductionmotor_slip_synchronousspeed_rotorfrequency_torque_owner_mismatch",
}


def _identity() -> dict[str, object]:
    generation = "jmag-public-v56-test"; generations = lambda fields: {field: generation for field in fields}
    speed = 3000.0; torque = 20.0; output = torque * speed * 2.0 * math.pi / 60.0
    losses = {"copper_w": 240.0, "iron_w": 110.0, "mechanical_w": 50.0}; input_power = output + sum(losses.values())
    frequency = 50.0; poles = 4; synchronous = 120.0 * frequency / poles; rotor = 1440.0; slip = (synchronous - rotor) / synchronous
    return {
        MAP: {"generation": generation, **generations(("speed_generation", "torque_generation", "input_generation", "output_generation", "loss_generation", "efficiency_generation", "owner_generation", "result_generation")), "speed_rpm": speed, "result_speed_rpm": speed, "torque_nm": torque, "result_torque_nm": torque, "input_power_w": input_power, "result_input_power_w": input_power, "output_power_w": output, "result_output_power_w": output, "loss_components_w": losses, "result_loss_components_w": losses, "efficiency": output / input_power, "result_efficiency": output / input_power, "result_owner": "result:map-v56", "accepted_result_owner": "result:map-v56", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64},
        INDUCTION: {"generation": generation, **generations(("frequency_generation", "pole_generation", "speed_generation", "slip_generation", "rotorfrequency_generation", "torque_generation", "owner_generation", "result_generation")), "supply_frequency_hz": frequency, "result_supply_frequency_hz": frequency, "pole_count": poles, "result_pole_count": poles, "synchronous_speed_rpm": synchronous, "result_synchronous_speed_rpm": synchronous, "rotor_speed_rpm": rotor, "result_rotor_speed_rpm": rotor, "slip": slip, "result_slip": slip, "rotor_frequency_hz": slip * frequency, "result_rotor_frequency_hz": slip * frequency, "torque_nm": 48.0, "result_torque_nm": 48.0, "torque_state": "motoring", "result_torque_state": "motoring", "result_owner": "result:induction-v56", "accepted_result_owner": "result:induction-v56", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
    }


def test_v56_positive_identity_is_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v56_frozen_result_mutations_are_rejected() -> None:
    identity = deepcopy(_identity()); identity[MAP]["result_efficiency"] = 2.0; identity[INDUCTION]["result_slip"] = -0.5
    assert not all(validate_public_identity(identity).values())


def test_v56_self_consistent_power_and_slip_errors_are_rejected() -> None:
    identity = deepcopy(_identity()); identity[MAP]["output_power_w"] = identity[MAP]["result_output_power_w"] = 9999.0; identity[INDUCTION]["slip"] = identity[INDUCTION]["result_slip"] = -0.5
    assert not all(validate_public_identity(identity).values())


def test_v56_malformed_losses_reject_without_raising() -> None:
    identity = deepcopy(_identity()); identity[MAP]["loss_components_w"] = []
    assert not all(validate_public_identity(identity).values())


def test_v56_numeric_digests_are_rejected() -> None:
    identity = deepcopy(_identity())
    numeric_digest = int("1" * 64)
    for contract_name in (MAP, INDUCTION):
        identity[contract_name]["result_sha256"] = numeric_digest
        identity[contract_name]["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_identity(identity).values())
