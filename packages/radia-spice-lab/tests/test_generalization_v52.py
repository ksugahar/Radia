from copy import deepcopy

from ltspice_converter.ltspice_v52_gates import MOS, TLINE, validate_ltspice_v52_identity


PROMOTED_CASE_IDS = {
    "v52_public_mosfet_operatingpoint_region_gm_gds_capacitance_temperature_owner_mismatch",
    "v52_public_transmissionline_delay_impedance_termination_reflection_event_owner_mismatch",
}


def _positive() -> dict[str, object]:
    generation = "ltspice-public-v52"
    generations = lambda names: {name: generation for name in names}
    gamma = 0.2
    return {
        MOS: {
            "generation_id": generation, **generations(("bias_generation_id", "region_generation_id", "small_signal_generation_id", "capacitance_generation_id", "temperature_generation_id", "owner_generation_id", "result_generation_id")),
            "vgs_v": 5.0, "result_vgs_v": 5.0, "vds_v": 12.0, "result_vds_v": 12.0,
            "vth_v": 2.5, "result_vth_v": 2.5, "id_a": 3.2, "result_id_a": 3.2,
            "region": "saturation", "result_region": "saturation", "gm_s": 2.1, "result_gm_s": 2.1,
            "gds_s": 0.025, "result_gds_s": 0.025,
            "capacitances": {"cgs_f": 1e-9, "cgd_f": 2e-10, "cds_f": 1e-10},
            "result_capacitances": {"cgs_f": 1e-9, "cgd_f": 2e-10, "cds_f": 1e-10},
            "model_temperature_c": 75.0, "result_model_temperature_c": 75.0,
            "device_owner": "device:M1", "result_device_owner": "device:M1",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        TLINE: {
            "generation_id": generation, **generations(("delay_generation_id", "impedance_generation_id", "termination_generation_id", "reflection_generation_id", "event_generation_id", "owner_generation_id", "result_generation_id")),
            "propagation_delay_s": 8e-9, "result_propagation_delay_s": 8e-9,
            "characteristic_impedance_ohm": 50.0, "result_characteristic_impedance_ohm": 50.0,
            "termination_ohm": 75.0, "result_termination_ohm": 75.0,
            "reflection_coefficient": gamma, "result_reflection_coefficient": gamma,
            "reflection_event": {"launch_time_s": 2e-9, "reflection_time_s": 18e-9, "incident_v": 1.0, "reflected_v": gamma},
            "result_reflection_event": {"launch_time_s": 2e-9, "reflection_time_s": 18e-9, "incident_v": 1.0, "reflected_v": gamma},
            "waveform_owner": "waveform:tline-v52", "result_waveform_owner": "waveform:tline-v52",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v52_positive_identity_is_accepted() -> None:
    assert validate_ltspice_v52_identity(_positive()) is True


def test_v52_frozen_public_counterfactuals_are_rejected() -> None:
    value = deepcopy(_positive())
    value[MOS].update({"result_region": "linear", "result_gm_s": 0.1, "result_model_temperature_c": 25.0, "result_device_owner": "device:foreign"})
    value[TLINE].update({"result_propagation_delay_s": 4e-9, "result_reflection_coefficient": -0.2, "result_reflection_event": {"launch_time_s": 0.0, "reflection_time_s": 4e-9, "incident_v": 1.0, "reflected_v": -0.2}, "result_waveform_owner": "waveform:foreign"})
    assert validate_ltspice_v52_identity(value) is False


def test_v52_self_consistent_wrong_region_and_reflection_are_rejected() -> None:
    value = deepcopy(_positive())
    value[MOS]["region"] = value[MOS]["result_region"] = "linear"
    value[TLINE]["reflection_coefficient"] = value[TLINE]["result_reflection_coefficient"] = -0.2
    value[TLINE]["reflection_event"] = value[TLINE]["result_reflection_event"] = {"launch_time_s": 2e-9, "reflection_time_s": 18e-9, "incident_v": 1.0, "reflected_v": -0.2}
    assert validate_ltspice_v52_identity(value) is False
