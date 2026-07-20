from copy import deepcopy

from radia_mcp.radia_ngsolve.wave_sar_identity_v54 import CUTOFF, SAR, validate_public_v54_identity


PROMOTED_CASE_IDS = {
    "v54_public_waveguide_cutoff_mode_normalization_power_impedance_port_owner_mismatch",
    "v54_public_sar_average_mass_density_voxel_frequency_field_owner_mismatch",
}
_C0 = 299792458.0
_ETA0 = 376.730313668


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _payload():
    width = 0.02286
    cutoff_hz = _C0 / (2.0 * width)
    frequency_hz = 10.0e9
    impedance = _ETA0 / (1.0 - (cutoff_hz / frequency_hz) ** 2) ** 0.5
    cutoff = {
        **_generations("cutoff-v54", ("cutoff_generation", "mode_generation", "normalization_generation", "power_generation", "impedance_generation", "owner_generation", "result_generation")),
        "mode_id": "TE10", "result_mode_id": "TE10",
        "waveguide_width_m": width, "result_waveguide_width_m": width,
        "cutoff_frequency_hz": cutoff_hz, "result_cutoff_frequency_hz": cutoff_hz,
        "sample_frequency_hz": frequency_hz, "result_sample_frequency_hz": frequency_hz,
        "field_normalization": "unit_modal_power_w", "result_field_normalization": "unit_modal_power_w",
        "modal_power_w": 1.0, "result_modal_power_w": 1.0,
        "mode_impedance_ohm": impedance, "result_mode_impedance_ohm": impedance,
        "port_owner": "port:v54", "result_port_owner": "port:v54",
        "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
    }
    voxels = [1.0e-6] * 10
    fields = [12.0 + value for value in range(10)]
    sar = {
        **_generations("sar-v54", ("mass_generation", "density_generation", "voxel_generation", "frequency_generation", "field_generation", "owner_generation", "result_generation")),
        "averaging_mass_kg": 0.01, "result_averaging_mass_kg": 0.01,
        "tissue_density_kg_m3": 1000.0, "result_tissue_density_kg_m3": 1000.0,
        "voxel_support_m3": voxels, "result_voxel_support_m3": voxels,
        "electric_field_rms_v_m": fields, "result_electric_field_rms_v_m": fields,
        "frequency_hz": 2.45e9, "result_frequency_hz": 2.45e9,
        "field_solution_sha256": "2" * 64, "result_field_solution_sha256": "2" * 64,
        "monitor_owner": "monitor:v54", "result_monitor_owner": "monitor:v54",
        "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
    }
    return {"runs": [{CUTOFF: cutoff, SAR: sar}]}


def test_v54_positive_public_artifacts_are_accepted():
    assert all(validate_public_v54_identity(_payload()).values())


def test_v54_frozen_counterfactuals_are_rejected():
    payload = deepcopy(_payload())
    payload["runs"][0][CUTOFF]["result_mode_impedance_ohm"] = 50.0
    payload["runs"][0][SAR]["result_field_solution_sha256"] = "9" * 64
    assert not all(validate_public_v54_identity(payload).values())


def test_v54_self_consistent_nonphysical_artifacts_are_rejected():
    payload = deepcopy(_payload())
    payload["runs"][0][CUTOFF]["modal_power_w"] = payload["runs"][0][CUTOFF]["result_modal_power_w"] = 2.0
    payload["runs"][0][SAR]["voxel_support_m3"] = payload["runs"][0][SAR]["result_voxel_support_m3"] = [1.0e-6]
    assert not all(validate_public_v54_identity(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload())
    payload["runs"][0][CUTOFF]["sample_frequency_hz"] = [10.0e9]
    payload["runs"][0][SAR]["electric_field_rms_v_m"] = [[12.0]]
    assert not all(validate_public_v54_identity(payload).values())
