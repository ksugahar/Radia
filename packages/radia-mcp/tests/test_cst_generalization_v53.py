from copy import deepcopy

from radia_mcp.radia_ngsolve.wave_port_identity_v53 import FARFIELD, WAVEGUIDE, validate_public_v53_identity


PROMOTED_CASE_IDS = {"v53_public_waveguide_mode_cutoff_normalization_port_referenceplane_owner_mismatch", "v53_public_farfield_realizedgain_polarization_basis_angulargrid_owner_mismatch"}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]: return {"generation": generation, **{field: generation for field in fields}}


def _payload():
    width = 0.02286; cutoff = 299792458.0 / (2.0 * width)
    waveguide = {**_generations("wave-v53", ("mode_generation", "cutoff_generation", "normalization_generation", "plane_generation", "owner_generation", "result_generation")), "mode_id": "TE10", "result_mode_id": "TE10", "waveguide_width_m": width, "result_waveguide_width_m": width, "cutoff_frequency_hz": cutoff, "result_cutoff_frequency_hz": cutoff, "sample_frequency_hz": 10.0e9, "result_sample_frequency_hz": 10.0e9, "normalization": "unit_incident_power_w", "result_normalization": "unit_incident_power_w", "reference_plane_offset_m": -0.01, "result_reference_plane_offset_m": -0.01, "port_owner": "port:v53", "result_port_owner": "port:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    theta = [0.0, 90.0, 180.0]; phi = [0.0, 90.0, 180.0, 270.0]; gain = [[1.0] * 4, [8.0] * 4, [1.0] * 4]
    farfield = {**_generations("far-v53", ("gain_generation", "polarization_generation", "grid_generation", "owner_generation", "result_generation")), "theta_deg": theta, "result_theta_deg": theta, "phi_deg": phi, "result_phi_deg": phi, "realized_gain_dbi": gain, "result_realized_gain_dbi": gain, "polarization_basis": "ludwig3", "result_polarization_basis": "ludwig3", "monitor_owner": "monitor:v53", "result_monitor_owner": "monitor:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {"runs": [{WAVEGUIDE: waveguide, FARFIELD: farfield}]}


def test_v53_positive_public_artifacts_are_accepted(): assert all(validate_public_v53_identity(_payload()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    payload = deepcopy(_payload()); payload["runs"][0][WAVEGUIDE]["result_mode_id"] = "TE20"; payload["runs"][0][FARFIELD]["result_polarization_basis"] = "spherical"; assert not all(validate_public_v53_identity(payload).values())


def test_v53_self_consistent_invalid_semantics_are_rejected():
    payload = deepcopy(_payload()); payload["runs"][0][WAVEGUIDE]["cutoff_frequency_hz"] = payload["runs"][0][WAVEGUIDE]["result_cutoff_frequency_hz"] = 1.0; payload["runs"][0][FARFIELD]["realized_gain_dbi"] = payload["runs"][0][FARFIELD]["result_realized_gain_dbi"] = [[1.0]]; assert not all(validate_public_v53_identity(payload).values())
