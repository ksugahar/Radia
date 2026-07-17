from __future__ import annotations

from test_femm_generalization_v19 import _gate
from test_femm_generalization_v27 import _identity_v27
from test_force_coenergy_gate import _quadratic_case


def _identity_v28(sample_count):
    identity = _identity_v27(sample_count)
    identity["axisymmetric_aphi_radial_weight_region_energy_force_mesh_solution_generation_identity"] = {
        "axisymmetric_generation": "axisym-force-151", "aphi_axisymmetric_generation": "axisym-force-151",
        "weight_axisymmetric_generation": "axisym-force-151", "region_axisymmetric_generation": "axisym-force-151",
        "energy_axisymmetric_generation": "axisym-force-151", "force_axisymmetric_generation": "axisym-force-151",
        "mesh_axisymmetric_generation": "axisym-force-151", "solution_axisymmetric_generation": "axisym-force-151",
        "result_axisymmetric_generation": "axisym-force-151", "formulation": "Aphi", "result_formulation": "Aphi",
        "radial_weighting": "2*pi*r", "result_radial_weighting": "2*pi*r",
        "region_ids": [4, 5], "result_region_ids": [4, 5],
        "displacement_m": [-0.0001, 0.0, 0.0001],
        "magnetic_coenergy_j": [0.4982, 0.5, 0.5018],
        "result_magnetic_coenergy_j": [0.4982, 0.5, 0.5018],
        "force_axis": "z", "result_force_axis": "z",
        "force_from_energy_n": 18.0, "result_force_from_energy_n": 18.0,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "solution_sha256": "2" * 64, "accepted_solution_sha256": "2" * 64,
    }
    identity["permanent_magnet_recoil_temperature_operating_point_demag_force_generation_identity"] = {
        "magnet_generation": "pm-operating-point-151", "recoil_magnet_generation": "pm-operating-point-151",
        "temperature_magnet_generation": "pm-operating-point-151", "operating_point_magnet_generation": "pm-operating-point-151",
        "frame_magnet_generation": "pm-operating-point-151", "demag_magnet_generation": "pm-operating-point-151",
        "force_magnet_generation": "pm-operating-point-151", "mesh_magnet_generation": "pm-operating-point-151",
        "result_magnet_generation": "pm-operating-point-151",
        "recoil_relative_permeability": 1.05, "result_recoil_relative_permeability": 1.05,
        "remanence_t": 1.2, "result_remanence_t": 1.2,
        "magnet_temperature_c": 80.0, "result_magnet_temperature_c": 80.0,
        "operating_point_bh": [0.82, -302000.0], "result_operating_point_bh": [0.82, -302000.0],
        "demag_margin_a_m": 95000.0, "result_demag_margin_a_m": 95000.0,
        "magnetization_frame_sha256": "3" * 64, "result_magnetization_frame_sha256": "3" * 64,
        "force_n": [12.0, 0.2], "result_force_n": [12.0, 0.2],
        "mesh_sha256": "4" * 64, "result_mesh_sha256": "4" * 64,
        "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
    }
    return identity


def test_v28_public_positive_axisymmetric_and_pm_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v28(len(positions)))["status"] == "ok"


def test_v28_public_axisymmetric_aphi_radial_weighting_region_force_energy_mesh_solution_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v28(len(positions))
    identity["axisymmetric_aphi_radial_weight_region_energy_force_mesh_solution_generation_identity"].update(
        {"aphi_axisymmetric_generation": "axisym-force-150", "result_formulation": "Az",
         "result_radial_weighting": "1", "result_region_ids": [5, 6],
         "result_magnetic_coenergy_j": [0.49, 0.5, 0.52], "result_force_axis": "r",
         "result_force_from_energy_n": 55.0, "accepted_solution_sha256": "a" * 64}
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["axisymmetric_force_uses_current_aphi_weight_regions_energy_axis_mesh_and_solution"]


def test_v28_public_permanent_magnet_recoil_line_temperature_operating_point_demag_margin_force_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v28(len(positions))
    identity["permanent_magnet_recoil_temperature_operating_point_demag_force_generation_identity"].update(
        {"recoil_magnet_generation": "pm-operating-point-150", "result_recoil_relative_permeability": 1.2,
         "result_remanence_t": 1.0, "result_magnet_temperature_c": 20.0,
         "result_operating_point_bh": [0.45, -500000.0], "result_demag_margin_a_m": -10000.0,
         "result_force_n": [8.0, -1.0], "accepted_result_sha256": "d" * 64}
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["permanent_magnet_force_uses_current_recoil_temperature_operating_point_frame_demag_and_mesh"]
