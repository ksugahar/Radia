from __future__ import annotations

import math

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v41 import _identity_v41


_FORCE = (
    "weightedstress_airgapcontour_virtualwork_force_direction_mesh_"
    "fieldowner_result_generation_identity"
)
_CONDUCTOR = (
    "harmonicconductor_skin_depth_complexcurrent_jouleloss_impedance_"
    "power_mesh_result_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v42_public_weightedstress_force_airgap_contour_independence_virtualwork_mesh_mismatch",
    "v42_public_harmonicconductor_skin_depth_complexcurrent_jouleloss_impedance_power_mismatch",
)


def _identity_v42():
    identity = _identity_v41()
    generation = "force-closure-842"
    force = [12.0, -3.0]
    contours = [[12.02, -2.99], [11.98, -3.01], [12.01, -3.0]]
    force_norm = math.hypot(*force)
    spread = max(
        math.hypot(sample[0] - force[0], sample[1] - force[1])
        for sample in contours
    ) / force_norm
    displacement = 1.0e-4
    energy_minus = 1.0012
    energy_plus = 0.9988
    virtual_force = -(energy_plus - energy_minus) / (2.0 * displacement)
    mesh_forces = [[11.9, -2.95], force]
    identity[_FORCE] = {
        "force_generation": generation,
        **{
            key: generation
            for key in (
                "weighted_stress_generation",
                "airgap_contour_generation",
                "virtual_work_generation",
                "direction_generation",
                "mesh_generation",
                "field_generation",
                "result_generation",
            )
        },
        "weighted_stress_force_n": force,
        "result_weighted_stress_force_n": force,
        "airgap_contour_forces_n": contours,
        "result_airgap_contour_forces_n": contours,
        "contour_independence_relative_spread": spread,
        "result_contour_independence_relative_spread": spread,
        "virtual_work_displacement_m": displacement,
        "result_virtual_work_displacement_m": displacement,
        "energy_minus_j": energy_minus,
        "result_energy_minus_j": energy_minus,
        "energy_plus_j": energy_plus,
        "result_energy_plus_j": energy_plus,
        "virtual_work_direction": [1.0, 0.0],
        "result_virtual_work_direction": [1.0, 0.0],
        "virtual_work_force_n": virtual_force,
        "result_virtual_work_force_n": virtual_force,
        "mesh_refinement_force_samples_n": mesh_forces,
        "result_mesh_refinement_force_samples_n": mesh_forces,
        "field_owner": "field:force-closure-842",
        "accepted_field_owner": "field:force-closure-842",
        "mesh_owner": "mesh:force-closure-842",
        "accepted_mesh_owner": "mesh:force-closure-842",
        "force_result_sha256": "1" * 64,
        "accepted_force_result_sha256": "1" * 64,
    }

    generation = "harmonic-conductor-842"
    frequency = 1000.0
    conductivity = 5.8e7
    permeability = 4.0e-7 * math.pi
    skin_depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * permeability * conductivity))
    current_density = [[1.0e6, -2.0e5], [0.8e6, -1.5e5]]
    integrated_j_squared = 2.0 * conductivity * 20.8
    identity[_CONDUCTOR] = {
        "conductor_generation": generation,
        **{
            key: generation
            for key in (
                "skin_depth_generation",
                "current_density_generation",
                "joule_generation",
                "impedance_generation",
                "power_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "conductivity_s_per_m": conductivity,
        "result_conductivity_s_per_m": conductivity,
        "absolute_permeability_h_per_m": permeability,
        "result_absolute_permeability_h_per_m": permeability,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "complex_current_density_a_per_m2": current_density,
        "result_complex_current_density_a_per_m2": current_density,
        "integrated_abs_current_density_sq_a2_per_m": integrated_j_squared,
        "result_integrated_abs_current_density_sq_a2_per_m": integrated_j_squared,
        "joule_loss_w": 20.8,
        "result_joule_loss_w": 20.8,
        "terminal_current_a": [10.0, -2.0],
        "result_terminal_current_a": [10.0, -2.0],
        "terminal_voltage_v": [4.4, 1.2],
        "result_terminal_voltage_v": [4.4, 1.2],
        "terminal_impedance_ohm": [0.4, 0.2],
        "result_terminal_impedance_ohm": [0.4, 0.2],
        "complex_power_va": [20.8, 10.4],
        "result_complex_power_va": [20.8, 10.4],
        "mesh_owner": "mesh:harmonic-conductor-842",
        "accepted_mesh_owner": "mesh:harmonic-conductor-842",
        "conductor_result_sha256": "2" * 64,
        "accepted_conductor_result_sha256": "2" * 64,
    }
    return identity


def test_v42_public_positive_force_and_harmonic_conductor_closure():
    assert _gate(_identity_v42())["status"] == "ok"


def test_v42_public_weighted_stress_force_mismatch():
    identity = _identity_v42()
    identity[_FORCE].update(
        {
            "airgap_contour_generation": "force-closure-841",
            "mesh_generation": "force-closure-840",
            "result_generation": "force-closure-839",
            "result_weighted_stress_force_n": [-12.0, 3.0],
            "result_airgap_contour_forces_n": [[30.0, 0.0]],
            "result_contour_independence_relative_spread": 2.0,
            "result_energy_plus_j": 1.0012,
            "result_virtual_work_direction": [-1.0, 0.0],
            "result_virtual_work_force_n": -12.0,
            "result_mesh_refinement_force_samples_n": [[1.0, 1.0], [-12.0, 3.0]],
            "accepted_field_owner": "stale:field",
            "accepted_mesh_owner": "stale:mesh",
            "accepted_force_result_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "weighted_stress_forces_close_airgap_contours_virtual_work_direction_mesh_field_owner_and_result"
    ]


def test_v42_public_harmonic_conductor_mismatch():
    identity = _identity_v42()
    identity[_CONDUCTOR].update(
        {
            "skin_depth_generation": "harmonic-conductor-841",
            "power_generation": "harmonic-conductor-840",
            "result_generation": "harmonic-conductor-839",
            "result_skin_depth_m": 0.1,
            "result_complex_current_density_a_per_m2": [[-1.0e6, 2.0e5]],
            "result_integrated_abs_current_density_sq_a2_per_m": -1.0,
            "result_joule_loss_w": -20.8,
            "result_terminal_impedance_ohm": [-0.4, -0.2],
            "result_complex_power_va": [-20.8, 10.4],
            "accepted_mesh_owner": "stale:mesh",
            "accepted_conductor_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_conductors_close_skin_depth_complex_current_joule_impedance_power_mesh_and_result"
    ]


def test_v42_public_rejects_self_consistent_wrong_contour_spread():
    identity = _identity_v42()
    record = identity[_FORCE]
    record["airgap_contour_forces_n"] = [[30.0, 0.0]] * 3
    record["result_airgap_contour_forces_n"] = [[30.0, 0.0]] * 3
    spread = math.hypot(18.0, 3.0) / math.hypot(12.0, -3.0)
    record["contour_independence_relative_spread"] = spread
    record["result_contour_independence_relative_spread"] = spread
    assert _gate(identity)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_complex_power():
    identity = _identity_v42()
    record = identity[_CONDUCTOR]
    record["complex_power_va"] = [20.8, -10.4]
    record["result_complex_power_va"] = [20.8, -10.4]
    assert _gate(identity)["status"] == "needs_attention"
