from __future__ import annotations

from copy import deepcopy
from math import log10, sqrt

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v41 import _summary_v41


_PROMOTED_CASE_IDS = (
    "v42_public_lowfrequency_bem_stabilized_kernel_condition_staticlimit_charge_energy_mismatch",
    "v42_public_duct_scattering_transmission_reflection_loss_power_modal_balance_mismatch",
)
_LOW_KEY = (
    "lowfrequency_bem_stabilized_kernel_staticlimit_condition_charge_energy_"
    "boundarymesh_result_identity"
)
_DUCT_KEY = (
    "duct_scattering_mode_pressure_transmission_reflection_loss_power_mesh_"
    "result_identity"
)


def _summary_v42() -> dict:
    payload = deepcopy(_summary_v41())
    generation = "low-frequency-bem-842"
    values = {
        "frequency_hz": [1.0e-3, 1.0e-2, 1.0e-1],
        "stabilized_kernel": "static_dynamic_split_p1",
        "dynamic_kernel_correction": [1.0e-6, 1.0e-5, 1.0e-4],
        "static_limit_residual": [1.0e-8, 1.0e-7, 1.0e-6],
        "condition_estimate": [100.0, 80.0, 60.0],
        "boundary_charge_c": [1.0e-9, -1.0e-9],
        "boundary_potential_v": [2.0, -2.0],
        "potential_energy_j": 2.0e-9,
    }
    payload[_LOW_KEY] = {
        "lowfrequency_bem_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "kernel_generation", "staticlimit_generation",
            "condition_generation", "charge_generation", "energy_generation",
            "mesh_generation", "result_generation",
        )},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "boundary_mesh_owner": "boundary-mesh:low-frequency-bem-842",
        "accepted_boundary_mesh_owner": "boundary-mesh:low-frequency-bem-842",
        "bem_result_sha256": "1" * 64,
        "accepted_bem_result_sha256": "1" * 64,
    }

    generation = "duct-scattering-842"
    values = {
        "frequency_hz": 1000.0,
        "incident_mode_pressure": [1.0, 0.0],
        "reflected_mode_pressure": [0.2, 0.0],
        "transmitted_mode_pressure": [sqrt(0.96), 0.0],
        "incident_modal_power_w": 1.0,
        "reflected_modal_power_w": 0.04,
        "transmitted_modal_power_w": 0.96,
        "dissipated_power_w": 0.0,
        "transmission_loss_db": -10.0 * log10(0.96),
        "modal_power_balance_residual_w": 0.0,
    }
    payload[_DUCT_KEY] = {
        "duct_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "mode_generation", "pressure_generation",
            "reflection_generation", "transmission_generation", "loss_generation",
            "power_generation", "mesh_generation", "result_generation",
        )},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:duct-scattering-842",
        "accepted_mesh_owner": "mesh:duct-scattering-842",
        "duct_result_sha256": "2" * 64,
        "accepted_duct_result_sha256": "2" * 64,
    }
    return payload


def test_v42_public_positive_lowfrequency_bem_and_duct_scattering_closure() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v42())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v42_public_lowfrequency_bem_stabilized_kernel_condition_staticlimit_charge_energy_mismatch() -> None:
    payload = _summary_v42()
    payload[_LOW_KEY].update({
        "kernel_generation": "low-frequency-bem-841",
        "energy_generation": "low-frequency-bem-840",
        "result_generation": "low-frequency-bem-839",
        "result_stabilized_kernel": "unstabilized",
        "result_static_limit_residual": [1.0],
        "result_condition_estimate": [-1.0],
        "result_boundary_charge_c": [1.0],
        "result_potential_energy_j": -1.0,
        "accepted_boundary_mesh_owner": "stale:mesh",
        "accepted_bem_result_sha256": "9" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "lowfrequency_bem_uses_current_stabilized_kernel_static_limit_condition_charge_energy_mesh_owner_and_result"
    ]


def test_v42_public_duct_scattering_transmission_reflection_loss_power_modal_balance_mismatch() -> None:
    payload = _summary_v42()
    payload[_DUCT_KEY].update({
        "mode_generation": "duct-scattering-841",
        "power_generation": "duct-scattering-840",
        "result_generation": "duct-scattering-839",
        "result_incident_mode_pressure": [0.0, 0.0],
        "result_reflected_mode_pressure": [2.0, 0.0],
        "result_transmitted_mode_pressure": [-1.0, 0.0],
        "result_transmission_loss_db": -10.0,
        "result_modal_power_balance_residual_w": 1.0,
        "accepted_mesh_owner": "stale:mesh",
        "accepted_duct_result_sha256": "a" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "duct_scattering_uses_current_modes_pressures_transmission_reflection_loss_power_mesh_owner_and_result"
    ]


def test_v42_public_rejects_self_consistent_non_neutral_lowfrequency_charge() -> None:
    payload = _summary_v42()
    payload[_LOW_KEY]["boundary_charge_c"] = [1.0e-9, 1.0e-9]
    payload[_LOW_KEY]["result_boundary_charge_c"] = [1.0e-9, 1.0e-9]
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_duct_power_imbalance() -> None:
    payload = _summary_v42()
    payload[_DUCT_KEY]["reflected_modal_power_w"] = 0.4
    payload[_DUCT_KEY]["result_reflected_modal_power_w"] = 0.4
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
