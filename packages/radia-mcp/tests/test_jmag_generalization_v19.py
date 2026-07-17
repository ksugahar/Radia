from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_pwm_controlled_motor_loss_gate import _payload, _with_artifact_identity


def _payload_v19():
    payload = _with_artifact_identity(_payload())
    identity = payload["artifact_identity"]
    identity["iron_loss_harmonic_frequency_coefficient_unit_basis_identity"] = {
        "loss_generation": "iron-loss-21",
        "harmonic_loss_generation": "iron-loss-21",
        "coefficient_loss_generation": "iron-loss-21",
        "frequency_unit": "Hz",
        "coefficient_frequency_unit": "Hz",
        "flux_density_unit": "T",
        "coefficient_flux_density_unit": "T",
        "harmonic_frequencies_hz": [50.0, 150.0, 250.0],
        "evaluated_harmonic_frequencies_hz": [50.0, 150.0, 250.0],
        "loss_coefficients": [1.0, 0.02, 0.001],
        "evaluated_loss_coefficients": [1.0, 0.02, 0.001],
        "loss_basis_sha256": "a" * 64,
        "evaluated_loss_basis_sha256": "a" * 64,
    }
    identity["demagnetization_temperature_current_phase_operating_point_identity"] = {
        "operating_point_generation": "operating-point-21",
        "temperature_operating_point_generation": "operating-point-21",
        "current_phase_operating_point_generation": "operating-point-21",
        "demag_margin_operating_point_generation": "operating-point-21",
        "magnet_temperature_c": 120.0,
        "demag_margin_temperature_c": 120.0,
        "current_phase_deg": 90.0,
        "demag_margin_current_phase_deg": 90.0,
        "operating_point_sha256": "b" * 64,
        "demag_margin_operating_point_sha256": "b" * 64,
    }
    return payload


def test_v19_public_positive_iron_loss_and_demag_operating_point_identity():
    result = pwm_controlled_motor_loss_gate(_payload_v19())
    assert result["status"] == "ok"
    assert result["checks"][
        "iron_loss_harmonics_and_coefficients_share_frequency_flux_units"
    ]
    assert result["checks"][
        "demag_margin_uses_current_temperature_and_current_phase_state"
    ]


def test_v19_public_iron_loss_harmonic_frequency_basis_coefficient_unit_mismatch():
    payload = _payload_v19()
    payload["artifact_identity"][
        "iron_loss_harmonic_frequency_coefficient_unit_basis_identity"
    ].update(
        {
            "coefficient_loss_generation": "iron-loss-20",
            "coefficient_frequency_unit": "kHz",
            "coefficient_flux_density_unit": "mT",
            "evaluated_harmonic_frequencies_hz": [0.05, 0.15, 0.25],
            "evaluated_loss_coefficients": [1000.0, 20.0, 1.0],
            "evaluated_loss_basis_sha256": "e" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "iron_loss_harmonics_and_coefficients_share_frequency_flux_units"
    ] is False


def test_v19_public_demagnetization_operating_point_temperature_phase_generation_mismatch():
    payload = _payload_v19()
    payload["artifact_identity"][
        "demagnetization_temperature_current_phase_operating_point_identity"
    ].update(
        {
            "temperature_operating_point_generation": "operating-point-20",
            "demag_margin_operating_point_generation": "operating-point-20",
            "demag_margin_temperature_c": 80.0,
            "demag_margin_current_phase_deg": 60.0,
            "demag_margin_operating_point_sha256": "e" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "demag_margin_uses_current_temperature_and_current_phase_state"
    ] is False
