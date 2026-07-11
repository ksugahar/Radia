from __future__ import annotations

import copy
import math

from radia_mcp.radia_ngsolve.harmonic_zero_net_circuit_gate import (
    harmonic_zero_net_circuit_gate,
)


def _complex(value: complex) -> dict[str, float]:
    return {"real": value.real, "imag": value.imag}


def _summary() -> dict:
    omega = 2.0 * math.pi * 5.0
    rows = []
    for amplitude, flux in (
        (0.5, 0.0009 - 0.0001j),
        (1.0, 0.0017 - 0.0004j),
        (1.5, 0.00235 - 0.00078j),
    ):
        rows.append(
            {
                "source_amplitude": amplitude,
                "source_current_a": 1290.32 * amplitude,
                "constrained_circuit_current_a": _complex(0.0j),
                "circuit_voltage_v": _complex(1j * omega * flux),
                "circuit_flux_linkage_wb_turn": _complex(flux),
                "conductive_losses_w": [
                    _complex((4.0 + amplitude) * amplitude**2 + 1.0e-18j),
                    _complex((8.0 + amplitude) * amplitude**2 - 1.0e-18j),
                ],
                "dc_force_components_n": [-40.0 * amplitude**1.8, 5.0 * amplitude**1.7],
                "two_x_force_phasors_n": [
                    _complex((-12.0 - 25.0j) * amplitude**1.8),
                    _complex((0.2 + 3.0j) * amplitude**1.7),
                ],
                "node_count": 1000,
                "element_count": 1900,
            }
        )
    return {
        "contract": {
            "frequency_hz": 5.0,
            "phasor_convention": "exp(+j*omega*t)",
            "faraday_identity": "V=+j*omega*flux_linkage",
            "circuit_constraint": "zero_net_current",
            "material_response": "nonlinear",
            "scaling_interpretation": "diagnostic_only",
            "force_components": ["dc_time_average", "two_x_phasor"],
            "force_method": "lorentz_volume_current_density",
            "force_scope": "current_density_contribution_not_total_ferromagnetic_force",
        },
        "rows": rows,
    }


def test_harmonic_zero_net_circuit_gate_accepts_consistent_phasors() -> None:
    result = harmonic_zero_net_circuit_gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["faraday_voltage_flux_identity"] is True
    assert result["checks"]["nonlinear_scaling_not_overclaimed"] is True
    assert result["metrics"]["max_loss_per_amplitude_squared_relative_span"] > 0.0


def test_harmonic_zero_net_circuit_gate_rejects_wrong_sign_and_current_leak() -> None:
    payload = copy.deepcopy(_summary())
    payload["rows"][1]["circuit_voltage_v"]["imag"] *= -1.0
    payload["rows"][1]["constrained_circuit_current_a"] = _complex(0.2 + 0.1j)
    result = harmonic_zero_net_circuit_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["faraday_voltage_flux_identity"] is False
    assert result["checks"]["zero_net_current_constraint"] is False


def test_harmonic_zero_net_circuit_gate_rejects_square_law_overclaim() -> None:
    payload = _summary()
    payload["contract"]["scaling_interpretation"] = "must_follow_square_law"
    result = harmonic_zero_net_circuit_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["nonlinear_scaling_not_overclaimed"] is False


def test_harmonic_zero_net_circuit_gate_rejects_total_force_overclaim() -> None:
    payload = _summary()
    payload["contract"]["force_scope"] = "total_ferromagnetic_force"
    result = harmonic_zero_net_circuit_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["force_scope_not_overclaimed"] is False
