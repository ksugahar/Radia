"""Validation-class acoustic boundary power from pressure and normal velocity.

Run:

    python examples/acoustic_bem/validation_acoustic_boundary_power.py

Given boundary pressure ``p`` and outward normal velocity ``v_n`` phasors, the
normal complex intensity is

    I_n = alpha p conj(v_n)

with ``alpha=0.5`` for peak phasors and ``alpha=1`` for RMS phasors.  This is
the smallest post-processing gate for acoustic FEM/BEM radiation power: active
power is the real part integrated over the boundary patch, while the imaginary
part records reactive near-field exchange.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.acoustics import acoustic_boundary_power_summary  # noqa: E402


OUT_JSON = HERE / "validation_acoustic_boundary_power_summary.json"


def _complex_record(value):
    if value is None:
        return None
    z = complex(value)
    return {"real": z.real, "imag": z.imag, "abs": abs(z)}


def _json_clean(value):
    if isinstance(value, complex):
        return _complex_record(value)
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, tuple):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-15) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    rho = 1.2041
    c = 343.0
    z0 = rho * c
    area = 0.32
    pressure = 2.0 + 0.0j
    velocity = pressure / z0

    plane_peak = acoustic_boundary_power_summary(pressure, velocity, area=area, amplitude="peak")
    plane_rms = acoustic_boundary_power_summary(pressure, velocity, area=area, amplitude="rms")

    reactive_velocity = 0.03 - 0.01j
    reactive_impedance = 1.0j * 120.0
    reactive = acoustic_boundary_power_summary(
        reactive_impedance * reactive_velocity,
        reactive_velocity,
        area=area,
        amplitude="peak",
    )

    mixed_impedance = z0 * (0.75 + 0.25j)
    mixed_velocity = 0.02 + 0.015j
    mixed = acoustic_boundary_power_summary(
        mixed_impedance * mixed_velocity,
        mixed_velocity,
        area=area,
        amplitude="peak",
    )

    expected_plane_intensity = abs(pressure) ** 2 / (2.0 * z0)
    expected_reactive_intensity = 0.5 * reactive_impedance.imag * abs(reactive_velocity) ** 2
    expected_mixed_active = 0.5 * mixed_impedance.real * abs(mixed_velocity) ** 2
    expected_mixed_reactive = 0.5 * mixed_impedance.imag * abs(mixed_velocity) ** 2

    _assert_close(plane_peak["active_intensity"], expected_plane_intensity)
    _assert_close(plane_peak["active_power"], area * expected_plane_intensity)
    _assert_close(plane_rms["active_power"], 2.0 * plane_peak["active_power"])
    _assert_close(reactive["active_power"], 0.0)
    _assert_close(reactive["reactive_intensity"], expected_reactive_intensity)
    _assert_close(mixed["active_intensity"], expected_mixed_active)
    _assert_close(mixed["reactive_intensity"], expected_mixed_reactive)

    return {
        "kind": "acoustic_boundary_power_validation",
        "validation_class": True,
        "rho": rho,
        "c": c,
        "area": area,
        "rows": {
            "plane_wave_peak": plane_peak,
            "plane_wave_rms": plane_rms,
            "pure_reactive_peak": reactive,
            "mixed_impedance_peak": mixed,
        },
        "checks": {
            "expected_plane_active_intensity": expected_plane_intensity,
            "plane_peak_active_power": plane_peak["active_power"],
            "plane_rms_over_peak_power": plane_rms["active_power"] / plane_peak["active_power"],
            "reactive_active_power_abs": abs(reactive["active_power"]),
            "expected_reactive_intensity": expected_reactive_intensity,
            "mixed_active_intensity": mixed["active_intensity"],
            "expected_mixed_active_intensity": expected_mixed_active,
            "mixed_reactive_intensity": mixed["reactive_intensity"],
            "expected_mixed_reactive_intensity": expected_mixed_reactive,
        },
    }


def main() -> int:
    summary = build_summary()
    OUT_JSON.write_text(json.dumps(_json_clean(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[acoustic boundary power]")
    print(f"  plane peak P={checks['plane_peak_active_power']:.12e} W")
    print(f"  RMS/peak power ratio={checks['plane_rms_over_peak_power']:.6g}")
    print(
        f"  reactive active |P|={checks['reactive_active_power_abs']:.3e}, "
        f"reactive intensity={checks['expected_reactive_intensity']:.12e} W/m^2"
    )
    print(
        f"  mixed active/reactive={checks['mixed_active_intensity']:.12e}/"
        f"{checks['mixed_reactive_intensity']:.12e} W/m^2"
    )
    print(f"  wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
