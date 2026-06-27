"""Validation-class acoustic impedance reflection and absorption gate.

Run:

    python validation_test/acoustic_bem/validation_acoustic_impedance_reflection.py

A locally reacting acoustic boundary can be read as a one-port load.  For a
propagating plane wave, the pressure reflection coefficient is

    Gamma = (Z_load - Z_n) / (Z_load + Z_n),

where Z_n = rho*c/cos(theta) is the normal characteristic impedance.  This
example records matched, mismatched, purely reactive, pressure-release, and
oblique-incidence cases.  The public gate is analytic only: no commercial model
or benchmark value is embedded.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.acoustics import acoustic_impedance_reflection_summary  # noqa: E402


OUT_JSON = HERE / "validation_acoustic_impedance_reflection_summary.json"


def _complex_record(value):
    z = complex(value)
    return {
        "real": 0.0 if z.real == 0.0 else z.real,
        "imag": 0.0 if z.imag == 0.0 else z.imag,
        "abs": abs(z),
    }


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
    p_inc = 2.0
    theta = math.radians(60.0)

    matched = acoustic_impedance_reflection_summary(z0, incident_pressure=p_inc, rho=rho, c=c)
    double_resistance = acoustic_impedance_reflection_summary(2.0 * z0, incident_pressure=p_inc, rho=rho, c=c)
    reactive = acoustic_impedance_reflection_summary(1.0j * z0, incident_pressure=p_inc, rho=rho, c=c)
    pressure_release = acoustic_impedance_reflection_summary(0.0, incident_pressure=p_inc, rho=rho, c=c)
    oblique_matched = acoustic_impedance_reflection_summary(
        z0 / math.cos(theta),
        incidence_angle_rad=theta,
        incident_pressure=p_inc,
        rho=rho,
        c=c,
    )

    _assert_close(abs(matched["pressure_reflection_coefficient"]), 0.0)
    _assert_close(matched["absorption_coefficient"], 1.0)
    _assert_close(double_resistance["pressure_reflection_coefficient"].real, 1.0 / 3.0)
    _assert_close(double_resistance["absorption_coefficient"], 8.0 / 9.0)
    _assert_close(abs(reactive["pressure_reflection_coefficient"]), 1.0)
    _assert_close(reactive["absorbed_intensity"], 0.0)
    _assert_close(pressure_release["pressure_reflection_coefficient"].real, -1.0)
    _assert_close(abs(pressure_release["total_boundary_pressure"]), 0.0)
    _assert_close(abs(oblique_matched["pressure_reflection_coefficient"]), 0.0)
    _assert_close(oblique_matched["absorption_coefficient"], 1.0)

    for name, row in {
        "matched": matched,
        "double_resistance": double_resistance,
        "reactive": reactive,
        "pressure_release": pressure_release,
        "oblique_matched": oblique_matched,
    }.items():
        if abs(row["power_balance_residual"]) > 1.0e-14:
            raise AssertionError(f"{name} power balance residual too large: {row['power_balance_residual']!r}")

    return {
        "kind": "acoustic_impedance_reflection_validation",
        "validation_class": True,
        "rho": rho,
        "c": c,
        "incident_pressure": p_inc,
        "rows": {
            "matched": matched,
            "double_resistance": double_resistance,
            "pure_reactive": reactive,
            "pressure_release": pressure_release,
            "oblique_60deg_matched": oblique_matched,
        },
        "checks": {
            "matched_absorption": matched["absorption_coefficient"],
            "double_resistance_gamma": double_resistance["pressure_reflection_coefficient"],
            "double_resistance_absorption": double_resistance["absorption_coefficient"],
            "reactive_absorbed_intensity": reactive["absorbed_intensity"],
            "reactive_boundary_reactive_intensity": reactive["boundary_reactive_intensity_into_load"],
            "pressure_release_total_pressure_abs": abs(pressure_release["total_boundary_pressure"]),
            "oblique_characteristic_normal_impedance": oblique_matched["characteristic_normal_impedance"],
            "oblique_absorption": oblique_matched["absorption_coefficient"],
        },
    }


def main() -> int:
    summary = build_summary()
    OUT_JSON.write_text(json.dumps(_json_clean(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    gamma = checks["double_resistance_gamma"]
    print("[acoustic impedance reflection]")
    print(f"  matched absorption={checks['matched_absorption']:.12g}")
    print(f"  Z=2Z0 gamma={gamma.real:.12g}{gamma.imag:+.12g}j absorption={checks['double_resistance_absorption']:.12g}")
    print(
        "  reactive absorbed="
        f"{checks['reactive_absorbed_intensity']:.3e}, reactive intensity="
        f"{checks['reactive_boundary_reactive_intensity']:.12e}"
    )
    print(f"  pressure-release |p_total|={checks['pressure_release_total_pressure_abs']:.3e}")
    print(
        "  oblique 60deg Zn="
        f"{checks['oblique_characteristic_normal_impedance']:.12e}, absorption="
        f"{checks['oblique_absorption']:.12g}"
    )
    print(f"  wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
