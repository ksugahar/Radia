"""Validation-class acoustic radiation example: uniformly pulsating sphere.

This is intentionally an example/validation problem, not just a unit test.  It
sweeps from low ``ka`` to high ``ka`` and records the anchors a scalar acoustic
FEM/BEM implementation should reproduce before moving to general geometry:

* radiation resistance scales like ``(ka)^2`` at low frequency;
* reactive specific impedance scales like ``ka`` at low frequency;
* active acoustic power is conserved on every spherical observation surface;
* outgoing pressure amplitude decays exactly as ``1/r``;
* high ``ka`` approaches the plane-wave specific resistance ``rho*c``.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.acoustics import pulsating_sphere_radiation


def _complex_record(value):
    z = complex(value)
    return {
        "real": z.real,
        "imag": z.imag,
        "abs": abs(z),
        "phase_rad": cmath.phase(z),
    }


def _case(name, ka, radius, surface_velocity, rho, c):
    frequency = ka * c / (2.0 * math.pi * radius)
    surface = pulsating_sphere_radiation(
        radius,
        frequency,
        surface_velocity,
        rho=rho,
        c=c,
        sample_radius=radius,
    )
    samples = []
    max_power_rel_error = 0.0
    for factor in (2.0, 5.0, 10.0, 20.0):
        out = pulsating_sphere_radiation(
            radius,
            frequency,
            surface_velocity,
            rho=rho,
            c=c,
            sample_radius=factor * radius,
        )
        rel = abs(out["sample_power"] - out["radiated_power"]) / max(out["radiated_power"], 1e-300)
        max_power_rel_error = max(max_power_rel_error, rel)
        samples.append(
            {
                "radius_over_a": factor,
                "pressure": _complex_record(out["sample_pressure"]),
                "radial_velocity": _complex_record(out["sample_radial_velocity"]),
                "specific_impedance": _complex_record(out["sample_specific_impedance"]),
                "intensity": out["sample_intensity"],
                "sample_power": out["sample_power"],
                "power_rel_error": rel,
            }
        )

    p10 = samples[2]["pressure"]["abs"]
    p20 = samples[3]["pressure"]["abs"]
    return {
        "name": name,
        "ka": surface["ka"],
        "frequency_hz": frequency,
        "specific_impedance": _complex_record(surface["specific_impedance"]),
        "radiation_efficiency": surface["radiation_efficiency"],
        "reactance_ratio": surface["reactance_ratio"],
        "volume_velocity_impedance": _complex_record(surface["volume_velocity_impedance"]),
        "radiated_power_w": surface["radiated_power"],
        "max_power_rel_error": max_power_rel_error,
        "pressure_20a_over_10a": p20 / p10,
        "samples": samples,
    }


def build_summary(quick=False):
    radius = 0.1
    surface_velocity = 0.01
    rho = 1.2041
    c = 343.0
    cases = [
        _case("deep_subwavelength", 0.025, radius, surface_velocity, rho, c),
        _case("low_frequency", 0.050, radius, surface_velocity, rho, c),
        _case("transition", 1.00, radius, surface_velocity, rho, c),
    ]
    if not quick:
        cases.extend(
            [
                _case("moderately_large", 5.00, radius, surface_velocity, rho, c),
                _case("plane_wave_limit", 10.00, radius, surface_velocity, rho, c),
            ]
        )

    low_a, low_b = cases[0], cases[1]
    resistance_scaling = (
        low_b["specific_impedance"]["real"] / low_a["specific_impedance"]["real"]
    )
    reactance_scaling = (
        low_b["specific_impedance"]["imag"] / low_a["specific_impedance"]["imag"]
    )
    max_power_rel_error = max(row["max_power_rel_error"] for row in cases)
    max_pressure_decay_error = max(abs(row["pressure_20a_over_10a"] - 0.5) for row in cases)

    checks = {
        "low_frequency_resistance_ratio_expected": 4.0,
        "low_frequency_resistance_ratio": resistance_scaling,
        "low_frequency_reactance_ratio_expected": 2.0,
        "low_frequency_reactance_ratio": reactance_scaling,
        "max_power_rel_error": max_power_rel_error,
        "max_pressure_decay_abs_error": max_pressure_decay_error,
        "largest_ka_radiation_efficiency": cases[-1]["radiation_efficiency"],
    }

    if abs(resistance_scaling - 4.0) > 2.0e-2:
        raise AssertionError(f"low-frequency resistance scaling drifted: {resistance_scaling}")
    if abs(reactance_scaling - 2.0) > 2.0e-2:
        raise AssertionError(f"low-frequency reactance scaling drifted: {reactance_scaling}")
    if max_power_rel_error > 1.0e-12:
        raise AssertionError(f"spherical power conservation drifted: {max_power_rel_error}")
    if max_pressure_decay_error > 1.0e-14:
        raise AssertionError(f"pressure 1/r decay drifted: {max_pressure_decay_error}")
    if not quick and cases[-1]["radiation_efficiency"] < 0.99:
        raise AssertionError("high-ka radiation efficiency did not approach rho*c")

    return {
        "problem": "uniformly pulsating sphere acoustic radiation",
        "time_convention": "peak phasors with exp(+i omega t)",
        "radius_m": radius,
        "surface_velocity_peak_m_per_s": surface_velocity,
        "rho_kg_per_m3": rho,
        "c_m_per_s": c,
        "checks": checks,
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="Run the three-case smoke sweep")
    parser.add_argument("--output", default=str(Path(__file__).with_name("validation_pulsating_sphere_radiation_summary.json")))
    args = parser.parse_args()

    summary = build_summary(quick=args.quick)
    out = Path(args.output)
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary["checks"], indent=2, sort_keys=True))
    print(f"wrote {os.fspath(out)}")


if __name__ == "__main__":
    main()
