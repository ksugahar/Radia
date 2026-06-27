"""Validation-class baffled circular-piston acoustic radiation example.

This example records the closed-form radiation impedance of a uniformly
vibrating circular piston in an infinite baffle.  It is a practical
FEM/BEM-style gate for flat acoustic boundaries such as speakers,
transducers, and duct openings:

* radiation resistance scales as ``(ka)^2 / 2`` at low frequency;
* radiation reactance scales as ``8 ka / (3 pi)`` at low frequency;
* high ``ka`` approaches the plane-wave specific resistance ``rho c``;
* active radiated power follows ``0.5 * S * Re(z) * |v0|^2``.

Run:

    python examples/acoustic_bem/validation_baffled_piston_radiation.py
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.acoustics import baffled_circular_piston_radiation  # noqa: E402


OUT_JSON = Path(__file__).with_name("validation_baffled_piston_radiation_summary.json")


def _complex_record(value):
    z = complex(value)
    return {
        "real": z.real,
        "imag": z.imag,
        "abs": abs(z),
        "phase_rad": cmath.phase(z),
    }


def _case(ka, radius, velocity, rho, c):
    frequency = ka * c / (2.0 * math.pi * radius)
    out = baffled_circular_piston_radiation(radius, frequency, velocity, rho=rho, c=c)
    return {
        "ka": out["ka"],
        "frequency_hz": out["frequency"],
        "specific_impedance": _complex_record(out["specific_impedance"]),
        "radiation_efficiency": out["radiation_efficiency"],
        "reactance_ratio": out["reactance_ratio"],
        "low_ka_resistance_asymptote": out["low_ka_resistance_asymptote"],
        "low_ka_reactance_asymptote": out["low_ka_reactance_asymptote"],
        "volume_velocity_impedance": _complex_record(out["volume_velocity_impedance"]),
        "radiated_power_w": out["radiated_power"],
        "power_formula_w": 0.5 * out["surface_area"] * out["specific_resistance"] * abs(velocity) ** 2,
    }


def build_summary():
    radius = 0.08
    velocity = 0.02
    rho = 1.2041
    c = 343.0
    kas = (0.025, 0.050, 0.10, 0.50, 1.0, 5.0, 10.0, 20.0)
    rows = [_case(ka, radius, velocity, rho, c) for ka in kas]
    low1, low2 = rows[0], rows[1]
    resistance_ratio = low2["radiation_efficiency"] / low1["radiation_efficiency"]
    reactance_ratio = low2["reactance_ratio"] / low1["reactance_ratio"]
    power_errors = [
        abs(row["radiated_power_w"] - row["power_formula_w"])
        / max(row["radiated_power_w"], 1.0e-300)
        for row in rows
    ]
    low_asymptote_errors = {
        "resistance_rel_error": abs(
            low1["radiation_efficiency"] - low1["low_ka_resistance_asymptote"]
        ) / low1["low_ka_resistance_asymptote"],
        "reactance_rel_error": abs(
            low1["reactance_ratio"] - low1["low_ka_reactance_asymptote"]
        ) / low1["low_ka_reactance_asymptote"],
    }
    checks = {
        "low_resistance_ratio_expected": 4.0,
        "low_resistance_ratio": resistance_ratio,
        "low_reactance_ratio_expected": 2.0,
        "low_reactance_ratio": reactance_ratio,
        "low_asymptote_errors": low_asymptote_errors,
        "max_power_formula_rel_error": max(power_errors),
        "ka20_radiation_efficiency": rows[-1]["radiation_efficiency"],
        "ka20_reactance_ratio": rows[-1]["reactance_ratio"],
    }
    assert abs(resistance_ratio - 4.0) < 5.0e-3
    assert abs(reactance_ratio - 2.0) < 5.0e-3
    assert low_asymptote_errors["resistance_rel_error"] < 1.0e-3
    assert low_asymptote_errors["reactance_rel_error"] < 1.0e-3
    assert checks["max_power_formula_rel_error"] < 1.0e-14
    assert checks["ka20_radiation_efficiency"] > 0.95

    return {
        "problem": "baffled circular piston acoustic radiation impedance",
        "time_convention": "peak phasors with exp(+i omega t)",
        "radius_m": radius,
        "surface_velocity_peak_m_per_s": velocity,
        "rho_kg_per_m3": rho,
        "c_m_per_s": c,
        "rows": rows,
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT_JSON)
    args = parser.parse_args()
    summary = build_summary()
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary["checks"], indent=2, sort_keys=True))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
