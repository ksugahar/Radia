"""Validation-class one-port reflection sweep momentum force.

This example turns a synthetic one-port reflection sweep into a compact force
audit: max/min force frequencies, mean force, reflectance, and passivity flags.

Run:

    python validation_test/rf_waveguide/validation_one_port_reflection_sweep_force.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from result_metadata import add_result_metadata  # noqa: E402

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    C0,
    one_port_reflection_sweep_momentum_force_summary,
)


OUT_JSON = HERE / "validation_one_port_reflection_sweep_force_summary.json"
POWER_W = 1.2
FREQUENCIES_HZ = [1.0e9, 2.0e9, 3.0e9, 4.0e9, 5.0e9]
S11_MAGNITUDES = [0.80, 0.45, 0.05, 0.55, 0.95]
S11_PHASES_DEG = [-150.0, -60.0, 0.0, 50.0, 140.0]


def _polar(magnitude: float, phase_deg: float) -> complex:
    phase = math.radians(phase_deg)
    return magnitude * complex(math.cos(phase), math.sin(phase))


def _assert_close(actual: float, expected: float, atol: float = 1.0e-15) -> float:
    error = abs(actual - expected)
    if error > atol:
        raise AssertionError(f"{actual!r} != {expected!r}; error={error!r}")
    return error


def build_summary() -> dict[str, object]:
    s11 = [_polar(mag, phase) for mag, phase in zip(S11_MAGNITUDES, S11_PHASES_DEG)]
    sweep = one_port_reflection_sweep_momentum_force_summary(
        FREQUENCIES_HZ,
        s11,
        power_incident_W=POWER_W,
        incident_direction=(0.0, 0.0, 1.0),
    )
    passivity_probe = one_port_reflection_sweep_momentum_force_summary(
        [6.0e9],
        [1.002 + 0.0j],
        power_incident_W=POWER_W,
        passivity_tolerance=1.0e-6,
    )

    expected_reflectance = [mag * mag for mag in S11_MAGNITUDES]
    expected_factors = [1.0 + value for value in expected_reflectance]
    expected_mean_reflectance = sum(expected_reflectance) / len(expected_reflectance)
    expected_mean_force = sum(expected_factors) * POWER_W / (len(expected_factors) * C0)
    expected_max_force = max(expected_factors) * POWER_W / C0
    expected_min_force = min(expected_factors) * POWER_W / C0

    errors = {
        "mean_reflectance_error": _assert_close(sweep["mean_reflectance"], expected_mean_reflectance),
        "mean_force_error_N": _assert_close(sweep["mean_force_magnitude_N"], expected_mean_force),
        "max_force_error_N": _assert_close(sweep["max_force_magnitude_N"], expected_max_force),
        "min_force_error_N": _assert_close(sweep["min_force_magnitude_N"], expected_min_force),
    }
    if sweep["max_force_frequency_Hz"] != FREQUENCIES_HZ[-1]:
        raise AssertionError("maximum force should occur at the largest reflection sample")
    if sweep["min_force_frequency_Hz"] != FREQUENCIES_HZ[2]:
        raise AssertionError("minimum force should occur at the reflection notch")
    if not sweep["passivity_ok"]:
        raise AssertionError("synthetic passive sweep should not violate passivity")
    if passivity_probe["passivity_ok"]:
        raise AssertionError("passivity probe should be flagged")

    checks = {
        "n_points": sweep["n_points"],
        "passed": True,
        "max_abs_error": max(errors.values()),
        "mean_reflectance": sweep["mean_reflectance"],
        "expected_mean_reflectance": expected_mean_reflectance,
        "mean_force_magnitude_N": sweep["mean_force_magnitude_N"],
        "expected_mean_force_magnitude_N": expected_mean_force,
        "max_force_frequency_Hz": sweep["max_force_frequency_Hz"],
        "max_force_magnitude_N": sweep["max_force_magnitude_N"],
        "expected_max_force_magnitude_N": expected_max_force,
        "min_force_frequency_Hz": sweep["min_force_frequency_Hz"],
        "min_force_magnitude_N": sweep["min_force_magnitude_N"],
        "expected_min_force_magnitude_N": expected_min_force,
        "passivity_ok": sweep["passivity_ok"],
        "passivity_probe_violation_count": passivity_probe["passivity_violation_count"],
        "passivity_probe_max_excess_magnitude": passivity_probe["max_passivity_excess_magnitude"],
    }

    return {
        "kind": "one_port_reflection_sweep_force_validation",
        "validation_class": True,
        "force_learning": (
            "one-port reflection sweeps map to force extrema through "
            "F=(1+|S11|^2)P_inc/c, while passivity is checked by |S11|<=1"
        ),
        "input": {
            "power_incident_W": POWER_W,
            "frequency_Hz": FREQUENCIES_HZ,
            "s11_magnitude": S11_MAGNITUDES,
            "s11_phase_deg": S11_PHASES_DEG,
        },
        "checks": checks,
        "errors": errors,
        "sweep": sweep,
        "passivity_probe": passivity_probe,
    }


def _json_clean(value):
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = _json_clean(build_summary())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary = add_result_metadata(summary, __file__)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[one-port reflection sweep force]")
    print(f"  points={checks['n_points']}, max_abs_error={checks['max_abs_error']:.3e}")
    print(f"  mean_reflectance={checks['mean_reflectance']:.6g}")
    print(f"  mean_force_N={checks['mean_force_magnitude_N']:.12g}")
    print(f"  max_force: f={checks['max_force_frequency_Hz']:.6g} Hz, F={checks['max_force_magnitude_N']:.12g} N")
    print(f"  min_force: f={checks['min_force_frequency_Hz']:.6g} Hz, F={checks['min_force_magnitude_N']:.12g} N")
    print(f"  passivity_probe_violations={checks['passivity_probe_violation_count']}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
