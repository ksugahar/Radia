"""Validation-class torque waveform harmonic health example.

Run:

    python validation_test/force_validation/validation_torque_waveform_health.py

The example uses an analytic torque-angle waveform with a mean component and
two ripple harmonics.  The health summary separates mean torque, RMS ripple,
dominant harmonic order, and the variance budget of the largest harmonics.
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

from radia_mcp.radia_ngsolve.solve import torque_angle_sweep_health_summary  # noqa: E402
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_torque_waveform_health_summary.json"
SAMPLES = 720


def _waveform(mean_torque_Nm: float, ripple6_Nm: float, ripple12_Nm: float) -> list[float]:
    values = []
    for idx in range(SAMPLES):
        theta = 2.0 * math.pi * idx / SAMPLES
        values.append(
            mean_torque_Nm
            + ripple6_Nm * math.cos(6 * theta)
            + ripple12_Nm * math.sin(12 * theta)
        )
    return values


def _assert_close(actual: float, expected: float, atol: float = 1.0e-12) -> float:
    error = abs(actual - expected)
    if error > atol:
        raise AssertionError(f"{actual!r} != {expected!r}; error={error!r}")
    return error


def build_summary() -> dict[str, object]:
    mean = 10.0
    ripple6 = 0.4
    ripple12 = 0.1
    torque = _waveform(mean, ripple6, ripple12)
    health = torque_angle_sweep_health_summary(
        torque,
        max_harmonic=18,
        max_ac_rms_over_mean=0.04,
        allowed_dominant_harmonics=[6],
        min_mean_abs_torque_Nm=9.0,
        top_harmonics=3,
    )
    top = {row["order"]: row for row in health["top_harmonic_rows"]}
    expected_ac_rms = math.sqrt((ripple6 * ripple6 + ripple12 * ripple12) / 2.0)
    expected_variance6 = ripple6 * ripple6 / (ripple6 * ripple6 + ripple12 * ripple12)
    expected_variance12 = ripple12 * ripple12 / (ripple6 * ripple6 + ripple12 * ripple12)
    errors = {
        "mean_error_Nm": _assert_close(health["mean_torque_Nm"], mean),
        "ac_rms_error_Nm": _assert_close(health["ac_rms_torque_Nm"], expected_ac_rms),
        "sixth_amplitude_error_Nm": _assert_close(top[6]["amplitude_Nm"], ripple6),
        "twelfth_amplitude_error_Nm": _assert_close(top[12]["amplitude_Nm"], ripple12),
        "sixth_variance_fraction_error": _assert_close(top[6]["ac_variance_fraction"], expected_variance6),
        "twelfth_variance_fraction_error": _assert_close(top[12]["ac_variance_fraction"], expected_variance12),
    }
    checks = {
        "samples": SAMPLES,
        "status": health["status"],
        "mean_torque_Nm": health["mean_torque_Nm"],
        "ac_rms_torque_Nm": health["ac_rms_torque_Nm"],
        "ac_rms_over_mean": health["ac_rms_over_mean"],
        "dominant_harmonic": health["dominant_harmonic"],
        "dominant_harmonic_amplitude_Nm": health["dominant_harmonic_amplitude_Nm"],
        "sixth_variance_fraction": top[6]["ac_variance_fraction"],
        "twelfth_variance_fraction": top[12]["ac_variance_fraction"],
        "max_abs_error": max(errors.values()),
    }

    assert checks["status"] == "ok"
    assert checks["dominant_harmonic"] == 6
    assert checks["max_abs_error"] < 1.0e-12

    return {
        "kind": "torque_waveform_health_validation",
        "validation_class": True,
        "learning_theme": (
            "torque-angle tables should separate mean torque, RMS ripple, "
            "dominant harmonic order, and harmonic variance budget"
        ),
        "checks": checks,
        "errors": errors,
        "health": health,
    }


def _json_clean(value):
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, tuple):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = add_result_metadata(_json_clean(build_summary()), __file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[torque waveform health]")
    print(
        f"  samples={checks['samples']} status={checks['status']} "
        f"mean={checks['mean_torque_Nm']:.12g} ac_rms={checks['ac_rms_torque_Nm']:.12g}"
    )
    print(
        f"  dominant_harmonic={checks['dominant_harmonic']} "
        f"amplitude={checks['dominant_harmonic_amplitude_Nm']:.12g}"
    )
    print(
        f"  variance fractions: h6={checks['sixth_variance_fraction']:.12g} "
        f"h12={checks['twelfth_variance_fraction']:.12g}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
