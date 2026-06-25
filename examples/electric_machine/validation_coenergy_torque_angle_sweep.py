"""Validation-class coenergy torque from a rotor-angle sweep.

For fixed currents, virtual work gives

    T(theta) = dW'(theta) / dtheta

where ``W'`` is magnetic coenergy.  This example uses an analytic periodic
coenergy table to validate the same finite-difference post-processing used for
motor angle sweeps.

Run:

    python examples/electric_machine/validation_coenergy_torque_angle_sweep.py
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

from radia_mcp.radia_ngsolve.force import coenergy_torque_summary  # noqa: E402


OUT_JSON = HERE / "validation_coenergy_torque_angle_sweep_summary.json"

SAMPLES = 1440
HARMONIC_ORDER = 3
COENERGY_OFFSET_J = 2.0
COENERGY_AMPLITUDE_J = 0.75


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-12) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    angles = [2.0 * math.pi * index / SAMPLES for index in range(SAMPLES)]
    coenergy = [
        COENERGY_OFFSET_J - COENERGY_AMPLITUDE_J * math.cos(HARMONIC_ORDER * angle)
        for angle in angles
    ]
    exact_torque = [
        COENERGY_AMPLITUDE_J * HARMONIC_ORDER * math.sin(HARMONIC_ORDER * angle)
        for angle in angles
    ]
    summary = coenergy_torque_summary(angles, coenergy, periodic=True)
    errors = [
        row["torque_Nm"] - exact
        for row, exact in zip(summary["rows"], exact_torque)
    ]
    checks = {
        "samples": SAMPLES,
        "harmonic_order": HARMONIC_ORDER,
        "coenergy_amplitude_J": COENERGY_AMPLITUDE_J,
        "expected_torque_peak_abs_Nm": COENERGY_AMPLITUDE_J * HARMONIC_ORDER,
        "computed_torque_peak_abs_Nm": summary["torque_peak_abs_Nm"],
        "torque_peak_rel_error": abs(
            summary["torque_peak_abs_Nm"] - COENERGY_AMPLITUDE_J * HARMONIC_ORDER
        ) / (COENERGY_AMPLITUDE_J * HARMONIC_ORDER),
        "max_torque_abs_error_Nm": max(abs(value) for value in errors),
        "torque_mean_Nm": summary["torque_mean_Nm"],
    }

    if checks["torque_peak_rel_error"] > 4.0e-5:
        raise AssertionError("coenergy torque peak is outside tolerance")
    if checks["max_torque_abs_error_Nm"] > 7.0e-5:
        raise AssertionError("coenergy torque waveform is outside tolerance")
    _assert_close(checks["torque_mean_Nm"], 0.0, atol=1.0e-13)

    sample_indices = [0, SAMPLES // 12, SAMPLES // 6, SAMPLES // 4, SAMPLES // 3, SAMPLES // 2]
    compact_summary = {key: value for key, value in summary.items() if key != "rows"}
    return {
        "kind": "coenergy_torque_angle_sweep",
        "validation_class": True,
        "force_learning": "fixed-current virtual work: torque is d(coenergy)/d(theta)",
        "coenergy_model": {
            "offset_J": COENERGY_OFFSET_J,
            "amplitude_J": COENERGY_AMPLITUDE_J,
            "harmonic_order": HARMONIC_ORDER,
            "formula": "Wprime(theta) = offset - amplitude*cos(harmonic*theta)",
        },
        "checks": checks,
        "summary": compact_summary,
        "sample_rows": [summary["rows"][index] for index in sample_indices],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[Coenergy torque angle sweep]")
    print(f"  expected_torque_peak_abs_Nm: {checks['expected_torque_peak_abs_Nm']:.12g}")
    print(f"  computed_torque_peak_abs_Nm: {checks['computed_torque_peak_abs_Nm']:.12g}")
    print(f"  torque_peak_rel_error: {checks['torque_peak_rel_error']:.3e}")
    print(f"  max_torque_abs_error_Nm: {checks['max_torque_abs_error_Nm']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
