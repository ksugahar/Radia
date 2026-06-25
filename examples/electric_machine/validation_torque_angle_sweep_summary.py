"""Validation-class torque-angle sweep summary.

Torque tables from rotor-position sweeps are easiest to compare when the
post-processing reports the same compact metrics: mean torque, peak-to-peak
ripple, AC RMS, harmonic amplitudes, and the dominant ripple order.  This
example validates that summary on an analytic torque waveform.

Run:

    python examples/electric_machine/validation_torque_angle_sweep_summary.py
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

from radia_mcp.radia_ngsolve.solve import torque_angle_sweep_summary  # noqa: E402


OUT_JSON = HERE / "validation_torque_angle_sweep_summary.json"

SAMPLES = 720
MEAN_TORQUE_NM = 10.0
RIPPLE6_NM = 0.4
RIPPLE12_NM = 0.1


def _torque_waveform():
    values = []
    for idx in range(SAMPLES):
        theta = 2.0 * math.pi * idx / SAMPLES
        values.append(
            MEAN_TORQUE_NM
            + RIPPLE6_NM * math.cos(6 * theta)
            + RIPPLE12_NM * math.sin(12 * theta)
        )
    return values


def _by_order(summary):
    return {row["order"]: row for row in summary["harmonic_rows"]}


def _assert_close(actual, expected, atol=1.0e-12):
    if abs(actual - expected) > atol:
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    torque = _torque_waveform()
    summary = torque_angle_sweep_summary(torque, max_harmonic=18)
    rows = _by_order(summary)
    expected_rms = math.sqrt((RIPPLE6_NM * RIPPLE6_NM + RIPPLE12_NM * RIPPLE12_NM) / 2.0)

    _assert_close(summary["mean_torque_Nm"], MEAN_TORQUE_NM)
    _assert_close(summary["ac_rms_torque_Nm"], expected_rms)
    _assert_close(rows[6]["amplitude_Nm"], RIPPLE6_NM)
    _assert_close(rows[12]["amplitude_Nm"], RIPPLE12_NM)
    if summary["dominant_harmonic"] != 6:
        raise AssertionError("dominant torque ripple harmonic drifted")

    return {
        "kind": "torque_angle_sweep_summary_validation",
        "validation_class": True,
        "force_learning": "periodic torque-angle tables reduce to mean, RMS, peak-to-peak, and harmonic amplitude metrics",
        "waveform": {
            "samples": SAMPLES,
            "mean_torque_Nm": MEAN_TORQUE_NM,
            "ripple6_amplitude_Nm": RIPPLE6_NM,
            "ripple12_amplitude_Nm": RIPPLE12_NM,
        },
        "checks": {
            "mean_torque_Nm": summary["mean_torque_Nm"],
            "ac_rms_torque_Nm": summary["ac_rms_torque_Nm"],
            "expected_ac_rms_torque_Nm": expected_rms,
            "dominant_harmonic": summary["dominant_harmonic"],
            "dominant_harmonic_amplitude_Nm": summary["dominant_harmonic_amplitude_Nm"],
            "sixth_harmonic_amplitude_Nm": rows[6]["amplitude_Nm"],
            "twelfth_harmonic_amplitude_Nm": rows[12]["amplitude_Nm"],
            "sixth_harmonic_over_mean": rows[6]["amplitude_over_mean"],
            "twelfth_harmonic_over_mean": rows[12]["amplitude_over_mean"],
        },
        "summary": summary,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[Torque-angle sweep summary]")
    print(f"  mean_torque_Nm: {checks['mean_torque_Nm']:.12g}")
    print(f"  ac_rms_torque_Nm: {checks['ac_rms_torque_Nm']:.12g}")
    print(f"  dominant_harmonic: {checks['dominant_harmonic']}")
    print(f"  sixth_harmonic_amplitude_Nm: {checks['sixth_harmonic_amplitude_Nm']:.12g}")
    print(f"  twelfth_harmonic_amplitude_Nm: {checks['twelfth_harmonic_amplitude_Nm']:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
