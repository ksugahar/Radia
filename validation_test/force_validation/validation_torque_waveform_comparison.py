"""Validation-class periodic torque waveform comparison.

This example uses an analytic before/after torque-angle table to check that
the comparison summary separates mean torque drift, sample-wise error, and
harmonic ripple changes.

Run:

    python validation_test/force_validation/validation_torque_waveform_comparison.py
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

from radia_mcp.radia_ngsolve.solve import torque_angle_sweep_comparison_summary  # noqa: E402
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_torque_waveform_comparison_summary.json"
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


def _by_order(rows: list[dict[str, float]]) -> dict[int, dict[str, float]]:
    return {int(row["order"]): row for row in rows}


def _assert_close(actual: float, expected: float, atol: float = 1.0e-12) -> float:
    error = abs(actual - expected)
    if error > atol:
        raise AssertionError(f"{actual!r} != {expected!r}; error={error!r}")
    return error


def build_summary() -> dict[str, object]:
    reference = _waveform(mean_torque_Nm=10.0, ripple6_Nm=0.4, ripple12_Nm=0.1)
    candidate = _waveform(mean_torque_Nm=10.2, ripple6_Nm=0.3, ripple12_Nm=0.12)
    comparison = torque_angle_sweep_comparison_summary(
        reference,
        candidate,
        max_harmonic=18,
        reference_label="analytic_reference",
        candidate_label="analytic_candidate",
    )
    rows = _by_order(comparison["harmonic_delta_rows"])

    expected_sample_delta_rms = math.sqrt(0.2 * 0.2 + (0.1 * 0.1 + 0.02 * 0.02) / 2.0)
    expected_delta_ac_rms = math.sqrt((0.1 * 0.1 + 0.02 * 0.02) / 2.0)
    errors = {
        "mean_delta_error_Nm": _assert_close(comparison["mean_delta_Nm"], 0.2),
        "sample_delta_rms_error_Nm": _assert_close(
            comparison["sample_delta_rms_Nm"],
            expected_sample_delta_rms,
        ),
        "delta_ac_rms_error_Nm": _assert_close(
            comparison["difference_summary"]["ac_rms_torque_Nm"],
            expected_delta_ac_rms,
        ),
        "sixth_harmonic_delta_error_Nm": _assert_close(rows[6]["amplitude_delta_Nm"], -0.1),
        "twelfth_harmonic_delta_error_Nm": _assert_close(rows[12]["amplitude_delta_Nm"], 0.02),
    }

    if comparison["dominant_harmonic_changed"]:
        raise AssertionError("dominant harmonic should remain unchanged")
    if comparison["worst_harmonic_order"] != 6:
        raise AssertionError("largest ripple amplitude change should be harmonic order 6")

    checks = {
        "samples": SAMPLES,
        "passed": True,
        "max_abs_error_Nm": max(errors.values()),
        "mean_delta_Nm": comparison["mean_delta_Nm"],
        "sample_delta_rms_Nm": comparison["sample_delta_rms_Nm"],
        "expected_sample_delta_rms_Nm": expected_sample_delta_rms,
        "delta_ac_rms_Nm": comparison["difference_summary"]["ac_rms_torque_Nm"],
        "expected_delta_ac_rms_Nm": expected_delta_ac_rms,
        "dominant_harmonic_changed": comparison["dominant_harmonic_changed"],
        "worst_harmonic_order": comparison["worst_harmonic_order"],
        "sixth_harmonic_amplitude_delta_Nm": rows[6]["amplitude_delta_Nm"],
        "twelfth_harmonic_amplitude_delta_Nm": rows[12]["amplitude_delta_Nm"],
    }

    return {
        "kind": "torque_waveform_comparison_validation",
        "validation_class": True,
        "learning_theme": (
            "periodic torque-angle comparisons should separate mean drift, "
            "sample-wise error, and harmonic ripple deltas"
        ),
        "checks": checks,
        "errors": errors,
        "comparison": comparison,
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
    print("[torque waveform comparison]")
    print(f"  samples={checks['samples']}, max_abs_error={checks['max_abs_error_Nm']:.3e}")
    print(f"  mean_delta_Nm={checks['mean_delta_Nm']:.12g}")
    print(f"  sample_delta_rms_Nm={checks['sample_delta_rms_Nm']:.12g}")
    print(f"  worst_harmonic_order={checks['worst_harmonic_order']}")
    print(f"  sixth_harmonic_amplitude_delta_Nm={checks['sixth_harmonic_amplitude_delta_Nm']:.12g}")
    print(f"  twelfth_harmonic_amplitude_delta_Nm={checks['twelfth_harmonic_amplitude_delta_Nm']:.12g}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
