"""Validation-class laminated-steel complex-permeability frequency sweep.

This is an example/validation run, not a pytest test.  It turns the compact
closed-form lamination model into a readable table that students can inspect:

* low-frequency limit -> fill-weighted static permeability
* thicker sheets or higher permeability enter flux exclusion earlier
* deep-skin high-frequency magnitude follows the expected 1/sqrt(f) law
* the imaginary part is negative under the solver's loss convention

Run:

    python examples/electric_machine/validation_lamination_mu_eff_sweep.py
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.solve import MU0, laminated_mu_eff  # noqa: E402


OUT_JSON = HERE / "validation_lamination_mu_eff_sweep_summary.json"
FREQUENCIES_HZ = (0.0, 1.0, 10.0, 50.0, 100.0, 1.0e3, 1.0e4, 1.0e5, 1.0e6)

CASES = [
    {
        "label": "mild_0p35mm_95fill",
        "mu_r": 2000.0,
        "sigma": 2.0e6,
        "d_lam": 0.35e-3,
        "fill": 0.95,
    },
    {
        "label": "thick_1p00mm_95fill",
        "mu_r": 2000.0,
        "sigma": 2.0e6,
        "d_lam": 1.00e-3,
        "fill": 0.95,
    },
    {
        "label": "high_mu_0p35mm_875fill",
        "mu_r": 5000.0,
        "sigma": 4.0e6,
        "d_lam": 0.35e-3,
        "fill": 0.875,
    },
]


def _omega(freq_hz: float) -> float:
    return 2.0 * math.pi * freq_hz


def _b_parameter(case: dict, freq_hz: float) -> complex:
    if freq_hz == 0.0:
        return 0.0j
    return (case["d_lam"] / 2.0) * cmath.sqrt(
        1j * _omega(freq_hz) * MU0 * case["mu_r"] * case["sigma"]
    )


def _static_mu_rel(case: dict) -> float:
    return case["fill"] * case["mu_r"] + (1.0 - case["fill"])


def _mu_rel(case: dict, freq_hz: float) -> complex:
    return laminated_mu_eff(
        case["mu_r"],
        case["sigma"],
        _omega(freq_hz),
        case["d_lam"],
        fill=case["fill"],
    ) / MU0


def _logspace(start: float, stop: float, count: int) -> list[float]:
    a = math.log10(start)
    b = math.log10(stop)
    return [10.0 ** (a + (b - a) * i / (count - 1)) for i in range(count)]


def _half_magnitude_frequency(case: dict, threshold: float = 0.5) -> float | None:
    static = _static_mu_rel(case)
    prev_f = 1.0e-3
    prev = abs(_mu_rel(case, prev_f)) / static
    for freq in _logspace(1.0e-3, 1.0e7, 2001):
        ratio = abs(_mu_rel(case, freq)) / static
        if ratio <= threshold:
            if ratio == prev:
                return freq
            # Log-frequency interpolation is enough for a reporting metric.
            x0 = math.log(prev_f)
            x1 = math.log(freq)
            frac = (threshold - prev) / (ratio - prev)
            return math.exp(x0 + frac * (x1 - x0))
        prev_f, prev = freq, ratio
    return None


def _case_record(case: dict) -> dict:
    static = _static_mu_rel(case)
    rows = []
    for freq in FREQUENCIES_HZ:
        mu = _mu_rel(case, freq)
        ratio = abs(mu) / static
        rows.append({
            "frequency_hz": freq,
            "b_abs": abs(_b_parameter(case, freq)),
            "mu_rel_real": mu.real,
            "mu_rel_imag": mu.imag,
            "mu_rel_abs": abs(mu),
            "abs_over_static": ratio,
            "loss_tangent_like": (-mu.imag / mu.real) if mu.real != 0.0 else None,
        })

    high_a = abs(_mu_rel(case, 1.0e5))
    high_b = abs(_mu_rel(case, 1.0e6))
    return {
        "label": case["label"],
        "parameters": case,
        "static_mu_rel": static,
        "half_magnitude_frequency_hz": _half_magnitude_frequency(case),
        "high_frequency_abs_ratio_100k_over_1M": high_a / high_b,
        "expected_high_frequency_ratio": math.sqrt(10.0),
        "rows": rows,
    }


def _validate(records: list[dict]) -> dict:
    by_label = {rec["label"]: rec for rec in records}
    static_errors = []
    imag_nonzero = []
    monotone_violations = []
    hf_errors = []
    for rec in records:
        rows = rec["rows"]
        static_errors.append(abs(rows[0]["mu_rel_real"] - rec["static_mu_rel"]))
        imag_nonzero.extend(row["mu_rel_imag"] for row in rows[1:])
        for left, right in zip(rows, rows[1:]):
            if right["abs_over_static"] > left["abs_over_static"] + 1.0e-12:
                monotone_violations.append((rec["label"], left["frequency_hz"], right["frequency_hz"]))
        hf = rec["high_frequency_abs_ratio_100k_over_1M"]
        hf_errors.append(abs(hf / rec["expected_high_frequency_ratio"] - 1.0))

    thin_half = by_label["mild_0p35mm_95fill"]["half_magnitude_frequency_hz"]
    thick_half = by_label["thick_1p00mm_95fill"]["half_magnitude_frequency_hz"]
    high_mu_half = by_label["high_mu_0p35mm_875fill"]["half_magnitude_frequency_hz"]
    checks = {
        "max_static_abs_error": max(static_errors),
        "max_high_frequency_ratio_rel_error": max(hf_errors),
        "monotone_violations": monotone_violations,
        "thick_half_frequency_below_thin": thick_half < thin_half,
        "high_mu_half_frequency_below_thin": high_mu_half < thin_half,
        "all_nonzero_imaginary_parts_negative": all(v < 0.0 for v in imag_nonzero),
    }
    assert checks["max_static_abs_error"] < 1.0e-12
    assert checks["max_high_frequency_ratio_rel_error"] < 0.01
    assert not checks["monotone_violations"]
    assert checks["thick_half_frequency_below_thin"]
    assert checks["high_mu_half_frequency_below_thin"]
    assert checks["all_nonzero_imaginary_parts_negative"]
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    records = [_case_record(case) for case in CASES]
    checks = _validate(records)
    summary = {
        "kind": "lamination_mu_eff_frequency_sweep_validation",
        "validation_class": True,
        "frequencies_hz": list(FREQUENCIES_HZ),
        "cases": records,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    for rec in records:
        print(
            f"[{rec['label']}] static mu/mu0={rec['static_mu_rel']:.6g}, "
            f"half-|mu| frequency ~= {rec['half_magnitude_frequency_hz']:.6g} Hz"
        )
        for row in rec["rows"]:
            print(
                f"  f={row['frequency_hz']:9.0f} Hz  |b|={row['b_abs']:8.4g}  "
                f"mu/mu0={row['mu_rel_real']:10.4g}{row['mu_rel_imag']:+10.4g}j  "
                f"|mu|/static={row['abs_over_static']:.6f}"
            )
        print(
            "  |mu|(100 kHz)/|mu|(1 MHz) = "
            f"{rec['high_frequency_abs_ratio_100k_over_1M']:.6f} "
            f"(sqrt(10)={math.sqrt(10.0):.6f})"
        )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
