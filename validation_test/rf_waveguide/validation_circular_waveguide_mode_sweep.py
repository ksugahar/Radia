"""Validation-class circular-waveguide cutoff and band sweep.

This is an example/validation run rather than a pytest test. It complements
the rectangular-guide table with the Bessel-zero circular guide:

* dominant TE11 from the first zero of J1';
* TM01 and TE21 as the next cutoffs;
* TE0n / TM1n degeneracy;
* single-mode band classification with angular degeneracy counted explicitly.

Run:

    python validation_test/rf_waveguide/validation_circular_waveguide_mode_sweep.py
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

from radia_mcp.radia_ngsolve.waveguide import (  # noqa: E402
    C0,
    circular_waveguide_band_summary,
    circular_waveguide_cutoff,
    circular_waveguide_mode_table,
    waveguide_dispersion,
    waveguide_evanescent_attenuation,
)


OUT_JSON = HERE / "validation_circular_waveguide_mode_sweep_summary.json"
RADIUS_M = 12.7e-3
FREQUENCIES_HZ = (6.0e9, 7.5e9, 10.0e9, 12.0e9, 15.0e9)


def _row_mode(row: dict) -> dict:
    return {
        **row,
        "cutoff_frequency_ghz": row["cutoff_frequency"] / 1.0e9,
    }


def _json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def build_summary() -> dict:
    table = [_row_mode(row) for row in circular_waveguide_mode_table(RADIUS_M, max_m=3, max_n=2)]
    bands = []
    for freq in FREQUENCIES_HZ:
        band = circular_waveguide_band_summary(RADIUS_M, freq, max_m=3, max_n=2)
        propagating_names = [row["mode"] for row in band["propagating_modes"]]
        row = {
            "frequency_hz": freq,
            "frequency_ghz": freq / 1.0e9,
            "below_dominant_cutoff": band["below_dominant_cutoff"],
            "single_mode": band["single_mode"],
            "n_propagating_rows": band["n_propagating_rows"],
            "n_propagating_with_degeneracy": band["n_propagating_with_degeneracy"],
            "propagating_modes": propagating_names,
            "next_mode": band["next_mode"]["mode"] if band["next_mode"] else None,
        }
        if band["propagating_modes"]:
            dominant = band["propagating_modes"][0]
            disp = waveguide_dispersion(freq, dominant["cutoff_frequency"])
            row["dominant_beta_per_m"] = disp["beta"]
            row["dominant_guide_wavelength_m"] = disp["lambda_g"]
            row["dominant_group_velocity_m_per_s"] = disp["v_group"]
        else:
            fc = band["dominant_mode"]["cutoff_frequency"]
            row["dominant_evanescent_alpha_per_m"] = waveguide_evanescent_attenuation(freq, fc)
        bands.append(row)

    return {
        "kind": "circular_waveguide_mode_sweep_validation",
        "validation_class": True,
        "radius_m": RADIUS_M,
        "mode_table": table,
        "band_summaries": bands,
    }


def validate(summary: dict) -> dict:
    table = summary["mode_table"]
    by_mode = {row["mode"]: row for row in table}
    bands = summary["band_summaries"]
    first_modes = [row["mode"] for row in table[:5]]
    te11 = by_mode["TE11"]["cutoff_frequency"]
    tm01 = by_mode["TM01"]["cutoff_frequency"]
    te21 = by_mode["TE21"]["cutoff_frequency"]
    degeneracy_error = abs(by_mode["TE01"]["cutoff_frequency"] - by_mode["TM11"]["cutoff_frequency"])
    radius_scaling_error = abs(circular_waveguide_cutoff(2.0 * RADIUS_M, "TE", 1, 1) - 0.5 * te11)
    checks = {
        "first_five_modes": first_modes,
        "TE11_cutoff_hz": te11,
        "TM01_cutoff_hz": tm01,
        "TE21_cutoff_hz": te21,
        "TE11_below_TM01": te11 < tm01,
        "TM01_below_TE21": tm01 < te21,
        "TE01_TM11_degeneracy_abs_error_hz": degeneracy_error,
        "radius_scaling_abs_error_hz": radius_scaling_error,
        "below_cutoff_6ghz": bands[0]["below_dominant_cutoff"],
        "single_mode_7p5ghz": bands[1]["single_mode"],
        "multi_mode_10ghz": not bands[2]["single_mode"] and bands[2]["n_propagating_rows"] >= 2,
        "single_mode_degeneracy_count": bands[1]["n_propagating_with_degeneracy"],
    }
    assert first_modes == ["TE11", "TM01", "TE21", "TE01", "TM11"]
    assert checks["TE11_below_TM01"]
    assert checks["TM01_below_TE21"]
    assert degeneracy_error < 1.0e-3
    assert radius_scaling_error < 1.0e-6
    assert checks["below_cutoff_6ghz"]
    assert checks["single_mode_7p5ghz"]
    assert checks["multi_mode_10ghz"]
    assert checks["single_mode_degeneracy_count"] == 2
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    checks = validate(summary)
    summary["checks"] = checks
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary = add_result_metadata(summary, __file__)
    args.out.write_text(json.dumps(summary, indent=2, default=_json_default), encoding="utf-8")

    print("[Circular waveguide cutoff table]")
    for row in summary["mode_table"][:8]:
        print(
            f"  {row['mode']:<4}  fc={row['cutoff_frequency_ghz']:9.6f} GHz  "
            f"deg={row['angular_degeneracy']}"
        )
    print("[band sweep]")
    for row in summary["band_summaries"]:
        print(
            f"  f={row['frequency_ghz']:5.1f} GHz  "
            f"modes={row['propagating_modes']}  "
            f"single={row['single_mode']}  next={row['next_mode']}"
        )
    print(
        "[checks] "
        f"TE11={checks['TE11_cutoff_hz'] / 1e9:.6f} GHz, "
        f"TM01={checks['TM01_cutoff_hz'] / 1e9:.6f} GHz, "
        f"TE01/TM11 degeneracy error={checks['TE01_TM11_degeneracy_abs_error_hz']:.3e} Hz"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
