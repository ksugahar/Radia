"""Validation-class rectangular-waveguide mode-table sweep.

This is an example/validation run, not a pytest test.  It keeps a compact
analytic checklist for RF port work:

* rectangular TE/TM cutoff modes and their degeneracies
* below-cutoff, single-mode, and multi-mode frequency regions
* guide wavelength / group velocity for propagating TE10
* evanescent attenuation for a below-cutoff drive

Run:

    python examples/rf_waveguide/validation_waveguide_mode_table_sweep.py
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

from radia_mcp.radia_ngsolve.waveguide import (  # noqa: E402
    C0,
    rectangular_waveguide_band_summary,
    rectangular_waveguide_cutoff,
    rectangular_waveguide_mode_table,
    waveguide_dispersion,
    waveguide_evanescent_attenuation,
)


OUT_JSON = HERE / "validation_waveguide_mode_table_sweep_summary.json"
WIDTH_A = 0.02286
HEIGHT_B = 0.01016
FREQUENCIES_HZ = [5.0e9, 8.2e9, 10.0e9, 12.4e9, 13.2e9, 15.0e9, 18.0e9]


def _mode_record(row: dict) -> dict:
    return {
        "mode": row["mode"],
        "family": row["family"],
        "m": row["m"],
        "n": row["n"],
        "cutoff_frequency_hz": row["cutoff_frequency"],
        "cutoff_frequency_ghz": row["cutoff_frequency"] / 1.0e9,
        "cutoff_wavenumber_per_m": row["cutoff_wavenumber"],
    }


def _band_record(frequency: float) -> dict:
    summary = rectangular_waveguide_band_summary(WIDTH_A, HEIGHT_B, frequency, max_m=3, max_n=3)
    propagating = summary["propagating_modes"]
    next_mode = summary["next_mode"]
    row = {
        "frequency_hz": frequency,
        "frequency_ghz": frequency / 1.0e9,
        "below_dominant_cutoff": summary["below_dominant_cutoff"],
        "single_mode": summary["single_mode"],
        "n_propagating": summary["n_propagating"],
        "propagating_modes": [_mode_record(mode) for mode in propagating],
        "next_mode": None if next_mode is None else _mode_record(next_mode),
    }
    if propagating:
        first = propagating[0]
        disp = waveguide_dispersion(frequency, first["cutoff_frequency"])
        row["dominant_mode_propagation"] = {
            "mode": first["mode"],
            "beta_per_m": disp["beta"],
            "guide_wavelength_m": disp["lambda_g"],
            "v_group_over_c": disp["v_group"] / C0,
            "v_phase_over_c": disp["v_phase"] / C0,
        }
    else:
        fc10 = rectangular_waveguide_cutoff(WIDTH_A, HEIGHT_B, 1, 0)
        row["te10_evanescent_alpha_per_m"] = waveguide_evanescent_attenuation(frequency, fc10)
    return row


def _validate(table: list[dict], rows: list[dict]) -> dict:
    first_modes = [row["mode"] for row in table[:5]]
    by_mode = {row["mode"]: row for row in table}
    by_freq = {row["frequency_hz"]: row for row in rows}

    checks = {
        "first_five_modes": first_modes,
        "te10_cutoff_hz": by_mode["TE10"]["cutoff_frequency_hz"],
        "te20_over_te10": by_mode["TE20"]["cutoff_frequency_hz"] / by_mode["TE10"]["cutoff_frequency_hz"],
        "te01_cutoff_hz": by_mode["TE01"]["cutoff_frequency_hz"],
        "te11_tm11_degeneracy_abs_hz": abs(
            by_mode["TE11"]["cutoff_frequency_hz"] - by_mode["TM11"]["cutoff_frequency_hz"]
        ),
        "single_mode_upper_cutoff_hz": by_mode["TE20"]["cutoff_frequency_hz"],
        "xband_lower_single_mode": by_freq[8.2e9]["single_mode"],
        "xband_center_single_mode": by_freq[10.0e9]["single_mode"],
        "xband_upper_single_mode": by_freq[12.4e9]["single_mode"],
        "five_ghz_propagating_count": by_freq[5.0e9]["n_propagating"],
        "thirteen_two_ghz_modes": [mode["mode"] for mode in by_freq[13.2e9]["propagating_modes"]],
        "fifteen_ghz_modes": [mode["mode"] for mode in by_freq[15.0e9]["propagating_modes"]],
        "eighteen_ghz_modes": [mode["mode"] for mode in by_freq[18.0e9]["propagating_modes"]],
        "ten_ghz_te10_guide_wavelength_m": by_freq[10.0e9]["dominant_mode_propagation"]["guide_wavelength_m"],
        "five_ghz_te10_alpha_per_m": by_freq[5.0e9]["te10_evanescent_alpha_per_m"],
    }

    assert first_modes == ["TE10", "TE20", "TE01", "TE11", "TM11"]
    assert math.isclose(checks["te20_over_te10"], 2.0, rel_tol=1.0e-12)
    assert checks["te11_tm11_degeneracy_abs_hz"] < 1.0e-6
    assert checks["xband_lower_single_mode"]
    assert checks["xband_center_single_mode"]
    assert checks["xband_upper_single_mode"]
    assert checks["five_ghz_propagating_count"] == 0
    assert checks["thirteen_two_ghz_modes"] == ["TE10", "TE20"]
    assert checks["fifteen_ghz_modes"] == ["TE10", "TE20", "TE01"]
    assert checks["eighteen_ghz_modes"] == ["TE10", "TE20", "TE01", "TE11", "TM11"]
    assert checks["ten_ghz_te10_guide_wavelength_m"] > C0 / 10.0e9
    assert checks["five_ghz_te10_alpha_per_m"] > 0.0
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    table = [_mode_record(row) for row in rectangular_waveguide_mode_table(WIDTH_A, HEIGHT_B, max_m=3, max_n=3)]
    rows = [_band_record(frequency) for frequency in FREQUENCIES_HZ]
    checks = _validate(table, rows)
    summary = {
        "kind": "waveguide_mode_table_sweep_validation",
        "validation_class": True,
        "guide": {
            "name": "WR-90-like",
            "width_a_m": WIDTH_A,
            "height_b_m": HEIGHT_B,
        },
        "mode_table": table,
        "frequency_rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[rectangular guide]")
    print(
        f"  TE10={checks['te10_cutoff_hz']/1e9:.6f} GHz, "
        f"TE20={checks['single_mode_upper_cutoff_hz']/1e9:.6f} GHz, "
        f"TE01={checks['te01_cutoff_hz']/1e9:.6f} GHz"
    )
    print("[band sweep]")
    for row in rows:
        modes = ",".join(mode["mode"] for mode in row["propagating_modes"]) or "none"
        next_mode = "none" if row["next_mode"] is None else row["next_mode"]["mode"]
        print(
            f"  f={row['frequency_ghz']:5.2f} GHz  "
            f"n={row['n_propagating']}  modes={modes:<24} next={next_mode}"
        )
    print(
        "[checks] "
        f"first={checks['first_five_modes']}, "
        f"10 GHz lambda_g={checks['ten_ghz_te10_guide_wavelength_m']:.9f} m, "
        f"5 GHz alpha={checks['five_ghz_te10_alpha_per_m']:.6f} 1/m"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
