"""Validation-class rectangular-waveguide conductor-loss sweep.

The example uses a WR-90-like TE10 guide with finite-conductivity walls:

* frequency sweep: wall loss rises strongly near cutoff;
* length sweep: insertion loss is linear in length in dB;
* conductivity sweep: attenuation scales as 1/sqrt(sigma);
* S21 power balance is checked for a matched lossy line section.

It is intentionally an example/validation run, not a pytest test.

Run:

    python validation_test/rf_waveguide/validation_waveguide_conductor_loss.py
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
    rectangular_waveguide_cutoff,
    rectangular_waveguide_te10_conductor_loss,
)


OUT_JSON = HERE / "validation_waveguide_conductor_loss_summary.json"

WIDTH_A_M = 0.02286
HEIGHT_B_M = 0.01016
COPPER_SIGMA_S_PER_M = 5.8e7
BASE_LENGTH_M = 0.10
BASE_FREQUENCY_HZ = 10.0e9


def _row_subset(row: dict) -> dict:
    return {
        "frequency_hz": row["frequency"],
        "f_over_fc": row["frequency"] / row["fc"],
        "surface_resistance_ohm": row["surface_resistance_ohm"],
        "skin_depth_m": row["skin_depth_m"],
        "alpha_np_per_m": row["alpha_np_per_m"],
        "alpha_db_per_m": row["alpha_db_per_m"],
        "length_m": row.get("length_m"),
        "S21_mag": row.get("S21_mag"),
        "insertion_loss_db": row.get("insertion_loss_db"),
        "power_loss_fraction": row.get("power_loss_fraction"),
    }


def build_frequency_rows() -> list[dict]:
    fc = rectangular_waveguide_cutoff(WIDTH_A_M, HEIGHT_B_M, 1, 0)
    rows = []
    for ratio in (1.02, 1.05, 1.10, 1.25, 1.50, 2.00, 3.00):
        rows.append(_row_subset(rectangular_waveguide_te10_conductor_loss(
            ratio * fc, WIDTH_A_M, HEIGHT_B_M, COPPER_SIGMA_S_PER_M, length=BASE_LENGTH_M)))
    return rows


def build_length_rows() -> list[dict]:
    rows = []
    for length in (0.0, 0.025, 0.050, 0.100, 0.200, 0.500):
        rows.append(_row_subset(rectangular_waveguide_te10_conductor_loss(
            BASE_FREQUENCY_HZ, WIDTH_A_M, HEIGHT_B_M, COPPER_SIGMA_S_PER_M, length=length)))
    return rows


def build_conductivity_rows() -> list[dict]:
    rows = []
    for scale in (0.25, 1.0, 4.0, 16.0):
        sigma = COPPER_SIGMA_S_PER_M * scale
        row = rectangular_waveguide_te10_conductor_loss(
            BASE_FREQUENCY_HZ, WIDTH_A_M, HEIGHT_B_M, sigma, length=BASE_LENGTH_M)
        selected = _row_subset(row)
        selected["sigma_s_per_m"] = sigma
        selected["sigma_over_copper"] = scale
        selected["alpha_times_sqrt_sigma_scale"] = row["alpha_np_per_m"] * math.sqrt(scale)
        rows.append(selected)
    return rows


def validate(frequency_rows: list[dict], length_rows: list[dict],
             conductivity_rows: list[dict]) -> dict:
    base = rectangular_waveguide_te10_conductor_loss(
        BASE_FREQUENCY_HZ, WIDTH_A_M, HEIGHT_B_M, COPPER_SIGMA_S_PER_M, length=BASE_LENGTH_M)
    db_factor = 20.0 * math.log10(math.e)

    length_errors = []
    power_balance_errors = []
    for row in length_rows:
        expected_loss = base["alpha_db_per_m"] * row["length_m"]
        if expected_loss:
            length_errors.append(abs(row["insertion_loss_db"] - expected_loss) / expected_loss)
        power_balance_errors.append(abs(row["S21_mag"] ** 2 + row["power_loss_fraction"] - 1.0))

    sigma_scaled = [
        row["alpha_times_sqrt_sigma_scale"]
        for row in conductivity_rows
    ]
    sigma_errors = [
        abs(value - base["alpha_np_per_m"]) / base["alpha_np_per_m"]
        for value in sigma_scaled
    ]
    db_np_errors = [
        abs(row["alpha_db_per_m"] - db_factor * row["alpha_np_per_m"])
        / row["alpha_db_per_m"]
        for row in frequency_rows + length_rows + conductivity_rows
    ]

    checks = {
        "te10_cutoff_hz": rectangular_waveguide_cutoff(WIDTH_A_M, HEIGHT_B_M, 1, 0),
        "base_frequency_hz": BASE_FREQUENCY_HZ,
        "base_alpha_np_per_m": base["alpha_np_per_m"],
        "base_alpha_db_per_m": base["alpha_db_per_m"],
        "base_insertion_loss_db_0p1m": base["insertion_loss_db"],
        "near_cutoff_alpha_over_10ghz": frequency_rows[0]["alpha_np_per_m"] / base["alpha_np_per_m"],
        "max_length_linearity_rel_error": max(length_errors),
        "max_power_balance_abs_error": max(power_balance_errors),
        "max_sigma_scaling_rel_error": max(sigma_errors),
        "max_db_np_conversion_rel_error": max(db_np_errors),
    }

    assert checks["near_cutoff_alpha_over_10ghz"] > 1.0
    assert checks["max_length_linearity_rel_error"] < 1.0e-12
    assert checks["max_power_balance_abs_error"] < 1.0e-15
    assert checks["max_sigma_scaling_rel_error"] < 1.0e-14
    assert checks["max_db_np_conversion_rel_error"] < 1.0e-14
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    frequency_rows = build_frequency_rows()
    length_rows = build_length_rows()
    conductivity_rows = build_conductivity_rows()
    checks = validate(frequency_rows, length_rows, conductivity_rows)

    summary = {
        "kind": "waveguide_conductor_loss_validation",
        "validation_class": True,
        "guide": {
            "name": "WR-90-like",
            "width_a_m": WIDTH_A_M,
            "height_b_m": HEIGHT_B_M,
            "te10_cutoff_hz": checks["te10_cutoff_hz"],
        },
        "material": {
            "reference": "copper-like good conductor",
            "sigma_s_per_m": COPPER_SIGMA_S_PER_M,
        },
        "frequency_rows": frequency_rows,
        "length_rows": length_rows,
        "conductivity_rows": conductivity_rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary = add_result_metadata(summary, __file__)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[waveguide conductor loss]")
    print(f"  TE10 cutoff = {checks['te10_cutoff_hz'] / 1e9:.6f} GHz")
    print(
        f"  10 GHz alpha = {checks['base_alpha_np_per_m']:.9e} Np/m "
        f"= {checks['base_alpha_db_per_m']:.9e} dB/m"
    )
    print(
        f"  0.10 m insertion loss = {checks['base_insertion_loss_db_0p1m']:.9e} dB"
    )
    print(
        "  near-cutoff alpha / 10 GHz alpha = "
        f"{checks['near_cutoff_alpha_over_10ghz']:.6f}"
    )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
