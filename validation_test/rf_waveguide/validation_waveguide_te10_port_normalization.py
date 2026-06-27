"""Validation-class TE10 rectangular-waveguide port normalization.

The example is intentionally closed-form and public-safe: it checks the power
normalization used by a TE10 port without requiring a commercial RF solve.
"""

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from result_metadata import add_result_metadata  # noqa: E402

from radia_mcp.radia_ngsolve.waveguide import rectangular_waveguide_te10_port_normalization


WR90_A = 0.02286
WR90_B = 0.01016


def main():
    cases = [
        {"label": "wr90_low_band_1w", "frequency": 8.2e9, "power_w": 1.0},
        {"label": "wr90_mid_band_1w", "frequency": 10.0e9, "power_w": 1.0},
        {"label": "wr90_high_band_1w", "frequency": 12.4e9, "power_w": 1.0},
        {"label": "wr90_mid_band_4w", "frequency": 10.0e9, "power_w": 4.0},
    ]
    rows = [
        {
            "label": case["label"],
            "normalization": rectangular_waveguide_te10_port_normalization(
                case["frequency"], WR90_A, WR90_B, power_w=case["power_w"]
            ),
        }
        for case in cases
    ]

    one = rows[1]["normalization"]
    four = rows[3]["normalization"]
    max_power_abs_error = max(row["normalization"]["poynting_abs_error_W"] for row in rows)
    max_hz_ratio_error = max(
        abs(
            row["normalization"]["H_z_over_H_x_peak"]
            - (math.pi / WR90_A) / row["normalization"]["beta"]
        )
        for row in rows
    )
    power_scaling_e_error = abs(
        four["E_y_peak_V_per_m"] / one["E_y_peak_V_per_m"] - math.sqrt(4.0)
    )
    power_scaling_h_error = abs(
        four["H_x_peak_A_per_m"] / one["H_x_peak_A_per_m"] - math.sqrt(4.0)
    )
    one_w_rows = [row["normalization"] for row in rows[:3]]
    monotone = (
        one_w_rows[0]["Z_TE_ohm"] > one_w_rows[1]["Z_TE_ohm"] > one_w_rows[2]["Z_TE_ohm"]
        and one_w_rows[0]["E_y_peak_V_per_m"]
        > one_w_rows[1]["E_y_peak_V_per_m"]
        > one_w_rows[2]["E_y_peak_V_per_m"]
        and one_w_rows[0]["v_group"] < one_w_rows[1]["v_group"] < one_w_rows[2]["v_group"]
    )

    checks = {
        "max_power_abs_error_W": max_power_abs_error,
        "max_hz_ratio_error": max_hz_ratio_error,
        "power_scaling_e_error": power_scaling_e_error,
        "power_scaling_h_error": power_scaling_h_error,
        "one_w_frequency_trends_ok": monotone,
        "midband_E_y_peak_V_per_m": one["E_y_peak_V_per_m"],
        "midband_H_x_peak_A_per_m": one["H_x_peak_A_per_m"],
        "midband_H_z_wall_peak_A_per_m": one["H_z_wall_peak_A_per_m"],
        "midband_Z_TE_ohm": one["Z_TE_ohm"],
    }
    assert max_power_abs_error < 1e-12
    assert max_hz_ratio_error < 1e-12
    assert power_scaling_e_error < 1e-12
    assert power_scaling_h_error < 1e-12
    assert monotone

    summary = {
        "kind": "waveguide_te10_port_normalization_validation",
        "validation_class": True,
        "guide": {"name": "WR-90", "width_a_m": WR90_A, "height_b_m": WR90_B},
        "rows": rows,
        "checks": checks,
    }
    out = Path(__file__).with_name("validation_waveguide_te10_port_normalization_summary.json")
    summary = add_result_metadata(summary, __file__)
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[TE10 port normalization]")
    for row in rows:
        n = row["normalization"]
        print(
            f"  {row['label']}: f={n['frequency']/1e9:.3f} GHz P={n['power_w']:.3g} W "
            f"E0={n['E_y_peak_V_per_m']:.6f} V/m Hx0={n['H_x_peak_A_per_m']:.6f} A/m "
            f"P_err={n['poynting_abs_error_W']:.3e} W"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    print(f"[OK] wrote {out}")


if __name__ == "__main__":
    main()
