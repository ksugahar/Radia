"""Validation-class air-gap Maxwell shear torque identities.

This example is the motor-torque counterpart of the air-gap pressure example:
radial and tangential air-gap flux components create tangential Maxwell shear
stress,

    tau = Br Bt / mu0,
    T = tau * (r * angle * L) * r.

Run:

    python validation_test/electric_machine/validation_air_gap_shear_torque.py
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

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    MU0,
    air_gap_shear_stress,
    air_gap_shear_torque,
    air_gap_shear_torque_summary,
    maxwell_traction_summary,
)


OUT_JSON = HERE / "validation_air_gap_shear_torque_summary.json"

BR_T = 0.8
BT_T = 0.1
RADIUS_M = 0.05
AXIAL_LENGTH_M = 0.1


CASES = (
    {"label": "full_machine_positive", "Bt_T": BT_T, "angle_rad": 2.0 * math.pi},
    {"label": "half_sector_positive", "Bt_T": BT_T, "angle_rad": math.pi},
    {"label": "full_machine_reverse", "Bt_T": -BT_T, "angle_rad": 2.0 * math.pi},
    {"label": "no_tangential_field", "Bt_T": 0.0, "angle_rad": 2.0 * math.pi},
)


def build_rows() -> list[dict]:
    rows = []
    for case in CASES:
        summary = air_gap_shear_torque_summary(
            BR_T,
            case["Bt_T"],
            RADIUS_M,
            axial_length_m=AXIAL_LENGTH_M,
            angle_rad=case["angle_rad"],
        )
        traction = maxwell_traction_summary(
            (BR_T, case["Bt_T"], 0.0),
            (1.0, 0.0, 0.0),
            area_m2=summary["surface_area_m2"],
        )
        rows.append({
            "label": case["label"],
            "summary": summary,
            "traction": traction,
            "shear_vs_traction_abs_error": abs(
                summary["shear_stress_Pa"] - traction["tangential_traction_Pa"][1]
            ),
            "torque_identity_abs_error": abs(
                summary["torque_Nm"]
                - summary["shear_stress_Pa"] * RADIUS_M * RADIUS_M * case["angle_rad"] * AXIAL_LENGTH_M
            ),
            "force_identity_abs_error": abs(
                summary["tangential_force_N"]
                - summary["shear_stress_Pa"] * summary["surface_area_m2"]
            ),
        })
    return rows


def _assert_close(actual: float, expected: float, rtol=1.0e-12, atol=1.0e-9) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def validate(rows: list[dict]) -> dict:
    by_label = {row["label"]: row for row in rows}
    full = by_label["full_machine_positive"]["summary"]["torque_Nm"]
    half = by_label["half_sector_positive"]["summary"]["torque_Nm"]
    reverse = by_label["full_machine_reverse"]["summary"]["torque_Nm"]
    zero = by_label["no_tangential_field"]["summary"]["torque_Nm"]
    shear = air_gap_shear_stress(BR_T, BT_T)
    checks = {
        "shear_stress_Pa": shear,
        "full_machine_torque_Nm": full,
        "half_sector_torque_Nm": half,
        "reverse_torque_Nm": reverse,
        "zero_tangential_torque_Nm": zero,
        "closed_form_full_torque_Nm": (
            BR_T * BT_T / MU0 * RADIUS_M * RADIUS_M * 2.0 * math.pi * AXIAL_LENGTH_M
        ),
        "max_shear_vs_traction_abs_error": max(row["shear_vs_traction_abs_error"] for row in rows),
        "max_torque_identity_abs_error": max(row["torque_identity_abs_error"] for row in rows),
        "max_force_identity_abs_error": max(row["force_identity_abs_error"] for row in rows),
        "half_over_full": half / full if full else math.inf,
        "reverse_over_full": reverse / full if full else math.inf,
    }

    _assert_close(checks["shear_stress_Pa"], BR_T * BT_T / MU0)
    _assert_close(checks["full_machine_torque_Nm"], 100.0)
    _assert_close(checks["closed_form_full_torque_Nm"], 100.0)
    _assert_close(checks["half_over_full"], 0.5)
    _assert_close(checks["reverse_over_full"], -1.0)
    _assert_close(checks["zero_tangential_torque_Nm"], 0.0)
    _assert_close(checks["max_shear_vs_traction_abs_error"], 0.0)
    _assert_close(checks["max_torque_identity_abs_error"], 0.0)
    _assert_close(checks["max_force_identity_abs_error"], 0.0)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = build_rows()
    checks = validate(rows)
    summary = {
        "kind": "air_gap_shear_torque_validation",
        "validation_class": True,
        "geometry": {
            "radius_m": RADIUS_M,
            "axial_length_m": AXIAL_LENGTH_M,
        },
        "field": {
            "B_radial_T": BR_T,
            "B_tangential_reference_T": BT_T,
        },
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[air-gap Maxwell shear torque]")
    for row in rows:
        item = row["summary"]
        print(
            f"  {row['label']}: tau={item['shear_stress_Pa']:.6g} Pa, "
            f"Ft={item['tangential_force_N']:.6g} N, "
            f"T={item['torque_Nm']:.6g} N.m"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
