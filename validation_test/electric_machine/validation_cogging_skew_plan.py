"""Validation-class cogging order and skew planning table.

This is an example/validation run rather than a pytest test. It turns the
slot/pole topology rules into a readable design table:

* cogging order is LCM(slots, poles);
* the minimum symmetry sector is GCD(slots, poles);
* one stator-slot-pitch skew cancels the cogging order;
* the same skew has a visible fundamental-EMF cost, especially in low-slot
  fractional-slot layouts.

Run:

    python validation_test/electric_machine/validation_cogging_skew_plan.py
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

from radia_mcp.radia_ngsolve.solve import cogging_skew_plan  # noqa: E402


OUT_JSON = HERE / "validation_cogging_skew_plan_summary.json"

CASES = [
    {"label": "integer_36s_4p", "slots": 36, "poles": 4},
    {"label": "fractional_12s_10p", "slots": 12, "poles": 10},
    {"label": "fractional_9s_8p", "slots": 9, "poles": 8},
    {"label": "distributed_24s_8p", "slots": 24, "poles": 8},
    {"label": "pathological_12s_12p", "slots": 12, "poles": 12},
]


def _row(case: dict, skew_slot_pitches: float) -> dict:
    plan = cogging_skew_plan(
        case["slots"],
        case["poles"],
        skew_slot_pitches=skew_slot_pitches,
        emf_harmonics=(1, 5, 7, 11, 13),
    )
    return {
        "label": case["label"],
        **plan,
    }


def build_rows() -> tuple[list[dict], list[dict]]:
    one_slot = [_row(case, 1.0) for case in CASES]
    half_slot = [_row(case, 0.5) for case in CASES]
    return one_slot, half_slot


def validate(one_slot: list[dict], half_slot: list[dict]) -> dict:
    by_label = {row["label"]: row for row in one_slot}
    half_by_label = {row["label"]: row for row in half_slot}
    max_one_slot_cogging = max(abs(row["cogging_skew_factor"]) for row in one_slot)
    min_fundamental_factor = min(row["fundamental_emf_skew_factor"] for row in one_slot)
    equal_slot_period = by_label["pathological_12s_12p"]["cogging_period_mech_deg"]
    equal_slot_fundamental = by_label["pathological_12s_12p"]["fundamental_emf_skew_factor"]
    fractional_period = by_label["fractional_12s_10p"]["cogging_period_mech_deg"]
    integer_period = by_label["integer_36s_4p"]["cogging_period_mech_deg"]
    half_vs_one = {
        label: abs(half_by_label[label]["cogging_skew_factor"]) >= abs(one["cogging_skew_factor"])
        for label, one in by_label.items()
    }

    checks = {
        "max_one_slot_cogging_skew_abs": max_one_slot_cogging,
        "min_one_slot_fundamental_emf_skew_factor": min_fundamental_factor,
        "fractional_12s_10p_period_below_integer_36s_4p": fractional_period < integer_period,
        "pathological_12s_12p_period_deg": equal_slot_period,
        "pathological_12s_12p_fundamental_emf_skew_factor": equal_slot_fundamental,
        "half_slot_cogging_factor_not_better_than_one_slot": half_vs_one,
        "one_slot_all_cogging_cancelled": all(abs(row["cogging_skew_factor"]) < 1.0e-12 for row in one_slot),
        "one_slot_fundamental_cost_bounded": min_fundamental_factor > 0.60,
    }
    assert checks["one_slot_all_cogging_cancelled"]
    assert checks["one_slot_fundamental_cost_bounded"]
    assert checks["fractional_12s_10p_period_below_integer_36s_4p"]
    assert math.isclose(equal_slot_period, 30.0, rel_tol=1.0e-12)
    assert equal_slot_fundamental < 0.65
    assert all(half_vs_one.values())
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    one_slot, half_slot = build_rows()
    checks = validate(one_slot, half_slot)
    summary = {
        "kind": "cogging_skew_plan_validation",
        "validation_class": True,
        "one_slot_pitch_skew": one_slot,
        "half_slot_pitch_skew": half_slot,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[Cogging order / skew planning]")
    for row in one_slot:
        sym = row["symmetry"]
        print(
            f"  {row['label']:<24} "
            f"order={row['cogging_order_per_rev']:3d}/rev  "
            f"period={row['cogging_period_mech_deg']:6.2f} deg  "
            f"sector={sym['sector_angle_deg']:6.1f} deg/{sym['boundary']:<13}  "
            f"skew={row['skew_angle_mech_deg']:6.2f} mech deg  "
            f"k_cog={row['cogging_skew_factor']:+.3e}  "
            f"k_emf1={row['fundamental_emf_skew_factor']:.6f}"
        )
    print(
        "[checks] "
        f"max |k_cog(one-slot)|={checks['max_one_slot_cogging_skew_abs']:.3e}, "
        f"min k_emf1={checks['min_one_slot_fundamental_emf_skew_factor']:.6f}, "
        f"12s/12p period={checks['pathological_12s_12p_period_deg']:.1f} deg"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
