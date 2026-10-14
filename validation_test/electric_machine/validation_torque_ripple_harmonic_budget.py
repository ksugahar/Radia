"""Validation-class three-phase torque-ripple harmonic budget.

This example turns back-EMF harmonic content into a readable torque-ripple
budget and then shows how winding pitch and skew move that budget:

* 5th + 7th EMF harmonics create 6th torque ripple;
* 11th + 13th create 12th torque ripple;
* triplen harmonics cancel in balanced three-phase power;
* short-pitching and skew reduce the harmonic budget with a small fundamental cost.

Run:

    python validation_test/electric_machine/validation_torque_ripple_harmonic_budget.py
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

from radia_mcp.radia_ngsolve.solve import (  # noqa: E402
    skewed_winding_factor,
    slot_pitch_skew_angle,
    three_phase_torque_ripple_harmonics,
    three_phase_torque_ripple_pair_table,
)


OUT_JSON = HERE / "validation_torque_ripple_harmonic_budget_summary.json"

SLOTS = 36
POLES = 4
POLE_PAIRS = POLES // 2
PHASES = 3
Q = SLOTS // (POLES * PHASES)
SLOT_ANGLE_ELEC_DEG = 180.0 * POLES / SLOTS
CURRENT_PEAK_A = 120.0
MECHANICAL_SPEED_RAD_PER_S = 2.0 * math.pi * 1500.0 / 60.0

# Normalized open-circuit air-gap EMF content before the winding factor is applied.
# The triplen term is deliberately included; it must not appear in the torque table.
RAW_EMF_HARMONICS = {
    1: 1.0,
    3: 0.20,
    5: 0.16,
    7: 0.11,
    11: 0.045,
    13: 0.035,
}

DESIGN_CASES = (
    {"label": "full_pitch_no_skew", "pitch_fraction": 1.0, "skew_slot_pitches": 0.0},
    {"label": "five_sixth_pitch_no_skew", "pitch_fraction": 5.0 / 6.0, "skew_slot_pitches": 0.0},
    {"label": "five_sixth_pitch_one_slot_skew", "pitch_fraction": 5.0 / 6.0, "skew_slot_pitches": 1.0},
    {"label": "full_pitch_two_slot_skew", "pitch_fraction": 1.0, "skew_slot_pitches": 2.0},
)


def winding_scaled_harmonics(pitch_fraction: float, skew_slot_pitches: float) -> tuple[dict, dict]:
    skew_angle = slot_pitch_skew_angle(SLOTS, POLE_PAIRS) * skew_slot_pitches
    kw1 = skewed_winding_factor(1, Q, SLOT_ANGLE_ELEC_DEG, skew_angle, pitch_fraction)
    harmonics = {}
    factors = {}
    for order, raw in RAW_EMF_HARMONICS.items():
        kw = skewed_winding_factor(order, Q, SLOT_ANGLE_ELEC_DEG, skew_angle, pitch_fraction)
        factors[order] = {
            "winding_skew_factor": kw,
            "relative_to_fundamental": kw / kw1 if kw1 != 0.0 else math.nan,
        }
        harmonics[order] = raw * factors[order]["relative_to_fundamental"]
    harmonics[1] = 1.0
    return harmonics, {
        "fundamental_winding_skew_factor": kw1,
        "skew_angle_elec_rad": skew_angle,
        "harmonic_factors": factors,
    }


def _sample_power(emf_harmonics: dict[int, float], samples: int = 4096) -> list[float]:
    shifts = (0.0, -2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0)
    vals = []
    for idx in range(samples):
        theta = 2.0 * math.pi * idx / samples
        power = 0.0
        for shift in shifts:
            angle = theta + shift
            emf = sum(
                complex(value * complex(math.cos(order * angle), math.sin(order * angle))).real
                for order, value in emf_harmonics.items()
            )
            current = CURRENT_PEAK_A * math.cos(angle)
            power += emf * current
        vals.append(power)
    return vals


def _fourier_amplitude(values: list[float], order: int) -> float:
    n = len(values)
    coeff = sum(
        values[idx] * complex(
            math.cos(-2.0 * math.pi * order * idx / n),
            math.sin(-2.0 * math.pi * order * idx / n),
        )
        for idx in range(n)
    ) / n
    return 2.0 * abs(coeff)


def build_case(case: dict) -> dict:
    harmonics, factors = winding_scaled_harmonics(
        case["pitch_fraction"],
        case["skew_slot_pitches"],
    )
    summary = three_phase_torque_ripple_harmonics(
        harmonics,
        current_peak=CURRENT_PEAK_A,
        mechanical_speed=MECHANICAL_SPEED_RAD_PER_S,
    )
    table = three_phase_torque_ripple_pair_table(
        harmonics,
        current_peak=CURRENT_PEAK_A,
        mechanical_speed=MECHANICAL_SPEED_RAD_PER_S,
    )
    waveform = _sample_power(harmonics)
    fourier_errors = {}
    for row in table:
        order = row["ripple_order"]
        fourier_errors[order] = abs(_fourier_amplitude(waveform, order) - row["power_ripple"])

    return {
        **case,
        **factors,
        "emf_harmonics": harmonics,
        "mean_power": summary["mean_power"],
        "mean_torque": summary["mean_torque"],
        "pair_table": table,
        "normalized_ripple": summary["normalized_ripple"],
        "fourier_power_ripple_abs_errors": fourier_errors,
    }


def validate(rows: list[dict]) -> dict:
    by_label = {row["label"]: row for row in rows}
    full = by_label["full_pitch_no_skew"]
    short = by_label["five_sixth_pitch_no_skew"]
    short_skew = by_label["five_sixth_pitch_one_slot_skew"]
    two_slot = by_label["full_pitch_two_slot_skew"]

    def ripple(row: dict, order: int) -> float:
        return row["normalized_ripple"].get(order, 0.0)

    all_orders = {
        row["label"]: [entry["ripple_order"] for entry in row["pair_table"]]
        for row in rows
    }
    all_fourier_errors = [
        err
        for row in rows
        for err in row["fourier_power_ripple_abs_errors"].values()
    ]

    checks = {
        "all_tables_have_only_6_and_12": all(orders == [6, 12] for orders in all_orders.values()),
        "triplen_absent_from_pair_tables": all(3 not in orders for orders in all_orders.values()),
        "short_pitch_reduces_6th": ripple(short, 6) < ripple(full, 6),
        "short_pitch_preserves_12th_for_this_layout": math.isclose(
            ripple(short, 12), ripple(full, 12), rel_tol=1.0e-12, abs_tol=1.0e-15
        ),
        "one_slot_skew_reduces_12th": ripple(short_skew, 12) < ripple(short, 12),
        "two_slot_skew_reduces_6th_vs_full": ripple(two_slot, 6) < ripple(full, 6),
        "two_slot_skew_reduces_12th_vs_full": ripple(two_slot, 12) < ripple(full, 12),
        "max_fourier_power_ripple_abs_error": max(all_fourier_errors),
        "full_pitch_6th_normalized": ripple(full, 6),
        "short_pitch_6th_normalized": ripple(short, 6),
        "short_pitch_one_slot_skew_6th_normalized": ripple(short_skew, 6),
        "full_pitch_two_slot_skew_6th_normalized": ripple(two_slot, 6),
        "fundamental_factor_cost_one_slot_short_pitch": (
            short_skew["fundamental_winding_skew_factor"]
            / full["fundamental_winding_skew_factor"]
        ),
    }
    assert checks["all_tables_have_only_6_and_12"]
    assert checks["triplen_absent_from_pair_tables"]
    assert checks["short_pitch_reduces_6th"]
    assert checks["short_pitch_preserves_12th_for_this_layout"]
    assert checks["one_slot_skew_reduces_12th"]
    assert checks["two_slot_skew_reduces_6th_vs_full"]
    assert checks["two_slot_skew_reduces_12th_vs_full"]
    assert checks["max_fourier_power_ripple_abs_error"] < 1.0e-10
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = [build_case(case) for case in DESIGN_CASES]
    checks = validate(rows)
    summary = {
        "kind": "torque_ripple_harmonic_budget_validation",
        "validation_class": True,
        "machine": {
            "slots": SLOTS,
            "poles": POLES,
            "phases": PHASES,
            "slots_per_pole_per_phase": Q,
            "slot_angle_elec_deg": SLOT_ANGLE_ELEC_DEG,
            "current_peak_A": CURRENT_PEAK_A,
            "mechanical_speed_rad_per_s": MECHANICAL_SPEED_RAD_PER_S,
        },
        "raw_emf_harmonics": RAW_EMF_HARMONICS,
        "cases": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[three-phase torque-ripple harmonic budget]")
    for row in rows:
        print(
            f"  {row['label']}: kw1={row['fundamental_winding_skew_factor']:.6f}  "
            f"r6={row['normalized_ripple'].get(6, 0.0):.6f}  "
            f"r12={row['normalized_ripple'].get(12, 0.0):.6f}"
        )
        for entry in row["pair_table"]:
            print(
                f"    order {entry['ripple_order']:2d}: "
                f"harmonics={entry['contributing_harmonics']}  "
                f"|Epair|={entry['emf_phasor_abs']:.6f}  "
                f"T={entry['torque_ripple']:.6f} N.m"
            )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
