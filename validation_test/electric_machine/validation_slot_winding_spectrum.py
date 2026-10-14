"""Validation-class slot-table winding spectrum example.

This script is intentionally an example/validation run rather than a pytest
test.  It records harmonic winding factors from explicit slot sign tables,
including fractional-slot layouts that the closed-form integral-slot helper
correctly refuses.

Run:

    python validation_test/electric_machine/validation_slot_winding_spectrum.py
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
    integral_slot_winding_factor,
    mmf_harmonic_direction,
    slot_table_winding_factor,
)


OUT_JSON = HERE / "validation_slot_winding_spectrum_summary.json"
HARMONICS = (1, 3, 5, 7, 9, 11, 13)


def _angular_distance_deg(a: float, b: float) -> float:
    return min((a - b) % 360.0, (b - a) % 360.0)


def phase_belt_signs(slots: int, poles: int, phase_axis_deg: float = 0.0,
                     belt_width_deg: float = 60.0) -> list[float]:
    """One-phase sign table from a simple 60-degree phase-belt rule."""
    pole_pairs = poles // 2
    half = belt_width_deg / 2.0
    signs = []
    for k in range(slots):
        angle = (360.0 * pole_pairs * k / slots) % 360.0
        plus_axis = phase_axis_deg % 360.0
        minus_axis = (phase_axis_deg + 180.0) % 360.0
        if _angular_distance_deg(angle, plus_axis) <= half:
            signs.append(+1.0)
        elif _angular_distance_deg(angle, minus_axis) <= half:
            signs.append(-1.0)
        else:
            signs.append(0.0)
    return signs


def run_layout(label: str, slots: int, poles: int, harmonics=HARMONICS) -> dict:
    signs = phase_belt_signs(slots, poles)
    rows = []
    integral_errors = []
    for h in harmonics:
        row = slot_table_winding_factor(signs, poles, harmonic=h)
        item = {
            "harmonic": h,
            "kw_abs": row["winding_factor_abs"],
            "kw_angle_deg": row["winding_factor_angle_deg"],
            "direction": mmf_harmonic_direction(h, 3),
        }
        try:
            closed = integral_slot_winding_factor(slots, poles, harmonic=h)
        except ValueError:
            item["integral_slot_closed_form"] = None
        else:
            item["integral_slot_closed_form"] = closed["winding_factor_abs"]
            err = abs(item["kw_abs"] - closed["winding_factor_abs"])
            item["integral_slot_abs_error"] = err
            integral_errors.append(err)
        rows.append(item)
    return {
        "label": label,
        "slots": slots,
        "poles": poles,
        "slots_per_pole_per_phase": slots / (poles * 3),
        "phase_a_signs": signs,
        "active_phase_a_slots": sum(abs(s) for s in signs),
        "harmonics": rows,
        "max_integral_slot_abs_error": max(integral_errors) if integral_errors else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    layouts = [
        ("integral_36s4p", 36, 4),
        ("fractional_12s10p", 12, 10),
        ("fractional_24s22p", 24, 22),
    ]
    records = [run_layout(label, slots, poles) for label, slots, poles in layouts]
    summary = {
        "kind": "slot_table_winding_spectrum_validation",
        "validation_class": True,
        "harmonics": list(HARMONICS),
        "layouts": records,
    }
    out = args.out
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    for rec in records:
        print(
            f"[{rec['label']}] q={rec['slots_per_pole_per_phase']:.6g} "
            f"active_A={rec['active_phase_a_slots']:.0f}"
        )
        for row in rec["harmonics"]:
            print(
                f"  n={row['harmonic']:2d}  |kw|={row['kw_abs']:.6f}  "
                f"dir={row['direction']:2d}"
            )
        if rec["max_integral_slot_abs_error"] is not None:
            print(f"  integral closed-form max abs error = {rec['max_integral_slot_abs_error']:.3e}")
    print(f"[OK] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
