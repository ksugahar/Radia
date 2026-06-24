"""Validation-class one-port reflection momentum force.

Run:

    python examples/rf_waveguide/validation_one_port_reflection_momentum_force.py

The helper converts a measured or simulated one-port reflection coefficient
``S11`` into the time-average momentum force on the termination:

    F = (1 + |S11|^2) P_inc k_inc / c

The matched load, perfect short/open, and partial-reflection cases give compact
checks for RF post-processing before a full Maxwell-stress integration is used.
"""

from __future__ import annotations

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
    C0,
    one_port_reflection_momentum_force_summary,
)


OUT_JSON = HERE / "validation_one_port_reflection_momentum_force_summary.json"
POWER_W = 2.5
DIRECTION = (0.0, 0.0, 1.0)

CASES = (
    {"label": "matched_load", "s11": 0.0 + 0.0j, "expected_factor": 1.0},
    {"label": "perfect_short", "s11": -1.0 + 0.0j, "expected_factor": 2.0},
    {"label": "partial_reflection_plus90deg", "s11": 0.0 + 0.5j, "expected_factor": 1.25},
    {"label": "partial_reflection_minus90deg", "s11": 0.0 - 0.5j, "expected_factor": 1.25},
)


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-18) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def _json_clean(value):
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def build_summary() -> dict:
    rows = []
    p_over_c = POWER_W / C0
    for case in CASES:
        summary = one_port_reflection_momentum_force_summary(
            POWER_W,
            case["s11"],
            incident_direction=DIRECTION,
        )
        expected_force = case["expected_factor"] * p_over_c
        _assert_close(summary["reflectance"], abs(case["s11"]) ** 2)
        _assert_close(summary["axial_force_along_incident_direction_N"], expected_force)
        _assert_close(summary["force_magnitude_N"], expected_force)
        rows.append({
            "label": case["label"],
            "s11_real": case["s11"].real,
            "s11_imag": case["s11"].imag,
            "expected_factor": case["expected_factor"],
            "expected_force_N": expected_force,
            "summary": summary,
            "force_abs_error_N": abs(summary["force_magnitude_N"] - expected_force),
        })

    by_label = {row["label"]: row for row in rows}
    phase_independent_error = abs(
        by_label["partial_reflection_plus90deg"]["summary"]["force_magnitude_N"]
        - by_label["partial_reflection_minus90deg"]["summary"]["force_magnitude_N"]
    )
    return {
        "kind": "one_port_reflection_momentum_force_validation",
        "validation_class": True,
        "force_learning": "one-port S11 force: F=(1+|S11|^2)P_inc k_inc/c",
        "incident_power_W": POWER_W,
        "incident_direction": list(DIRECTION),
        "p_over_c_N": p_over_c,
        "rows": rows,
        "checks": {
            "matched_load_force_N": by_label["matched_load"]["summary"]["force_magnitude_N"],
            "perfect_short_force_N": by_label["perfect_short"]["summary"]["force_magnitude_N"],
            "partial_reflection_force_N": by_label["partial_reflection_plus90deg"]["summary"]["force_magnitude_N"],
            "phase_independent_force_error_N": phase_independent_error,
            "max_force_abs_error_N": max(row["force_abs_error_N"] for row in rows),
        },
    }


def main() -> int:
    summary = build_summary()
    OUT_JSON.write_text(json.dumps(_json_clean(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[One-port reflection momentum force]")
    for row in summary["rows"]:
        item = row["summary"]
        phase = math.degrees(math.atan2(row["s11_imag"], row["s11_real"]))
        print(
            f"  {row['label']}: |S11|={item['s11_magnitude']:.3g}, "
            f"phase={phase:.1f} deg, "
            f"factor={row['expected_factor']:.6g}, "
            f"F={item['force_magnitude_N']:.6g} N"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
