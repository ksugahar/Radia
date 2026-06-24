"""Validation-class PM load-line and irreversible-demag margin sweep.

This is an example/validation run rather than a pytest test. It turns the
compact PM circuit formulas into readable tables:

* larger air gaps lower the load-line permeance coefficient and gap flux;
* longer magnets move the operating point away from the knee;
* hotter/weaker knee assumptions can flip the same circuit from safe to unsafe;
* simple shape demagnetizing factors explain why thin magnets need care.

Run:

    python examples/electric_machine/validation_pm_loadline_demag_sweep.py
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
    MU0,
    demag_margin,
    demag_operating_field,
    pm_circuit_loadline_operating_point,
)


OUT_JSON = HERE / "validation_pm_loadline_demag_sweep_summary.json"

BASE = {
    "Br_T": 1.20,
    "magnet_len_m": 6.0e-3,
    "gap_m": 6.0e-3,
    "iron_path_m": 0.12,
    "mu_r_iron": 2000.0,
    "mu_rec": 1.05,
    "H_knee_A_per_m": -6.0e5,
}

GAPS_M = (0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3)
MAGNET_LENGTHS_M = (3.0e-3, 6.0e-3, 12.0e-3)
THERMAL_KNEE_CASES = (
    {"label": "cool_strong_knee", "Br_T": 1.25, "H_knee_A_per_m": -9.0e5},
    {"label": "room_reference", "Br_T": 1.20, "H_knee_A_per_m": -6.0e5},
    {"label": "hot_weak_knee", "Br_T": 1.08, "H_knee_A_per_m": -3.90e5},
)
SHAPE_CASES = (
    {"label": "long_axis_magnet", "demag_factor": 0.10},
    {"label": "sphere", "demag_factor": 1.0 / 3.0},
    {"label": "transverse_cylinder", "demag_factor": 0.50},
    {"label": "thin_slab_like", "demag_factor": 0.80},
)


def _loadline_row(label: str, **params) -> dict:
    op = pm_circuit_loadline_operating_point(
        params["Br_T"],
        params["magnet_len_m"],
        params["gap_m"],
        params["iron_path_m"],
        params["mu_r_iron"],
        params["mu_rec"],
        H_knee=params["H_knee_A_per_m"],
    )
    return {
        "label": label,
        **params,
        **op,
    }


def build_gap_sweep() -> list[dict]:
    rows = []
    for gap in GAPS_M:
        params = dict(BASE)
        params["gap_m"] = gap
        rows.append(_loadline_row(f"gap_{1e3 * gap:.1f}mm", **params))
    return rows


def build_magnet_length_sweep() -> list[dict]:
    rows = []
    for length in MAGNET_LENGTHS_M:
        params = dict(BASE)
        params["magnet_len_m"] = length
        rows.append(_loadline_row(f"lm_{1e3 * length:.1f}mm", **params))
    return rows


def build_thermal_knee_sweep() -> list[dict]:
    rows = []
    for case in THERMAL_KNEE_CASES:
        params = dict(BASE)
        params["Br_T"] = case["Br_T"]
        params["H_knee_A_per_m"] = case["H_knee_A_per_m"]
        rows.append(_loadline_row(case["label"], **params))
    return rows


def build_shape_sweep() -> list[dict]:
    rows = []
    for case in SHAPE_CASES:
        h_op = demag_operating_field(BASE["Br_T"], case["demag_factor"], BASE["mu_rec"])
        margin = demag_margin(h_op, BASE["H_knee_A_per_m"])
        rows.append({
            "label": case["label"],
            "Br_T": BASE["Br_T"],
            "mu_rec": BASE["mu_rec"],
            "demag_factor": case["demag_factor"],
            "H_operating_A_per_m": h_op,
            "H_knee_A_per_m": BASE["H_knee_A_per_m"],
            "demag_margin_A_per_m": margin,
            "safe_against_knee": margin >= 0.0,
            "H_over_Hc": h_op / (BASE["Br_T"] / (MU0 * BASE["mu_rec"])),
        })
    return rows


def _is_strict_decreasing(values: list[float]) -> bool:
    return all(a > b for a, b in zip(values, values[1:]))


def _is_strict_increasing(values: list[float]) -> bool:
    return all(a < b for a, b in zip(values, values[1:]))


def validate(gap_rows: list[dict], length_rows: list[dict],
             thermal_rows: list[dict], shape_rows: list[dict]) -> dict:
    gap_b = [row["B_gap_T"] for row in gap_rows]
    gap_h = [row["H_m_A_per_m"] for row in gap_rows]
    gap_pc = [row["permeance_coefficient"] for row in gap_rows]
    length_b = [row["B_gap_T"] for row in length_rows]
    length_h = [row["H_m_A_per_m"] for row in length_rows]
    shape_h = [row["H_operating_A_per_m"] for row in shape_rows]
    identity_errors = [
        row["B_identity_abs_error_T"]
        for row in [*gap_rows, *length_rows, *thermal_rows]
    ]

    checks = {
        "gap_B_monotone_decreasing": _is_strict_decreasing(gap_b),
        "gap_Hm_monotone_decreasing": _is_strict_decreasing(gap_h),
        "gap_permeance_monotone_decreasing": _is_strict_decreasing(gap_pc),
        "length_B_monotone_increasing": _is_strict_increasing(length_b),
        "length_Hm_monotone_increasing": _is_strict_increasing(length_h),
        "length_safe_flags": [row["safe_against_knee"] for row in length_rows],
        "thermal_safe_flags": [row["safe_against_knee"] for row in thermal_rows],
        "shape_H_monotone_decreasing": _is_strict_decreasing(shape_h),
        "shape_safe_flags": [row["safe_against_knee"] for row in shape_rows],
        "max_B_identity_abs_error_T": max(identity_errors),
    }
    assert checks["gap_B_monotone_decreasing"]
    assert checks["gap_Hm_monotone_decreasing"]
    assert checks["gap_permeance_monotone_decreasing"]
    assert checks["length_B_monotone_increasing"]
    assert checks["length_Hm_monotone_increasing"]
    assert checks["length_safe_flags"] == [False, True, True]
    assert checks["thermal_safe_flags"] == [True, True, False]
    assert checks["shape_H_monotone_decreasing"]
    assert checks["shape_safe_flags"] == [True, True, True, False]
    assert checks["max_B_identity_abs_error_T"] < 1.0e-12
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    gap_rows = build_gap_sweep()
    length_rows = build_magnet_length_sweep()
    thermal_rows = build_thermal_knee_sweep()
    shape_rows = build_shape_sweep()
    checks = validate(gap_rows, length_rows, thermal_rows, shape_rows)

    summary = {
        "kind": "pm_loadline_demag_margin_validation",
        "validation_class": True,
        "base": BASE,
        "gap_sweep": gap_rows,
        "magnet_length_sweep": length_rows,
        "thermal_knee_sweep": thermal_rows,
        "shape_demag_sweep": shape_rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[PM load-line / demag margin gap sweep]")
    for row in gap_rows:
        print(
            f"  g={1e3 * row['gap_m']:4.1f} mm  "
            f"Pc={row['permeance_coefficient']:7.3f}  "
            f"B={row['B_gap_T']:.6f} T  "
            f"Hm={row['H_m_A_per_m'] / 1e3:8.3f} kA/m  "
            f"margin={row['demag_margin_A_per_m'] / 1e3:8.3f} kA/m"
        )

    print("[magnet-length knee check]")
    for row in length_rows:
        print(
            f"  lm={1e3 * row['magnet_len_m']:4.1f} mm  "
            f"B={row['B_gap_T']:.6f} T  "
            f"Hm={row['H_m_A_per_m'] / 1e3:8.3f} kA/m  "
            f"safe={row['safe_against_knee']}"
        )

    print("[shape demag check]")
    for row in shape_rows:
        print(
            f"  {row['label']:<20} N={row['demag_factor']:.6f}  "
            f"H={row['H_operating_A_per_m'] / 1e3:8.3f} kA/m  "
            f"safe={row['safe_against_knee']}"
        )

    print(
        "[checks] "
        f"length safe flags={checks['length_safe_flags']}, "
        f"thermal safe flags={checks['thermal_safe_flags']}, "
        f"shape safe flags={checks['shape_safe_flags']}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
