"""Validation-class MTPA saliency sweep for dq machine models.

This example is intentionally outside the unit-test path. It compares the
closed-form maximum-torque-per-ampere operating point with an independent
current-angle sweep across non-salient PM, IPM, SynRM, and salient-PM cases.

Run:

    python validation_test/electric_machine/validation_mtpa_saliency_sweep.py
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
    dq_torque,
    dq_torque_components,
    mtpa_operating_point,
)


OUT_JSON = HERE / "validation_mtpa_saliency_sweep_summary.json"

CASES = [
    {
        "label": "surface_pm_nonsalient",
        "lambda_m": 0.12,
        "Ld": 3.0e-3,
        "Lq": 3.0e-3,
        "current": 20.0,
        "pole_pairs": 4,
    },
    {
        "label": "ipm_strong_saliency",
        "lambda_m": 0.12,
        "Ld": 2.0e-3,
        "Lq": 6.0e-3,
        "current": 20.0,
        "pole_pairs": 4,
    },
    {
        "label": "synrel_no_magnet",
        "lambda_m": 0.0,
        "Ld": 6.0e-3,
        "Lq": 2.0e-3,
        "current": 20.0,
        "pole_pairs": 2,
    },
    {
        "label": "salient_pm_Ld_gt_Lq",
        "lambda_m": 0.08,
        "Ld": 3.6e-3,
        "Lq": 2.4e-3,
        "current": 18.0,
        "pole_pairs": 3,
    },
]


def numeric_argmax(case: dict, n: int = 20001) -> dict:
    lm = case["lambda_m"]
    ld = case["Ld"]
    lq = case["Lq"]
    current = case["current"]
    p = case["pole_pairs"]
    best_gamma = None
    best_torque = None
    for i in range(n):
        gamma = -0.5 * math.pi + math.pi * i / (n - 1)
        id_ = -current * math.sin(gamma)
        iq = current * math.cos(gamma)
        torque = dq_torque(lm, ld, lq, id_, iq, p)
        if best_torque is None or torque > best_torque:
            best_gamma = gamma
            best_torque = torque
    return {
        "gamma_rad": best_gamma,
        "gamma_deg": math.degrees(best_gamma),
        "torque_Nm": best_torque,
    }


def run_case(case: dict) -> dict:
    gamma, id_, iq, torque = mtpa_operating_point(
        case["lambda_m"],
        case["Ld"],
        case["Lq"],
        case["current"],
        case["pole_pairs"],
    )
    numeric = numeric_argmax(case)
    t_mag, t_rel, t_total = dq_torque_components(
        case["lambda_m"],
        case["Ld"],
        case["Lq"],
        id_,
        iq,
        case["pole_pairs"],
    )
    pure_q = dq_torque(
        case["lambda_m"],
        case["Ld"],
        case["Lq"],
        0.0,
        case["current"],
        case["pole_pairs"],
    )
    torque_gain = torque / pure_q if abs(pure_q) > 0.0 else None
    return {
        "label": case["label"],
        "case": case,
        "closed_form": {
            "gamma_rad": gamma,
            "gamma_deg": math.degrees(gamma),
            "id_A": id_,
            "iq_A": iq,
            "torque_Nm": torque,
            "magnet_torque_Nm": t_mag,
            "reluctance_torque_Nm": t_rel,
            "components_total_Nm": t_total,
            "pure_q_torque_Nm": pure_q,
            "torque_gain_vs_pure_q": torque_gain,
        },
        "numeric_argmax": numeric,
        "errors": {
            "gamma_abs_error_deg": abs(math.degrees(gamma) - numeric["gamma_deg"]),
            "torque_abs_error_Nm": abs(torque - numeric["torque_Nm"]),
            "torque_rel_error": (
                abs(torque - numeric["torque_Nm"]) / max(1.0e-30, abs(numeric["torque_Nm"]))
            ),
            "component_sum_abs_error_Nm": abs(t_total - torque),
            "current_magnitude_error_A": abs(math.hypot(id_, iq) - case["current"]),
        },
    }


def validate(rows: list[dict]) -> dict:
    by_label = {row["label"]: row for row in rows}
    max_gamma_error = max(row["errors"]["gamma_abs_error_deg"] for row in rows)
    max_torque_rel_error = max(row["errors"]["torque_rel_error"] for row in rows)
    max_component_error = max(row["errors"]["component_sum_abs_error_Nm"] for row in rows)
    checks = {
        "max_gamma_abs_error_deg": max_gamma_error,
        "max_torque_rel_error": max_torque_rel_error,
        "max_component_sum_abs_error_Nm": max_component_error,
        "surface_pm_gamma_deg": by_label["surface_pm_nonsalient"]["closed_form"]["gamma_deg"],
        "ipm_gamma_deg": by_label["ipm_strong_saliency"]["closed_form"]["gamma_deg"],
        "ipm_torque_gain_vs_pure_q": (
            by_label["ipm_strong_saliency"]["closed_form"]["torque_gain_vs_pure_q"]
        ),
        "synrel_gamma_deg": by_label["synrel_no_magnet"]["closed_form"]["gamma_deg"],
        "salient_pm_Ld_gt_Lq_gamma_deg": (
            by_label["salient_pm_Ld_gt_Lq"]["closed_form"]["gamma_deg"]
        ),
    }
    assert max_gamma_error < 0.02
    assert max_torque_rel_error < 1.0e-7
    assert max_component_error < 1.0e-12
    assert abs(checks["surface_pm_gamma_deg"]) < 1.0e-12
    assert checks["ipm_gamma_deg"] > 0.0
    assert checks["ipm_torque_gain_vs_pure_q"] > 1.05
    assert abs(checks["synrel_gamma_deg"] + 45.0) < 1.0e-9
    assert checks["salient_pm_Ld_gt_Lq_gamma_deg"] < 0.0
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = [run_case(case) for case in CASES]
    checks = validate(rows)
    summary = {
        "kind": "mtpa_saliency_sweep_validation",
        "validation_class": True,
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[MTPA saliency sweep]")
    for row in rows:
        cf = row["closed_form"]
        err = row["errors"]
        gain = cf["torque_gain_vs_pure_q"]
        gain_text = "n/a" if gain is None else f"{gain:.6f}"
        print(
            f"  {row['label']}: gamma={cf['gamma_deg']:+8.4f} deg  "
            f"id={cf['id_A']:+8.4f} A  iq={cf['iq_A']:+8.4f} A  "
            f"T={cf['torque_Nm']:+9.5f} Nm  gain={gain_text}  "
            f"num_err={err['torque_rel_error']:.3e}"
        )
    print(
        "[checks] "
        f"max gamma error={checks['max_gamma_abs_error_deg']:.6f} deg, "
        f"max torque rel error={checks['max_torque_rel_error']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
