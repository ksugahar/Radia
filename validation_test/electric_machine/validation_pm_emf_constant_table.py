"""Validation-class PM back-EMF / torque-constant table.

This PM-machine teaching example turns PM flux linkage into the constants that
normally appear in motor reports:

* phase-peak, phase-RMS, and line-line-RMS back-EMF constants;
* q-axis peak/RMS torque constants;
* no-load back-EMF at selected speeds;
* the dq power identity, T*omega_mech = (3/2) e_phase_peak iq_peak.

Run:

    python validation_test/electric_machine/validation_pm_emf_constant_table.py
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
    dq_voltages,
    pm_flux_linkage_constants,
    pm_no_load_back_emf,
)


OUT_JSON = HERE / "validation_pm_emf_constant_table_summary.json"

CASES = [
    {"label": "baseline_4pp", "lambda_m": 0.10, "pole_pairs": 4, "iq_peak": 12.0},
    {"label": "compact_3pp", "lambda_m": 0.045, "pole_pairs": 3, "iq_peak": 18.0},
    {"label": "high_flux_5pp", "lambda_m": 0.075, "pole_pairs": 5, "iq_peak": 10.0},
]
OMEGA_MECH = (50.0, 150.0, 300.0)


def run_case(case: dict) -> dict:
    constants = pm_flux_linkage_constants(case["lambda_m"], case["pole_pairs"])
    speed_rows = []
    for omega_mech in OMEGA_MECH:
        emf = pm_no_load_back_emf(case["lambda_m"], omega_mech, case["pole_pairs"])
        vd, vq = dq_voltages(
            R=0.0,
            Ld=1.0e-3,
            Lq=1.0e-3,
            lambda_m=case["lambda_m"],
            id_=0.0,
            iq=0.0,
            omega_e=emf["omega_e_rad_per_s"],
        )
        torque = dq_torque(
            case["lambda_m"],
            Ld=1.0e-3,
            Lq=1.0e-3,
            id_=0.0,
            iq=case["iq_peak"],
            pole_pairs=case["pole_pairs"],
        )
        mechanical_power = torque * omega_mech
        electrical_power = 1.5 * emf["phase_peak_V"] * case["iq_peak"]
        speed_rows.append({
            "omega_mech_rad_per_s": omega_mech,
            "omega_e_rad_per_s": emf["omega_e_rad_per_s"],
            "phase_peak_V": emf["phase_peak_V"],
            "phase_rms_V": emf["phase_rms_V"],
            "line_line_rms_V": emf["line_line_rms_V"],
            "dq_vd_no_load": vd,
            "dq_vq_no_load": vq,
            "torque_Nm": torque,
            "mechanical_power_W": mechanical_power,
            "electrical_power_W": electrical_power,
            "dq_emf_abs_error_V": abs(vq - emf["phase_peak_V"]) + abs(vd),
            "power_abs_error_W": abs(mechanical_power - electrical_power),
        })
    return {
        "label": case["label"],
        "case": case,
        "constants": constants,
        "speed_rows": speed_rows,
    }


def validate(rows: list[dict]) -> dict:
    all_speed_rows = [speed for row in rows for speed in row["speed_rows"]]
    checks = {
        "max_dq_emf_abs_error_V": max(speed["dq_emf_abs_error_V"] for speed in all_speed_rows),
        "max_power_abs_error_W": max(speed["power_abs_error_W"] for speed in all_speed_rows),
        "max_Kt_rms_over_Ke_ll_rms_error": max(
            abs(row["constants"]["Kt_rms_over_line_line_rms_Ke"] - math.sqrt(3.0))
            for row in rows
        ),
        "max_Kt_peak_over_phase_peak_Ke_error": max(
            abs(row["constants"]["Kt_peak_over_phase_peak_Ke"] - 1.5)
            for row in rows
        ),
        "baseline_line_line_rms_Ke": rows[0]["constants"][
            "back_emf_constant_line_line_rms_V_per_rad_per_s_mech"
        ],
        "baseline_Kt_peak": rows[0]["constants"]["torque_constant_Nm_per_Aq_peak"],
    }
    assert checks["max_dq_emf_abs_error_V"] == 0.0
    assert checks["max_power_abs_error_W"] < 1.0e-12
    assert checks["max_Kt_rms_over_Ke_ll_rms_error"] < 1.0e-15
    assert checks["max_Kt_peak_over_phase_peak_Ke_error"] < 1.0e-15
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = [run_case(case) for case in CASES]
    checks = validate(rows)
    summary = {
        "kind": "pm_emf_constant_table_validation",
        "validation_class": True,
        "omega_mech_rad_per_s": list(OMEGA_MECH),
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[PM back-EMF / torque constants]")
    for row in rows:
        c = row["constants"]
        print(
            f"  {row['label']}: "
            f"Ke_ll_rms={c['back_emf_constant_line_line_rms_V_per_rad_per_s_mech']:.9f} "
            f"Kt_peak={c['torque_constant_Nm_per_Aq_peak']:.9f} "
            f"Kt_rms/Ke_ll={c['Kt_rms_over_line_line_rms_Ke']:.9f}"
        )
        for speed in row["speed_rows"]:
            print(
                f"    omega={speed['omega_mech_rad_per_s']:6.1f} rad/s "
                f"E_ll_rms={speed['line_line_rms_V']:.6f} V "
                f"P_err={speed['power_abs_error_W']:.3e} W"
            )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
