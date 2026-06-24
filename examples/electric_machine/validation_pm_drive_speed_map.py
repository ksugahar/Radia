"""Validation-class PM machine drive and short-circuit speed map.

This example stitches together the readable dq helpers used in JMAG/FEMM-style
motor studies:

* MTPA -> field weakening -> MTPV region selection
* finite-speed cutoff when characteristic current exceeds Imax
* short-circuit current / demagnetising fraction / braking torque trend

It is intentionally an example/validation run, not a pytest test.

Run:

    python examples/electric_machine/validation_pm_drive_speed_map.py
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

from radia_mcp.radia_ngsolve.solve import (  # noqa: E402
    pm_drive_speed_sweep,
    short_circuit_operating_point,
)


OUT_JSON = HERE / "validation_pm_drive_speed_map_summary.json"


MACHINES = [
    {
        "label": "wide_cpsr_ipm",
        "lambda_m": 0.1,
        "Ld": 8.0e-3,
        "Lq": 16.0e-3,
        "Imax": 20.0,
        "Vmax": 120.0,
        "pole_pairs": 4,
    },
    {
        "label": "finite_speed_pm",
        "lambda_m": 0.1,
        "Ld": 0.8e-3,
        "Lq": 1.6e-3,
        "Imax": 20.0,
        "Vmax": 120.0,
        "pole_pairs": 4,
    },
]


SPEED_MULTIPLES = (0.5, 1.0, 2.0, 5.0, 10.0, 100.0)
SHORT_CIRCUIT_HZ = (1.0, 20.0, 100.0, 1.0e6)


def drive_map(machine: dict) -> dict:
    sweep = pm_drive_speed_sweep(
        machine["lambda_m"],
        machine["Ld"],
        machine["Lq"],
        machine["Imax"],
        machine["Vmax"],
        machine["pole_pairs"],
        speed_multiples=SPEED_MULTIPLES,
        R=0.0,
    )
    return {
        "label": machine["label"],
        **sweep,
    }


def short_circuit_map(machine: dict, resistance: float = 0.05) -> dict:
    rows = []
    for f in SHORT_CIRCUIT_HZ:
        omega_e = 2.0 * math.pi * f
        row = short_circuit_operating_point(
            resistance,
            machine["Ld"],
            machine["Lq"],
            machine["lambda_m"],
            omega_e,
            machine["pole_pairs"],
        )
        row["frequency_hz"] = f
        rows.append(row)
    return {"resistance": resistance, "rows": rows}


def main() -> int:
    drive = [drive_map(m) for m in MACHINES]
    short = {m["label"]: short_circuit_map(m) for m in MACHINES}
    summary = {
        "kind": "pm_drive_speed_map_validation",
        "validation_class": True,
        "speed_multiples": list(SPEED_MULTIPLES),
        "drive_maps": drive,
        "short_circuit_maps": short,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    for rec in drive:
        cap = rec["speed_capability"]
        print(
            f"[{rec['label']}] Ich={rec['speed_capability']['characteristic_current']:.6g} A, "
            f"Imax={rec['parameters']['Imax']:.6g} A, "
            f"infinite_speed={cap['infinite_speed_possible']}"
        )
        for row in rec["rows"]:
            if row["region"] == "infeasible":
                print(f"  {row['speed_multiple']:6.1f} x base: infeasible")
            else:
                print(
                    f"  {row['speed_multiple']:6.1f} x base: {row['region']:4s} "
                    f"id={row['id']: .6g} iq={row['iq']: .6g} "
                    f"T={row['torque']: .6g} "
                    f"I/Imax={row['current_utilization']:.6f} "
                    f"V/Vmax={row['voltage_utilization']:.6f}"
                )
    for label, rec in short.items():
        high = rec["rows"][-1]
        print(
            f"[{label} short] {high['frequency_hz']:.0f} Hz: "
            f"|I|/Ich={high['current_ratio_to_characteristic']:.6f}, "
            f"demag={high['d_axis_demag_fraction']:.6f}, "
            f"T={high['torque']:.6g}"
        )
    print(f"[OK] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
