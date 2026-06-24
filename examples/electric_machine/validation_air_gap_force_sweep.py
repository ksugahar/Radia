"""Validation-class air-gap Maxwell pressure / holding-force sweep.

This is a lightweight FEMM-style teaching example: solve a readable series
magnetic circuit for the gap flux density, then convert that field to a pole
face pressure and holding force using Maxwell stress,

    p = B^2 / (2 mu0),   F = p A n_faces.

Run:

    python examples/electric_machine/validation_air_gap_force_sweep.py
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
    air_gap_force_summary,
    air_gap_maxwell_pressure,
)
from radia_mcp.radia_ngsolve.solve import magnetic_circuit_bh_operating_summary  # noqa: E402


OUT_JSON = HERE / "validation_air_gap_force_sweep_summary.json"

IRON_PATH_M = 0.0156
POLE_AREA_M2 = 2.5e-4
ACTIVE_FACES = 2
MMF_A_TURN = 96.0
GAP_SWEEP_M = (0.0, 0.1e-3, 0.2e-3, 0.5e-3, 1.0e-3, 2.0e-3)


def h_of_b(B: float) -> float:
    return 51.0 * B + 2.5 * B ** 15


def dh_db(B: float) -> float:
    return 51.0 + 37.5 * B ** 14


def build_rows() -> list[dict]:
    rows = []
    for gap in GAP_SWEEP_M:
        circuit = magnetic_circuit_bh_operating_summary(
            MMF_A_TURN,
            IRON_PATH_M,
            h_of_b,
            gap=gap,
            dh_db=dh_db,
        )
        force = air_gap_force_summary(circuit["B_T"], POLE_AREA_M2, faces=ACTIVE_FACES)
        pressure_ref = air_gap_maxwell_pressure(circuit["B_T"])
        force_ref = pressure_ref * POLE_AREA_M2 * ACTIVE_FACES
        rows.append({
            "gap_m": gap,
            "B_T": circuit["B_T"],
            "H_iron_A_per_m": circuit["H_iron_A_per_m"],
            "incremental_mu_r": circuit["incremental_mu_r"],
            "gap_mmf_fraction": circuit["gap_mmf_fraction"],
            "pressure_Pa": force["pressure_Pa"],
            "force_N": force["force_N"],
            "pressure_identity_abs_error": abs(force["pressure_Pa"] - pressure_ref),
            "force_identity_abs_error": abs(force["force_N"] - force_ref),
        })
    return rows


def _strictly_decreasing(values: list[float]) -> bool:
    return all(a > b for a, b in zip(values, values[1:]))


def _strictly_increasing(values: list[float]) -> bool:
    return all(a < b for a, b in zip(values, values[1:]))


def validate(rows: list[dict]) -> dict:
    B = [row["B_T"] for row in rows]
    force = [row["force_N"] for row in rows]
    pressure = [row["pressure_Pa"] for row in rows]
    gap_fraction = [row["gap_mmf_fraction"] for row in rows[1:]]
    checks = {
        "B_monotone_decreasing_with_gap": _strictly_decreasing(B),
        "pressure_monotone_decreasing_with_gap": _strictly_decreasing(pressure),
        "force_monotone_decreasing_with_gap": _strictly_decreasing(force),
        "gap_mmf_fraction_monotone_increasing": _strictly_increasing(gap_fraction),
        "max_pressure_identity_abs_error": max(row["pressure_identity_abs_error"] for row in rows),
        "max_force_identity_abs_error": max(row["force_identity_abs_error"] for row in rows),
        "force_2mm_over_closed": rows[-1]["force_N"] / rows[0]["force_N"],
        "pressure_at_1T_Pa": air_gap_maxwell_pressure(1.0),
        "force_at_1T_for_case_N": air_gap_maxwell_pressure(1.0) * POLE_AREA_M2 * ACTIVE_FACES,
    }

    assert checks["B_monotone_decreasing_with_gap"]
    assert checks["pressure_monotone_decreasing_with_gap"]
    assert checks["force_monotone_decreasing_with_gap"]
    assert checks["gap_mmf_fraction_monotone_increasing"]
    assert checks["max_pressure_identity_abs_error"] == 0.0
    assert checks["max_force_identity_abs_error"] == 0.0
    assert checks["force_2mm_over_closed"] < 0.01
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = build_rows()
    checks = validate(rows)
    summary = {
        "kind": "air_gap_maxwell_force_validation",
        "validation_class": True,
        "iron_path_m": IRON_PATH_M,
        "pole_area_m2": POLE_AREA_M2,
        "active_faces": ACTIVE_FACES,
        "mmf_A_turn": MMF_A_TURN,
        "gap_sweep_m": list(GAP_SWEEP_M),
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[air-gap Maxwell pressure / force sweep]")
    for row in rows:
        print(
            f"  g={1e3 * row['gap_m']:4.1f} mm  "
            f"B={row['B_T']:.6f} T  "
            f"p={row['pressure_Pa']:.3f} Pa  "
            f"F={row['force_N']:.6f} N  "
            f"gap_mmf={100.0 * row['gap_mmf_fraction']:6.2f}%"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
