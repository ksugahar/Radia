"""Validation-class nonlinear B-H winding-inductance sweep.

This is an example/validation run rather than a pytest test. It turns the
nonlinear series magnetic circuit into winding quantities:

* flux linkage and secant inductance at a driven current;
* incremental inductance from the B-H tangent;
* the constant-mu limit, where secant and incremental inductance collapse to the
  usual reluctance formula.

Run:

    python examples/electric_machine/validation_nonlinear_bh_inductance_sweep.py
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
    magnetic_circuit_bh_inductance_summary,
    magnetic_circuit_inductance,
)


OUT_JSON = HERE / "validation_nonlinear_bh_inductance_sweep_summary.json"

IRON_PATH_M = 0.0156
TURNS = 40.0
AREA_M2 = 2.5e-4
CURRENT_SWEEP_A = (0.2, 0.4, 0.8, 1.6, 3.2)
GAP_SWEEP_M = (0.0, 0.1e-3, 0.3e-3, 1.0e-3, 2.0e-3)
GAP_SWEEP_CURRENT_A = 2.4


def h_of_b(B: float) -> float:
    """Smooth monotone saturating iron curve, H [A/m] from B [T]."""
    return 51.0 * B + 2.5 * B ** 15


def dh_db(B: float) -> float:
    """Analytic derivative of :func:`h_of_b`."""
    return 51.0 + 37.5 * B ** 14


def linear_h_of_b(mu_r: float):
    return lambda B: B / (MU0 * mu_r)


def linear_dh_db(mu_r: float):
    return lambda B: 1.0 / (MU0 * mu_r)


def build_current_sweep() -> list[dict]:
    return [
        magnetic_circuit_bh_inductance_summary(
            TURNS,
            current,
            AREA_M2,
            IRON_PATH_M,
            h_of_b,
            gap=0.0,
            dh_db=dh_db,
        )
        for current in CURRENT_SWEEP_A
    ]


def build_gap_sweep() -> list[dict]:
    return [
        magnetic_circuit_bh_inductance_summary(
            TURNS,
            GAP_SWEEP_CURRENT_A,
            AREA_M2,
            IRON_PATH_M,
            h_of_b,
            gap=gap,
            dh_db=dh_db,
        )
        for gap in GAP_SWEEP_M
    ]


def build_constant_mu_rows() -> list[dict]:
    rows = []
    for mu_r, current, gap in (
        (300.0, 0.5, 0.0),
        (850.0, 1.25, 0.6e-3),
        (1500.0, 2.0, 1.2e-3),
    ):
        row = magnetic_circuit_bh_inductance_summary(
            TURNS,
            current,
            AREA_M2,
            IRON_PATH_M,
            linear_h_of_b(mu_r),
            gap=gap,
            dh_db=linear_dh_db(mu_r),
        )
        exact = magnetic_circuit_inductance(TURNS, AREA_M2, gap, IRON_PATH_M, mu_r)
        row.update({
            "mu_r_input": mu_r,
            "L_exact_linear_H": exact,
            "secant_L_exact_abs_error_H": abs(row["secant_inductance_H"] - exact),
            "incremental_L_exact_abs_error_H": abs(row["incremental_inductance_H"] - exact),
        })
        rows.append(row)
    return rows


def _strictly_increasing(values: list[float]) -> bool:
    return all(a < b for a, b in zip(values, values[1:]))


def _strictly_decreasing(values: list[float]) -> bool:
    return all(a > b for a, b in zip(values, values[1:]))


def validate(current_rows: list[dict], gap_rows: list[dict],
             constant_mu_rows: list[dict]) -> dict:
    current_B = [row["B_T"] for row in current_rows]
    current_Lsec = [row["secant_inductance_H"] for row in current_rows]
    current_Linc = [row["incremental_inductance_H"] for row in current_rows]
    gap_B = [row["B_T"] for row in gap_rows]
    gap_Lsec = [row["secant_inductance_H"] for row in gap_rows]
    gap_fraction = [row["gap_mmf_fraction"] for row in gap_rows[1:]]

    flux_identity_errors = [
        abs(row["flux_linkage_Wb_turn"] - row["turns"] * row["flux_Wb"])
        for row in current_rows + gap_rows + constant_mu_rows
    ]
    checks = {
        "current_B_monotone_increasing": _strictly_increasing(current_B),
        "current_secant_L_monotone_decreasing": _strictly_decreasing(current_Lsec),
        "current_incremental_L_monotone_decreasing": _strictly_decreasing(current_Linc),
        "current_incremental_L_below_secant": all(
            row["incremental_inductance_H"] < row["secant_inductance_H"]
            for row in current_rows
        ),
        "high_current_incremental_over_secant": current_rows[-1]["incremental_over_secant"],
        "gap_B_monotone_decreasing": _strictly_decreasing(gap_B),
        "gap_secant_L_monotone_decreasing": _strictly_decreasing(gap_Lsec),
        "gap_mmf_fraction_monotone_increasing": _strictly_increasing(gap_fraction),
        "gap_2mm_fraction": gap_rows[-1]["gap_mmf_fraction"],
        "max_flux_linkage_identity_error": max(flux_identity_errors),
        "constant_mu_max_secant_L_exact_abs_error_H": max(
            row["secant_L_exact_abs_error_H"] for row in constant_mu_rows
        ),
        "constant_mu_max_incremental_L_exact_abs_error_H": max(
            row["incremental_L_exact_abs_error_H"] for row in constant_mu_rows
        ),
    }

    assert checks["current_B_monotone_increasing"]
    assert checks["current_secant_L_monotone_decreasing"]
    assert checks["current_incremental_L_monotone_decreasing"]
    assert checks["current_incremental_L_below_secant"]
    assert checks["high_current_incremental_over_secant"] < 0.2
    assert checks["gap_B_monotone_decreasing"]
    assert checks["gap_secant_L_monotone_decreasing"]
    assert checks["gap_mmf_fraction_monotone_increasing"]
    assert checks["gap_2mm_fraction"] > 0.95
    assert checks["max_flux_linkage_identity_error"] < 1.0e-18
    assert checks["constant_mu_max_secant_L_exact_abs_error_H"] < 1.0e-18
    assert checks["constant_mu_max_incremental_L_exact_abs_error_H"] < 1.0e-17
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    current_rows = build_current_sweep()
    gap_rows = build_gap_sweep()
    constant_mu_rows = build_constant_mu_rows()
    checks = validate(current_rows, gap_rows, constant_mu_rows)

    summary = {
        "kind": "nonlinear_bh_winding_inductance_validation",
        "validation_class": True,
        "turns": TURNS,
        "area_m2": AREA_M2,
        "iron_path_m": IRON_PATH_M,
        "current_sweep_A": list(CURRENT_SWEEP_A),
        "gap_sweep_m": list(GAP_SWEEP_M),
        "gap_sweep_current_A": GAP_SWEEP_CURRENT_A,
        "current_sweep": current_rows,
        "gap_sweep": gap_rows,
        "constant_mu_rows": constant_mu_rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[nonlinear B-H winding-inductance current sweep]")
    for row in current_rows:
        print(
            f"  I={row['current_A']:4.1f} A  "
            f"B={row['B_T']:.6f} T  "
            f"Lsec={1e3 * row['secant_inductance_H']:.6f} mH  "
            f"Linc={1e3 * row['incremental_inductance_H']:.6f} mH  "
            f"Linc/Lsec={row['incremental_over_secant']:.4f}"
        )

    print("[fixed-current gap sweep]")
    for row in gap_rows:
        print(
            f"  g={1e3 * row['gap_m']:4.1f} mm  "
            f"B={row['B_T']:.6f} T  "
            f"Lsec={1e3 * row['secant_inductance_H']:.6f} mH  "
            f"gap_mmf={100.0 * row['gap_mmf_fraction']:6.2f}%"
        )

    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
