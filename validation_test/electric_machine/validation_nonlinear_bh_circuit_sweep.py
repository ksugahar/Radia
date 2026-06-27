"""Validation-class nonlinear B-H magnetic-circuit sweep.

This is an example/validation run rather than a pytest test. It turns a compact
series magnetic circuit into readable operating-point tables:

* a closed iron path moves into saturation as MMF rises;
* a small air gap lowers B but moves the iron back toward high incremental mu;
* the constant-mu limit collapses exactly to the linear reluctance formula.

Run:

    python validation_test/electric_machine/validation_nonlinear_bh_circuit_sweep.py
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
    magnetic_circuit_bh_operating_summary,
)


OUT_JSON = HERE / "validation_nonlinear_bh_circuit_sweep_summary.json"

IRON_PATH_M = 0.0156
MMF_SWEEP_A_TURN = (8.0, 16.0, 32.0, 64.0, 128.0, 256.0)
GAP_SWEEP_M = (0.0, 0.1e-3, 0.2e-3, 0.5e-3, 1.0e-3, 2.0e-3)
GAP_SWEEP_MMF_A_TURN = 96.0


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


def build_closed_path_sweep() -> list[dict]:
    return [
        magnetic_circuit_bh_operating_summary(
            mmf,
            IRON_PATH_M,
            h_of_b,
            gap=0.0,
            dh_db=dh_db,
        )
        for mmf in MMF_SWEEP_A_TURN
    ]


def build_gap_sweep() -> list[dict]:
    return [
        magnetic_circuit_bh_operating_summary(
            GAP_SWEEP_MMF_A_TURN,
            IRON_PATH_M,
            h_of_b,
            gap=gap,
            dh_db=dh_db,
        )
        for gap in GAP_SWEEP_M
    ]


def build_constant_mu_rows() -> list[dict]:
    rows = []
    for mu_r, mmf, gap in (
        (300.0, 40.0, 0.0),
        (750.0, 80.0, 0.8e-3),
        (1500.0, 120.0, 1.5e-3),
    ):
        row = magnetic_circuit_bh_operating_summary(
            mmf,
            IRON_PATH_M,
            linear_h_of_b(mu_r),
            gap=gap,
            dh_db=linear_dh_db(mu_r),
        )
        exact = mmf / (IRON_PATH_M / (MU0 * mu_r) + gap / MU0)
        row.update({
            "mu_r_input": mu_r,
            "B_exact_linear_T": exact,
            "B_exact_abs_error_T": abs(row["B_T"] - exact),
        })
        rows.append(row)
    return rows


def _strictly_increasing(values: list[float]) -> bool:
    return all(a < b for a, b in zip(values, values[1:]))


def _strictly_decreasing(values: list[float]) -> bool:
    return all(a > b for a, b in zip(values, values[1:]))


def validate(closed_rows: list[dict], gap_rows: list[dict],
             constant_mu_rows: list[dict]) -> dict:
    closed_B = [row["B_T"] for row in closed_rows]
    closed_mu_inc = [row["incremental_mu_r"] for row in closed_rows]
    closed_mu_sec = [row["apparent_mu_r"] for row in closed_rows]
    gap_B = [row["B_T"] for row in gap_rows]
    gap_mu_inc = [row["incremental_mu_r"] for row in gap_rows]
    gap_fraction = [row["gap_mmf_fraction"] for row in gap_rows[1:]]

    checks = {
        "closed_B_monotone_increasing": _strictly_increasing(closed_B),
        "closed_incremental_mu_monotone_decreasing": _strictly_decreasing(closed_mu_inc),
        "closed_incremental_mu_below_secant": all(
            inc < sec for inc, sec in zip(closed_mu_inc, closed_mu_sec)
        ),
        "closed_max_abs_mmf_residual_A_turn": max(
            abs(row["mmf_residual_A_turn"]) for row in closed_rows
        ),
        "closed_256A_B_over_32A_B": closed_rows[-1]["B_T"] / closed_rows[2]["B_T"],
        "gap_B_monotone_decreasing": _strictly_decreasing(gap_B),
        "gap_incremental_mu_monotone_increasing": _strictly_increasing(gap_mu_inc),
        "gap_mmf_fraction_monotone_increasing": _strictly_increasing(gap_fraction),
        "gap_2mm_fraction": gap_rows[-1]["gap_mmf_fraction"],
        "constant_mu_max_B_exact_abs_error_T": max(
            row["B_exact_abs_error_T"] for row in constant_mu_rows
        ),
        "constant_mu_max_mu_error": max(
            max(
                abs(row["apparent_mu_r"] - row["mu_r_input"]),
                abs(row["incremental_mu_r"] - row["mu_r_input"]),
            )
            for row in constant_mu_rows
        ),
    }

    assert checks["closed_B_monotone_increasing"]
    assert checks["closed_incremental_mu_monotone_decreasing"]
    assert checks["closed_incremental_mu_below_secant"]
    assert checks["closed_max_abs_mmf_residual_A_turn"] < 1.0e-9
    assert checks["closed_256A_B_over_32A_B"] < 1.4
    assert checks["gap_B_monotone_decreasing"]
    assert checks["gap_incremental_mu_monotone_increasing"]
    assert checks["gap_mmf_fraction_monotone_increasing"]
    assert checks["gap_2mm_fraction"] > 0.95
    assert checks["constant_mu_max_B_exact_abs_error_T"] < 1.0e-12
    assert checks["constant_mu_max_mu_error"] < 1.0e-9
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    closed_rows = build_closed_path_sweep()
    gap_rows = build_gap_sweep()
    constant_mu_rows = build_constant_mu_rows()
    checks = validate(closed_rows, gap_rows, constant_mu_rows)

    summary = {
        "kind": "nonlinear_bh_magnetic_circuit_validation",
        "validation_class": True,
        "iron_path_m": IRON_PATH_M,
        "mmf_sweep_A_turn": list(MMF_SWEEP_A_TURN),
        "gap_sweep_m": list(GAP_SWEEP_M),
        "gap_sweep_mmf_A_turn": GAP_SWEEP_MMF_A_TURN,
        "closed_path_sweep": closed_rows,
        "gap_sweep": gap_rows,
        "constant_mu_rows": constant_mu_rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[nonlinear B-H closed-path sweep]")
    for row in closed_rows:
        print(
            f"  NI={row['mmf_A_turn']:6.1f} A-turn  "
            f"B={row['B_T']:.6f} T  "
            f"H={row['H_iron_A_per_m']:9.3f} A/m  "
            f"mu_sec={row['apparent_mu_r']:8.1f}  "
            f"mu_inc={row['incremental_mu_r']:8.1f}"
        )

    print("[fixed-MMF gap sweep]")
    for row in gap_rows:
        print(
            f"  g={1e3 * row['gap_m']:4.1f} mm  "
            f"B={row['B_T']:.6f} T  "
            f"gap_mmf={100.0 * row['gap_mmf_fraction']:6.2f}%  "
            f"mu_inc={row['incremental_mu_r']:8.1f}"
        )

    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
